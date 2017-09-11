/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010-2014 Joel Andersson, Joris Gillis, Moritz Diehl,
 *                            K.U. Leuven. All rights reserved.
 *    Copyright (C) 2011-2014 Greg Horn
 *
 *    CasADi is free software; you can redistribute it and/or
 *    modify it under the terms of the GNU Lesser General Public
 *    License as published by the Free Software Foundation; either
 *    version 3 of the License, or (at your option) any later version.
 *
 *    CasADi is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *    Lesser General Public License for more details.
 *
 *    You should have received a copy of the GNU Lesser General Public
 *    License along with CasADi; if not, write to the Free Software
 *    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 *
 */

#include "optistack.hpp"
#include "nlpsol.hpp"

using namespace std;
namespace casadi {

class InternalOptiCallback : public Callback {
  public:
  InternalOptiCallback(OptiStack* sol);

  ~InternalOptiCallback();

  // Initialize the object
  void init() override;

  // Number of inputs and outputs
  int get_n_in() override;

  Sparsity get_sparsity_in(int i) override;

  // Evaluate numerically
  std::vector<DM> eval(const std::vector<DM>& arg) override;

  void set_sol(OptiStack* sol) { i=0; sol_= sol;}

  private:
    OptiStack* sol_;
    int i;
};


InternalOptiCallback::InternalOptiCallback(OptiStack* sol) : sol_(sol) {
}

InternalOptiCallback::~InternalOptiCallback() {  }

// Initialize the object
void InternalOptiCallback::init() {
}

// Number of inputs and outputs
int InternalOptiCallback::get_n_in() { return nlpsol_n_out();}

Sparsity InternalOptiCallback::get_sparsity_in(int i) {
  std::string n = nlpsol_out(i);
  int size = 0;
  if (n=="f") {
    size = 1;
  } else if (n=="lam_x" || n=="x") {
    size = sol_->nx();
  } else if (n=="lam_g" || n=="g") {
    size = sol_->ng();
  } else if (n=="p" || n=="lam_p") {
    size = sol_->np();
    if (size==0) return Sparsity::dense(0, 0);
  } else {
    return Sparsity::dense(0, 0);
  }
  return Sparsity::dense(size, 1);
}

// Evaluate numerically
std::vector<DM> InternalOptiCallback::eval(const std::vector<DM>& arg) {
  DMDict r;

  for (int i=0;i<nlpsol_n_out();++i) {
    r[nlpsol_out(i)] = arg[i];
  }

  sol_->res(r);


  for (OptiCallback* cb : sol_->callbacks_)
    cb->call(i);
  i+=1;
  return {0};
}

void OptiStack::callback_class(OptiCallback* callback) {
  callbacks_ = {callback};
}

void OptiStack::callback_class() {
  callbacks_ = {};
}

OptiStack::OptiStack() : count_(0), count_var_(0), count_par_(0) {
  solver_has_callback_ = false;
  internal_callback_ = 0;
  f_ = 0;
  instance_number_ = instance_count_++;
  mark_problem_dirty();
}

MX OptiStack::variable(int n, int m, const std::string& attribute) {

  // Prepare metadata
  MetaVar meta_data;
  meta_data.attribute = attribute;
  meta_data.n = n;
  meta_data.m = m;
  meta_data.type = OPTI_VAR;
  meta_data.count = count_++;
  meta_data.i = count_var_++;

  MX symbol, ret;

  if (attribute=="symmetric") {
    casadi_assert_message(n==m, "You specified attribute 'symmetric', "
      "while matrix is not even square, but " + std::to_string(n) + "-by-" + std::to_string(m) + ".");
    symbol = MX::sym(name_prefix() + "x_" + std::to_string(count_var_), n*(n+1)/2);
    ret = tril2symm(MX(Sparsity::lower(n), symbol));
  } else if (attribute=="full") {
    symbol = MX::sym(name_prefix() + "x_" + std::to_string(count_var_), n, m);
    ret = symbol;
  } else {
    casadi_error("Unknown attribute '" + attribute + "'. Choose from 'full' or 'symmetric'.");
  }

  // Store the symbol; preventing it from going ut of scope
  symbols_.push_back(symbol);
  initial_.push_back(DM::zeros(symbol.sparsity()));
  latest_.push_back(DM::nan(symbol.sparsity()));

  set_meta(symbol, meta_data);
  return ret;
}

std::string OptiStack::name_prefix() const {
  return "opti" + std::to_string(instance_number_) + "_";
}

void OptiStack::register_dual(MetaCon& c) {

  // Prepare metadata
  MetaVar meta_data;
  meta_data.attribute = "full";
  meta_data.n = c.canon.size1();
  meta_data.m = c.canon.size2();;
  meta_data.type = OPTI_DUAL_G;
  meta_data.count = count_++;
  meta_data.i = count_dual_++;

  MX symbol, ret;
  if (c.type==OPTI_PSD) {
    symbol = MX();
    ret = MX();
  } else {
    symbol = MX::sym(name_prefix()+"lam_g_"+std::to_string(count_dual_), c.canon.sparsity());

    casadi_assert(c.canon.is_dense());

    Sparsity ret_sp = repmat(c.original.sparsity(), 1, c.n);

    int N = c.canon.sparsity().nnz();

    MX flat = vec(symbol);
    if (c.type==OptiStack::OPTI_DOUBLE_INEQUALITY) {
      MX v = MX::sym("v");
      MX decide_left_right = vertcat(if_else_zero(v<0, -v), if_else_zero(v>=0, v));
      Function sign = Function("sign", {v}, {decide_left_right});
      Function sign_map = sign.map(c.canon.sparsity().nnz());
      ret = MX(ret_sp, sign_map((c.flipped ? -1 : 1)*flat)[0].T());
    } else {
      int block_size = N / c.n;
      std::vector<MX> original_blocks = vertsplit(fabs(flat), block_size);
      std::vector<MX> blocks(N);
      for (int i=0;i<c.n;++i) {
        int p = c.flipped? c.n-i-1: i;
        blocks[p] = original_blocks[i];
      }
      ret = MX(ret_sp, vertcat(blocks));
    }
  }

  symbols_.push_back(symbol);
  initial_duals_.push_back(DM::zeros(symbol.sparsity()));
  latest_duals_.push_back(DM::nan(symbol.sparsity()));

  c.dual = ret;
  c.dual_canon = symbol;

  set_meta(symbol, meta_data);
}

MX OptiStack::parameter(int n, int m, const std::string& attribute) {
  casadi_assert(attribute=="full");

  // Prepare metadata
  MetaVar meta_data;
  meta_data.attribute = attribute;
  meta_data.n = n;
  meta_data.m = m;
  meta_data.type = OPTI_PAR;
  meta_data.count = count_++;
  meta_data.i = count_par_++;

  MX symbol = MX::sym(name_prefix() + "p_" + std::to_string(count_par_), n, m);
  symbols_.push_back(symbol);
  values_.push_back(DM::nan(symbol.sparsity()));

  set_meta(symbol, meta_data);
  return symbol;
}

Dict OptiStack::stats() const {
  assert_solved();
  return solver_.stats();
}

std::string OptiStack::return_status() const {
  Dict mystats;
  try {
    mystats = stats();
  } catch (...) {
    //
  }
  if (mystats.find("return_status")!=mystats.end())
    return mystats.at("return_status");
  return "unknown";
}

Function OptiStack::casadi_solver() const {
  return solver_;
}

void OptiStack::set_meta(const MX& m, const MetaVar& meta) {
  meta_[m.get()] = meta;
}

void OptiStack::set_meta_con(const MX& m, const MetaCon& meta) {
  meta_con_[m.get()] = meta;
}

MX OptiStack::dual(const MX& m) const {
  return meta_con(m).dual;
}

const OptiStack::MetaVar& OptiStack::meta(const MX& m) const {
  assert_has(m);
  auto find = meta_.find(m.get());
  return find->second;
}

const OptiStack::MetaCon& OptiStack::meta_con(const MX& m) const {
  assert_has_con(m);
  auto find = meta_con_.find(m.get());
  return find->second;
}

OptiStack::MetaVar& OptiStack::meta(const MX& m) {
  assert_has(m);
  auto find = meta_.find(m.get());
  return find->second;
}

OptiStack::MetaCon& OptiStack::meta_con(const MX& m) {
  assert_has_con(m);
  auto find = meta_con_.find(m.get());
  return find->second;
}

OptiStack::MetaVar OptiStack::get_meta(const MX& m) const {
  return meta(m);
}

OptiStack::MetaCon OptiStack::get_meta_con(const MX& m) const {
  return meta_con(m);
}

bool OptiStack::has(const MX& m) const {
  return meta_.find(m.get())!=meta_.end();
}

bool OptiStack::has_con(const MX& m) const {
  return meta_con_.find(m.get())!=meta_con_.end();
}

void OptiStack::assert_has(const MX& m) const {
  if (!has(m)) {
    VariableType vt;
    if (parse_opti_name(m.name(), vt)) {
      casadi_error("Symbol '" + m.name() + "' (a " + variable_type_to_string(vt) + ") "
        "belongs to a different instance of Opti"
        "(this instance is #" + std::to_string(instance_count_) + ").");
    } else {
      casadi_error("MX symbol '" + m.name() + "' not declared with opti.variable/opti.parameter.\n"
        "Note: you cannot use a raw MX.sym in your Opti problem, "
        " only if you package it in a CasADi Function.");
    }
  }
}

void OptiStack::assert_has_con(const MX& m) const {
  casadi_assert_message(has_con(m), "Constraint not present in Opti stack.");
}

int OptiStack::instance_count_ = 0;

bool OptiStack::parse_opti_name(const std::string& name, VariableType& vt) const {
  int i = name.find("opti");
  if (i!=0) return false;

  i = name.find("_");
  i++;
  if (i==std::string::npos) return false;
  if (name.substr(i, 1)=="x") {
    vt = OPTI_VAR;
    return true;
  } else if (name.substr(i, 1)=="p") {
    vt = OPTI_PAR;
    return true;
  } else if (name.substr(i, 5)=="lam_g") {
    vt = OPTI_DUAL_G;
    return true;
  }

  return false;
}

std::string OptiStack::variable_type_to_string(VariableType vt) const {
  auto it = VariableType2String_.find(vt);
  if (it==VariableType2String_.end()) return "unknown variable type";
  return it->second;

}
std::map<OptiStack::VariableType, std::string> OptiStack::VariableType2String_ =
  {{OPTI_VAR, "decision variable"}, {OPTI_PAR, "parameter"}, {OPTI_DUAL_G, "dual variable"}};

std::vector<MX> OptiStack::initial() const {
  std::vector<MX> ret;
  for (const auto& e : symvar()) {
    if (meta(e).type==OPTI_VAR)
      ret.push_back(e==initial_[meta(e).i]);
    if (meta(e).type==OPTI_DUAL_G)
      ret.push_back(e==initial_duals_[meta(e).i]);
  }
  return ret;
}

std::vector<MX> OptiStack::value_variables() const {
  std::vector<MX> ret;
  for (const auto& e : symvar()) {
    if (meta(e).type==OPTI_VAR)
      ret.push_back(e==latest_[meta(e).i]);
  }
  return ret;
}

std::vector<MX> OptiStack::value_parameters() const {
  std::vector<MX> ret;
  for (const auto& e : symvar()) {
    if (meta(e).type==OPTI_PAR)
      ret.push_back(e==values_[meta(e).i]);
  }
  return ret;
}

void OptiStack::internal_bake() {
  casadi_assert_message(!f_.is_empty() || !g_.empty(),
    "You need to specify at least an objective (y calling 'minimize'), "
    "or a constraint (by calling 'subject_to').");

  symbol_active_.clear();
  symbol_active_.resize(symbols_.size());

  // Gather all expressions
  MX total_expr = vertcat(f_, veccat(g_));

  // Categorize the symbols appearing in those expressions
  for (const auto& d : symvar(total_expr))
    symbol_active_[meta(d).count] = true;

  std::vector<MX> x = active_symvar(OptiStack::OPTI_VAR);
  int offset = 0;
  for (const auto& v : x) {
    meta(v).start = offset;
    offset+= v.nnz();
    meta(v).stop = offset;
  }
  std::vector<MX> p = active_symvar(OptiStack::OPTI_PAR);

  // Fill the nlp definition
  nlp_["x"] = veccat(x);
  nlp_["p"] = veccat(p);

  nlp_["f"] = f_;

  offset = 0;
  for (int i=0;i<g_.size();++i) {
    OptiStack::MetaCon& r = meta_con(g_[i]);
    OptiStack::MetaVar& r2 = meta(r.dual_canon);
    symbol_active_[r2.count] = true;

    // Compute offsets for this constraint:
    // location into the global constraint variable
    r.start = offset;
    offset+= r.canon.nnz();
    r.stop = offset;

    r2.start = r.start;
    r2.stop  = r.stop;

  }

  lam_ = veccat(active_symvar(OptiStack::OPTI_DUAL_G));

  // Collect bounds and canonical form of constraints
  std::vector<MX> g_all;
  std::vector<MX> lbg_all;
  std::vector<MX> ubg_all;
  for (const auto& g : g_) {
    g_all.push_back(meta_con(g).canon);
    lbg_all.push_back(meta_con(g).lb);
    ubg_all.push_back(meta_con(g).ub);
  }

  nlp_["g"] = veccat(g_all);

  // Create bounds helper function
  MXDict bounds;
  bounds["p"] = nlp_["p"];
  bounds["lbg"] = veccat(lbg_all);
  bounds["ubg"] = veccat(ubg_all);

  bounds_ = Function("bounds", bounds, {"p"}, {"lbg", "ubg"});
  mark_problem_dirty(false);
}

void OptiStack::solver(const std::string& solver_name, const Dict& solver_options) {
  solver_name_ = solver_name;
  solver_options_ = solver_options;
  mark_solver_dirty();
}

std::vector<MX> OptiStack::sort(const std::vector<MX>& v) const {
  // We exploit the fact that std::map is ordered

  // Populate map
  std::map<int, MX> unordered;
  for (const auto& d : v)
    unordered[meta(d).count] = d;

  // Loop over map (ordered)
  std::vector<MX> ret;
  for (auto const &e : unordered)
    ret.push_back(e.second);
  return ret;
}

std::vector<MX> OptiStack::symvar() const {
  return symbols_;
}

std::vector<MX> OptiStack::symvar(const MX& expr) const {
  return sort(MX::symvar(expr));
}

std::vector<MX> OptiStack::ineq_unchain(const MX& a, bool& flipped) {
  flipped = false;
  casadi_assert(a.is_op(OP_LE) || a.is_op(OP_LT));

  // Is there inequalities in the left or right leaf?
  bool left  = a.dep(0).is_op(OP_LE) || a.dep(0).is_op(OP_LT);
  bool right = a.dep(1).is_op(OP_LE) || a.dep(1).is_op(OP_LT);
  casadi_assert(!left || !right);

  if (!left && !right)
    return {a.dep(0), a.dep(1)}; // Simple inequality

  // We have a chain of inequalities
  bool ineq = !left;
  std::vector<MX> ret = {a.dep(!ineq)};
  MX e = a.dep(ineq);
  while (e.is_op(OP_LE) || e.is_op(OP_LT)) {
    casadi_assert(!e.is_op(OP_EQ));
    casadi_assert(!e.dep(!ineq).is_op(OP_LE) && !e.dep(!ineq).is_op(OP_LT));
    ret.push_back(e.dep(!ineq));
    e = e.dep(ineq);
  }
  ret.push_back(e);
  if (left) std::reverse(ret.begin(), ret.end());
  flipped = !left;

  return ret;
}

void OptiStack::assert_only_opti_symbols(const MX& e) const {
  std::vector<MX> symbols = MX::symvar(e);
  for (const auto& s : symbols) assert_has(s);
}

void OptiStack::assert_only_opti_nondual(const MX& e) const {
  std::vector<MX> symbols = MX::symvar(e);
  for (const auto& s : symbols) {
    assert_has(s);
    casadi_assert_message(meta(s).type!=OPTI_DUAL_G, "Dual variables forbidden in this context.")
  }
}

bool OptiStack::is_parametric(const MX& expr) const {
  return symvar(expr, OptiStack::OPTI_VAR).empty();
}

OptiStack::MetaCon OptiStack::canon_expr(const MX& expr) const {
  MX c = expr;

  MetaCon con;
  con.original = expr;

  if (c.is_op(OP_LE) || c.is_op(OP_LT)) { // Inequalities
    std::vector<MX> ret;
    bool flipped;
    std::vector<MX> args = ineq_unchain(c, flipped);
    std::vector<bool> parametric;
    for (auto &a : args) parametric.push_back(is_parametric(a));

    if (args.size()==2 && (parametric[0] || parametric[1])) {
      // case: g(x,p) <= bound(p)
      MX e = args[0]-args[1];
      if (e.is_vector()) {
        casadi_assert(!parametric[0] || !parametric[1]);
        con.type = OPTI_INEQUALITY;
        if (parametric[0]) {
          con.lb = args[0]*DM::ones(e.sparsity());
          con.ub = inf*DM::ones(e.sparsity());
          con.canon = args[1];
        } else {
          con.lb = -inf*DM::ones(e.sparsity());
          con.ub = args[1]*DM::ones(e.sparsity());
          con.canon = args[0];
        }
        return con;
      }
      // Fall through to generic inequalities
    } else if (args.size()==3 && (parametric[0] || parametric[2])) {
      // lb(p) <= g(x,p) <= ub(p)
      con.type = OPTI_DOUBLE_INEQUALITY;
      con.lb = args[0]*DM::ones(args[1].sparsity());
      con.ub = args[2]*DM::ones(args[1].sparsity());
      con.canon = args[1];
      con.flipped = flipped;
      con.n = 2;
      return con;
    }

    bool type_known = false;
    for (int j=0;j<args.size()-1;++j) {
      MX e = args[j]-args[j+1];
      if (e.is_vector()) {
        // g1(x,p) <= g2(x,p)
        ret.push_back(e);
        casadi_assert(!type_known || con.type==OPTI_GENERIC_INEQUALITY);
        type_known = true;
        con.type = OPTI_GENERIC_INEQUALITY;
        con.flipped = flipped;
      } else {
        // A(x,p) >= b(p)
        e = args[j+1]-args[j];
        casadi_assert_message(e.size1()==e.size2(),
          "Matrix inequalities must be square. Did you mean element-wise inequality instead?");

        ret.push_back(e);
        casadi_assert(!type_known || con.type==OPTI_PSD);
        type_known = true;
        con.type = OPTI_PSD;
      }
    }

    if (con.type==OPTI_GENERIC_INEQUALITY) {
      con.canon = veccat(ret);
      con.lb = -inf*DM::ones(con.canon.sparsity());
      con.ub = DM::zeros(con.canon.sparsity());
      con.n = ret.size();
    } else {
      con.canon = diagcat(ret);
      con.n = ret.size();
    }
    return con;
  } else if (c.is_op(OP_EQ)) { // Inequalities

    MX e = c.dep(0)-c.dep(1);
    if (is_parametric(c.dep(0))) {
      con.canon = c.dep(1);
      con.lb = c.dep(0)*DM::ones(e.sparsity());
      con.type = OPTI_EQUALITY;
    } else if (is_parametric(c.dep(1))) {
      con.canon = c.dep(0);
      con.lb = c.dep(1)*DM::ones(e.sparsity());
      con.type = OPTI_EQUALITY;
    } else {
      con.lb = DM::zeros(e.sparsity());
      con.canon = e;
      con.type = OPTI_GENERIC_EQUALITY;
    }
    con.ub = con.lb;
    return con;
  } else { // Something else
    con.type = OPTI_UNKNOWN;
    con.canon = c;
    return con;
  }

}

void OptiStack::assert_solved() const {
  casadi_assert_message(solved(),
    "This action is forbidden since you have not solved the Opti stack yet "
    "(with calling 'solve').");
}

void OptiStack::assert_baked() const {
  casadi_assert_message(!problem_dirty(),
    "This action is forbidden since you have not baked the Opti stack yet "
    "(with calling 'solve').");
}

void OptiStack::assert_empty() const {
  casadi_assert(g_.empty());
  casadi_assert(f_.is_empty());
}

void OptiStack::minimize(const MX& f) {
  assert_only_opti_nondual(f);
  mark_problem_dirty();
  f_ = f;
}

void OptiStack::subject_to(const std::vector<MX>& g) {
  for (const auto& gs : g) subject_to(gs);
}

void OptiStack::subject_to(const MX& g) {
  assert_only_opti_nondual(g);
  mark_problem_dirty();
  g_.push_back(g);

  // Store the meta-data
  set_meta_con(g, canon_expr(g));

  register_dual(meta_con(g));
}

void OptiStack::subject_to() {
  mark_problem_dirty();
  g_.clear();
  latest_duals_.clear();
  initial_duals_.clear();
  count_dual_ = 0;
}

std::vector<MX> OptiStack::symvar(const MX& expr, VariableType type) const {
  std::vector<MX> ret;
  for (const auto& d : symvar(expr)) {
    if (meta(d).type==type) ret.push_back(d);
  }

  return ret;
}

void OptiStack::res(const DMDict& res) {
  const std::vector<double> & x_v = res.at("x").nonzeros();
  for (const auto &v : active_symvar(OPTI_VAR)) {
    int i = meta(v).i;
    std::vector<double> & data_v = latest_[i].nonzeros();
    std::copy(x_v.begin()+meta(v).start, x_v.begin()+meta(v).stop, data_v.begin());
  }
  if (res.find("lam_g")!=res.end()) {
    const std::vector<double> & lam_v = res.at("lam_g").nonzeros();
    for (const auto &v : active_symvar(OPTI_DUAL_G)) {
      int i = meta(v).i;
      std::vector<double> & data_v = latest_duals_[i].nonzeros();
      std::copy(lam_v.begin()+meta(v).start, lam_v.begin()+meta(v).stop, data_v.begin());
    }
  }
  res_ = res;
  mark_solved();
}

// Solve the problem
OptiSol OptiStack::solve() {

  if (problem_dirty()) {
    internal_bake();
  }
  // Verify the constraint types
  for (const auto& g : g_) {
    if (meta_con(g).type==OptiStack::OPTI_PSD)
      casadi_error("Psd constraints not implemented yet. "
      "Perhaps you intended an element-wise inequality? "
      "In that case, make sure that the matrix is flattened (e.g. mat(:)).");
  }

  bool need_callback = !callbacks_.empty();

  bool solver_update =  solver_dirty() ||  need_callback!=  solver_has_callback_;

  if (solver_update) {
    Dict opts = solver_options_;

    // Handle callbacks
    if (need_callback) {
      if (!internal_callback_) delete internal_callback_;
      internal_callback_ = new InternalOptiCallback(this);
      internal_callback_->construct("InternalOptiCallback", Dict());
      callback_ = *internal_callback_;
      internal_callback_->transfer_ownership();

      //callback_ = InternalOptiCallback::create();

      opts["iteration_callback"] = callback_;
    } else {
      if (!internal_callback_) delete internal_callback_;
      internal_callback_ = 0;
    }
    solver_has_callback_ = need_callback;

    casadi_assert_message(solver_name_!="",
      "You must call 'solver' on the Opti stack to select a solver.");
    solver_ = nlpsol("solver", solver_name_, nlp_, opts);
    mark_solver_dirty(false);
  }

  solve_prepare();
  res(solve_actual(arg_));

  std::string ret = return_status();

  bool success = ret=="Solve_Succeeded" || ret=="Solved_To_Acceptable_Level";

  casadi_assert_message(success,
    "Solver failed. You may use opti.debug.value to investigate the latest values of variables."
    " return_status is '" + ret + "'");

  return copy();
}

// Solve the problem
void OptiStack::solve_prepare() {


  // Verify the constraint types
  for (const auto& g : g_) {
    if (meta_con(g).type==OptiStack::OPTI_UNKNOWN)
     casadi_error("Constraint type unknown. Use ==, >= or <= .");
  }

  if (internal_callback_)
    internal_callback_->set_sol(this);

  // Get initial guess and parameter values
  arg_["x0"]     = veccat(active_values(OptiStack::OPTI_VAR));
  arg_["p"]      = veccat(active_values(OptiStack::OPTI_PAR));
  arg_["lam_g0"] = veccat(active_values(OptiStack::OPTI_DUAL_G));
  casadi_assert_message(arg_["p"].is_regular(),
    "You have forgotten to assign a value to a parameter ('set_value'), "
    "or have set it to NaN/Inf.");

  // Evaluate bounds for given parameter values
  DMDict arg;
  arg["p"] = arg_["p"];
  DMDict res = bounds_(arg);
  arg_["lbg"] = res["lbg"];
  arg_["ubg"] = res["ubg"];

}

DMDict OptiStack::solve_actual(const DMDict& arg) {
  return solver_(arg);
}

DM OptiStack::value(const MX& expr, const std::vector<MX>& values) const {
  std::vector<MX> x   = symvar(expr, OPTI_VAR);
  std::vector<MX> p   = symvar(expr, OPTI_PAR);
  std::vector<MX> lam = symvar(expr, OPTI_DUAL_G);

  Function helper = Function("helper", std::vector<MX>{veccat(x), veccat(p), veccat(lam)}, {expr});
  if (helper.has_free())
    casadi_error("This expression has symbols that are not defined "
      "within Opti using variable/parameter.");


  std::map<int, MX> temp;
  for (const auto& v : values) {
    casadi_assert(v.is_op(OP_EQ));
    int i = meta(v.dep(1)).i;
    casadi_assert(v.dep(0).is_constant());
    temp[i] = v.dep(0);
  }

  bool undecided_vars = false;
  std::vector<DM> x_num;
  for (const auto& e : x) {
    int i = meta(e).i;
    x_num.push_back(latest_[i]);

    // Override when values are supplied
    auto it = temp.find(i);
    if (it==temp.end()) {
      undecided_vars = true;
    } else {
      Slice all;
      DM t = static_cast<DM>(it->second);
      x_num.back().set(t, false, all, all);
    }

  }

  if (undecided_vars) {
    assert_solved();
    for (const auto& e : x)
      casadi_assert_message(symbol_active_[meta(e).count],
        "This expression has symbols that do not appear in the constraints and objective.");
  }

  std::vector<DM> p_num;
  for (const auto& e : p) {
    p_num.push_back(values_[meta(e).i]);
  }

  std::vector<DM> lam_num;
  if (lam.size()>0) {
    assert_solved();
    for (const auto& e : lam) {
      casadi_assert_message(symbol_active_[meta(e).count],
        "This expression has a dual for a constraint that is not given to Opti.");
      lam_num.push_back(latest_duals_[meta(e).i]);
    }
  }

  std::vector<DM> arg = helper(std::vector<DM>{veccat(x_num), veccat(p_num), veccat(lam_num)});
  return arg[0];
}

void OptiStack::assert_active_symbol(const MX& m) const {
  assert_has(m);
  assert_baked();
  casadi_assert_message(symbol_active_[meta(m).count], "Opti symbol is not used in Solver."
    " It does not make sense to assign a value to it.");
}

void OptiStack::set_initial(const std::vector<MX>& assignments) {
  for (const auto& v : assignments) {
    casadi_assert(v.is_op(OP_EQ));
    casadi_assert(v.dep(0).is_constant());
    if (has(v.dep(1)))
      set_initial(v.dep(1), static_cast<DM>(v.dep(0)));
  }
}

void OptiStack::set_value(const std::vector<MX>& assignments) {
  for (const auto& v : assignments) {
    casadi_assert(v.is_op(OP_EQ));
    casadi_assert(v.dep(0).is_constant());
    if (has(v.dep(1)))
      set_value(v.dep(1), static_cast<DM>(v.dep(0)));
  }
}

void OptiStack::set_value_internal(const MX& x, const DM& v, std::vector<DM>& store) {
  casadi_assert(v.is_regular());
  if (x.is_symbolic()) {
    DM& target = store[meta(x).i];
    Slice all;
    mark_solver_dirty();
    target.set(v, false, all, all);
    return;
  }

  // Obtain symbolic primitives
  std::vector<MX> symbols = MX::symvar(x);
  MX symbols_cat = veccat(symbols);

  std::string failmessage = "You cannot set initial/value of an arbitrary expression. "
    "Use symbols or simple mappings of symbols.";

  // Assert x is linear in its symbolic primitives
  for (bool b : which_depends(x, symbols_cat, 2, false)) casadi_assert_message(!b, failmessage);

  // Evaluate jacobian of expr wrt symbols
  Dict opts = {{"compact", true}};
  Function Jf("Jf", std::vector<MX>{}, std::vector<MX>{jacobian(x, veccat(symbols), opts)});
  DM J = Jf(std::vector<DM>{})[0];
  Sparsity sp_JT = J.T().sparsity();

  Function Ff("Ff", symbols, {x});
  DM E = Ff(std::vector<DM>(symbols.size(), 0))[0];
  std::vector<double>& e = E.nonzeros();

  // Cast the v input into the expected sparsity
  Slice all;
  DM value(x.sparsity());
  value.set(v, false, all, all);

  // Purge empty rows
  std::vector<int> filled_rows = sum2(J).get_row();
  J = J(filled_rows, all);

  // Get rows and columns of the mapping
  std::vector<int> row, col;
  J.sparsity().get_triplet(row, col);
  const std::vector<double>& scaling = J.nonzeros();
  const std::vector<double>& data_original = value.nonzeros();

  std::vector<double> data; data.reserve(value.nnz());
  for (int i=0;i<value.nnz();++i) {
    double v = data_original[i];
    int nz = sp_JT.colind()[i+1]-sp_JT.colind()[i];
    casadi_assert_message(nz<=1, failmessage);
    if (nz) {
      data.push_back(v);
    } else {
      casadi_assert_message(v==e[i], "In initial/value assignment: inconsistent numerical values. "
        "At nonzero " + std::to_string(i) + ", lhs has " + std::to_string(e[i]) + ", while rhs has " + std::to_string(v) + ".");
    }
  }

  // Contiguous workspace for nonzeros of all involved symbols
  std::vector<double> temp(symbols_cat.nnz(), casadi::nan);
  for (int k=0;k<data.size();++k) {
    double& lhs = temp[col[k]];
    double rhs = data[row[k]]/scaling[row[k]];
    if (std::isnan(lhs)) {
      // Assign in the workspace
      lhs = rhs;
    } else {
      casadi_assert_message(lhs==rhs, "Initial/value assignment with mapping is ambiguous.");
    }
  }

  int offset = 0;
  for (const auto & s : symbols) {
    DM& target = store[meta(s).i];
    std::vector<double>& data = target.nonzeros();
    // Loop over nonzeros in each symbol
    for (int i=0;i<s.nnz();++i) {
      // Copy from the workspace (barring fields that were not set)
      double v = temp[offset+i];
      if (!std::isnan(v)) data[i] = v;
    }
    offset+=s.nnz();
  }

}

void OptiStack::set_initial(const MX& x, const DM& v) {
  for (const auto & s : MX::symvar(x)) {
    casadi_assert_message(meta(s).type!=OPTI_PAR,
      "You cannot set an initial value for a parameter. Did you mean 'set_value'?");
    set_value_internal(x, v, meta(s).type==OPTI_VAR ? initial_ : initial_duals_);
  }
}

void OptiStack::set_value(const MX& x, const DM& v) {
  for (const auto & s : MX::symvar(x))
    casadi_assert_message(meta(s).type!=OPTI_VAR,
      "You cannot set a value for a variable. Did you mean 'set_initial'?");
  set_value_internal(x, v, values_);
}

std::vector<MX> OptiStack::active_symvar(OptiStack::VariableType type) const {
  if (symbol_active_.empty()) return std::vector<MX>{};
  std::vector<MX> ret;
  for (const auto& s : symbols_) {
    if (symbol_active_[meta(s).count] && meta(s).type==type)
      ret.push_back(s);
  }
  return ret;
}

std::vector<DM> OptiStack::active_values(OptiStack::VariableType type) const {
  if (symbol_active_.empty()) return std::vector<DM>{};
  std::vector<DM> ret;
  for (const auto& s : symbols_) {
    if (symbol_active_[meta(s).count] && meta(s).type==type) {
      if (type==OPTI_VAR) {
        ret.push_back(initial_[meta(s).i]);
      } else if (type==OPTI_PAR) {
        ret.push_back(values_[meta(s).i]);
      } else if (type==OPTI_DUAL_G) {
        ret.push_back(initial_duals_[meta(s).i]);
      }
    }
  }
  return ret;
}

void OptiStack::disp(ostream &stream, bool more) const {
  stream << "Opti {" << std::endl;
  stream << "  instance #" << instance_number_ << std::endl;
  OptiStack mycopy = copy();
  if (problem_dirty()) mycopy.internal_bake();
  stream << "  #variables: " << mycopy.active_symvar(OPTI_VAR).size()
    << " (nx = " << mycopy.nx() << ")" <<  std::endl;
  stream << "  #parameters: " << mycopy.active_symvar(OPTI_PAR).size()
    << " (np = " << mycopy.np() << ")" << std::endl;
  stream << "  #constraints: " << mycopy.active_symvar(OPTI_DUAL_G).size()
    << " (ng = " << mycopy.ng() << ")" << std::endl;
  if (solver_dirty()) {
    stream << "  CasADi solver needs updating." << std::endl;
  } else {
    stream << "  CasADi solver allocated." << std::endl;
  }
  if (solved()) {
    stream << "  CasADi solver was called: " << return_status() << std::endl;
  }
  stream << "}";
}

} // namespace casadi
