/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010-2023 Joel Andersson, Joris Gillis, Moritz Diehl,
 *                            KU Leuven. All rights reserved.
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

#include "optistack_internal.hpp"
#include "nlpsol.hpp"
#include "conic.hpp"
#include "function_internal.hpp"
#include "global_options.hpp"

namespace casadi {

class InternalOptiCallback : public FunctionInternal {
  public:

  InternalOptiCallback(OptiNode& sol) : FunctionInternal(class_name()), sol_(sol) {}

  ~InternalOptiCallback() override {
    clear_mem();
  }

  /** \brief Get type name */
  std::string class_name() const override {return "InternalOptiCallback";}

  // Number of inputs and outputs
  size_t get_n_in() override { return nlpsol_n_out();}

  Sparsity get_sparsity_in(casadi_int i) override {
    std::string n = nlpsol_out(i);
    casadi_int size = 0;
    if (n=="f") {
      size = 1;
    } else if (n=="lam_x" || n=="x") {
      size = sol_.nx();
    } else if (n=="lam_g" || n=="g") {
      size = sol_.ng();
    } else if (n=="p" || n=="lam_p") {
      size = sol_.np();
      if (size==0) return Sparsity::dense(0, 0);
    } else {
      return Sparsity::dense(0, 0);
    }
    return Sparsity::dense(size, 1);
  }

  void reset() { i=0; }

  // eval_dm has been defined instead of eval
  bool has_eval_dm() const override { return true;}

  /// Evaluate the function numerically
  std::vector<DM> eval_dm(const std::vector<DM>& arg) const override {
    DMDict r;

    for (casadi_int i=0;i<nlpsol_n_out();++i) {
      r[nlpsol_out(i)] = arg[i];
    }

    sol_.res(r);

    if (sol_.user_callback_) sol_.user_callback_->call(i);

    i+=1;
    return {0};
  }

  bool associated_with(const OptiNode* o) { return &sol_==o; }

  private:
    OptiNode& sol_;
    mutable casadi_int i;
};

OptiNode* OptiNode::create(const std::string& problem_type) {
return new OptiNode(problem_type);
}


void OptiNode::callback_class(OptiCallback* callback) {
  user_callback_ = callback;
}

void OptiNode::callback_class() {
  user_callback_ = nullptr;
}

bool OptiNode::has_callback_class() const {
  return user_callback_ != 0;
}

std::string OptiNode::format_stacktrace(const Dict& stacktrace, casadi_int indent) {
  std::string s_indent;
  for (casadi_int i=0;i<indent;++i) {
    s_indent+= "  ";
  }
  std::string description;
  std::string filename = stacktrace.at("file").as_string();
  casadi_int line = stacktrace.at("line").as_int();
  description += "called from " + filename +":"+str(line);
  std::string name = stacktrace.at("name").as_string();
  if (name!="Unknown" && name!= "<module>")
    description += " in " + stacktrace.at("name").as_string();
  try {
    std::ifstream file(filename);
    for (casadi_int i=0; i<line-1; ++i) {
      file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    }
    std::string contents; std::getline(file, contents);
    auto it = contents.find_first_not_of(" \n");
    if (it!=std::string::npos) {
      description += "\n" + s_indent + contents.substr(it);
    }
  } catch(...) {
    // pass
  }
  return description;
}

std::string OptiNode::describe(const MX& expr, casadi_int indent, const Dict& opts) const {
  if (problem_dirty()) return baked_copy().describe(expr, indent, opts);
  casadi_int max_stacktrace_depth = 1;
  for (const auto& op : opts) {
    if (op.first=="max_stacktrace_depth") {
      max_stacktrace_depth = op.second.as_int();
    } else {
      casadi_warning("Unknown option '" + op.first + "'");
    }
  }
  std::string s_indent;
  for (casadi_int i=0;i<indent;++i) {
    s_indent+= "  ";
  }
  std::string description = s_indent;
  if (expr.is_symbolic()) {
    if (has(expr)) {
      description += "Opti " + variable_type_to_string(meta(expr).type) + " '" + expr.name() +
        "' of shape " + expr.dim();
      const Dict& extra = meta(expr).extra;
      auto it = extra.find("stacktrace");
      if (it!=extra.end()) {
        for (const Dict& stacktrace : it->second.as_dict_vector()) {
          description += ", " + format_stacktrace(stacktrace, indent+1);
          if (--max_stacktrace_depth==0) break;
        }
      }
    } else {
      VariableType vt;
      if (parse_opti_name(expr.name(), vt)) {
        description += "Opti " + variable_type_to_string(vt) + " '" + expr.name() +
          "' of shape " + expr.dim()+
          ", belonging to a different instance of Opti.";
      } else {
        description += "MX symbol '" + expr.name() + "' of shape " + expr.dim();
        description += ", declared outside of Opti.";
      }
    }
  } else {
    if (has_con(expr)) {
      description = "Opti constraint of shape " + expr.dim();
      const Dict& extra = meta_con(expr).extra;
      auto it = extra.find("stacktrace");
      if (it!=extra.end()) {
        for (const Dict& stacktrace : it->second.as_dict_vector()) {
          description += ", " + format_stacktrace(stacktrace, indent+1);
          if (--max_stacktrace_depth==0) break;
        }
      }
    } else {
      std::vector<MX> s = symvar(expr);
      if (s.empty()) {
        description+= "Constant epxression.";
      } else {
        description+= "General expression, dependent on " + str(s.size()) + " symbols:";
        for (casadi_int i=0;i<s.size();++i) {
          description+= "\n"+describe(s[i], indent+1);
          if (i>5) {
            description+= "\n...";
            break;
          }
        }
      }
    }
  }

  return description;
}

std::string OptiNode::g_describe(casadi_int i, const Dict& opts) const {
  if (problem_dirty()) return baked_copy().g_describe(i, opts);
  MX expr = g_lookup(i);
  casadi_int local_i = i-meta_con(expr).start + GlobalOptions::start_index;
  std::string description = describe(expr, 0, opts);
  if (expr.numel()>1)
    description += "\nAt nonzero " + str(local_i % expr.numel()) +
                   ", part " + str((casadi_int) local_i / expr.numel()) + ".";
  return description;
}

std::string OptiNode::x_describe(casadi_int i, const Dict& opts) const {
  if (problem_dirty()) return baked_copy().x_describe(i, opts);
  MX symbol = x_lookup(i);
  casadi_int local_i = i-meta(symbol).start + GlobalOptions::start_index;
  std::string description = describe(symbol, 0, opts);
  if (symbol.numel()>1)
    description += "\nAt nonzero " + str(local_i) + ".";
  return description;
}

MX OptiNode::x_lookup(casadi_int i) const {
  if (problem_dirty()) return baked_copy().x_lookup(i);
  casadi_assert_dev(i>=0);
  casadi_assert_dev(i<nx());
  std::vector<MX> x = active_symvar(OPTI_VAR);
  for (const auto& e : x) {
    const MetaVar& m = meta(e);
    if (i>=m.start && i<m.stop) return e;
  }
  casadi_error("Internal error");
  return MX();
}

MX OptiNode::g_lookup(casadi_int i) const {
  if (problem_dirty()) return baked_copy().g_lookup(i);
  casadi_assert_dev(i>=0);
  casadi_assert_dev(i<ng());
  for (const auto& e : g_) {
    const MetaCon& m = meta_con(e);
    if (i>=m.start && i<m.stop) return e;
  }
  casadi_error("Internal error");
  return MX();
}

OptiNode::OptiNode(const std::string& problem_type) :
    count_(0), count_var_(0), count_par_(0), count_dual_(0) {
  f_ = 0;
  f_linear_scale_ = 1;
  instance_number_ = instance_count_++;
  user_callback_ = nullptr;
  store_initial_[OPTI_VAR] = {};
  store_initial_[OPTI_PAR] = {};
  store_initial_[OPTI_DUAL_G] = {};
  store_latest_[OPTI_VAR] = {};
  store_latest_[OPTI_DUAL_G] = {};
  casadi_assert(problem_type=="nlp" || problem_type=="conic",
    "Specified problem type '" + problem_type + "'unknown. "
    "Choose 'nlp' (default) or 'conic'.");
  problem_type_ = problem_type;
  mark_problem_dirty();
}

OptiNode::~OptiNode() {
}

MX OptiNode::variable(casadi_int n, casadi_int m, const std::string& attribute) {

  // Prepare metadata
  MetaVar meta_data;
  meta_data.attribute = attribute;
  meta_data.n = n;
  meta_data.m = m;
  meta_data.type = OPTI_VAR;
  meta_data.count = count_++;
  meta_data.i = count_var_++;
  meta_data.domain = OPTI_DOMAIN_REAL;

  MX symbol, ret;

  if (attribute=="symmetric") {
    casadi_assert(n==m, "You specified attribute 'symmetric', "
      "while matrix is not even square, but " + str(n) + "-by-" + str(m) + ".");
    symbol = MX::sym(name_prefix() + "x_" + str(count_var_), n*(n+1)/2);
    ret = tril2symm(MX(Sparsity::lower(n), symbol));
  } else if (attribute=="full") {
    symbol = MX::sym(name_prefix() + "x_" + str(count_var_), n, m);
    ret = symbol;
  } else {
    casadi_error("Unknown attribute '" + attribute + "'. Choose from 'full' or 'symmetric'.");
  }

  // Store the symbol; preventing it from going ut of scope
  symbols_.push_back(symbol);
  store_initial_[OPTI_VAR].push_back(DM::zeros(symbol.sparsity()));
  store_latest_[OPTI_VAR].push_back(DM::nan(symbol.sparsity()));
  store_linear_scale_[OPTI_VAR].push_back(DM::ones(symbol.sparsity()));
  store_linear_scale_offset_[OPTI_VAR].push_back(DM::zeros(symbol.sparsity()));

  set_meta(symbol, meta_data);
  return ret;
}

MX OptiNode::variable(const MX& symbol, const std::string& attribute) {

  // Prepare metadata
  MetaVar meta_data;
  meta_data.attribute = attribute;
  meta_data.n = symbol.size1();
  meta_data.m = symbol.size2();
  meta_data.type = OPTI_VAR;
  meta_data.count = count_++;
  meta_data.i = count_var_++;

  // Store the symbol; preventing it from going ut of scope
  symbols_.push_back(symbol);
  store_initial_[OPTI_VAR].push_back(DM::zeros(symbol.sparsity()));
  store_latest_[OPTI_VAR].push_back(DM::nan(symbol.sparsity()));
  store_linear_scale_[OPTI_VAR].push_back(DM::ones(symbol.sparsity()));
  store_linear_scale_offset_[OPTI_VAR].push_back(DM::zeros(symbol.sparsity()));

  set_meta(symbol, meta_data);
  return symbol;
}

MX OptiNode::variable(const Sparsity& sp, const std::string& attribute) {

  // Prepare metadata
  MetaVar meta_data;
  meta_data.attribute = attribute;
  meta_data.n = sp.size1();
  meta_data.m = sp.size2();
  meta_data.type = OPTI_VAR;
  meta_data.count = count_++;
  meta_data.i = count_var_++;

  MX symbol = MX::sym(name_prefix() + "x_" + str(count_var_), sp);

  // Store the symbol; preventing it from going ut of scope
  symbols_.push_back(symbol);
  store_initial_[OPTI_VAR].push_back(DM::zeros(symbol.sparsity()));
  store_latest_[OPTI_VAR].push_back(DM::nan(symbol.sparsity()));
  store_linear_scale_[OPTI_VAR].push_back(DM::ones(symbol.sparsity()));
  store_linear_scale_offset_[OPTI_VAR].push_back(DM::zeros(symbol.sparsity()));

  set_meta(symbol, meta_data);
  return symbol;
}

std::string OptiNode::name_prefix() const {
  return "opti" + str(instance_number_) + "_";
}

Opti OptiNode::copy() const {
    return Opti::create(new OptiNode(*this));
}

void OptiNode::register_dual(MetaCon& c) {

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
    symbol = MX::sym(name_prefix()+"lam_g_"+str(count_dual_), c.canon.sparsity());

    casadi_assert_dev(c.canon.is_dense());

    Sparsity ret_sp = repmat(c.original.sparsity(), 1, c.n);

    casadi_int N = c.canon.sparsity().nnz();

    MX flat = vec(symbol);
    if (c.type==OPTI_DOUBLE_INEQUALITY) {
      MX v = MX::sym("v");
      MX decide_left_right = vertcat(if_else_zero(v<0, -v), if_else_zero(v>=0, v));
      Function sign = Function("sign", {v}, {decide_left_right});
      Function sign_map = sign.map(c.canon.sparsity().nnz());
      ret = MX(ret_sp, sign_map((c.flipped ? -1 : 1)*flat)[0].T());
    } else {
      casadi_int block_size = N / c.n;
      std::vector<MX> original_blocks = vertsplit(fabs(flat), block_size);
      std::vector<MX> blocks(N);
      for (casadi_int i=0;i<c.n;++i) {
        casadi_int p = c.flipped? c.n-i-1: i;
        blocks[p] = original_blocks[i];
      }
      ret = MX(ret_sp, vertcat(blocks));
    }
  }

  symbols_.push_back(symbol);
  store_initial_[OPTI_DUAL_G].push_back(DM::zeros(symbol.sparsity()));
  store_latest_[OPTI_DUAL_G].push_back(DM::nan(symbol.sparsity()));

  c.dual = ret;
  c.dual_canon = symbol;

  set_meta(symbol, meta_data);
}

MX OptiNode::parameter(const MX& symbol, const std::string& attribute) {
  casadi_assert_dev(attribute=="full");

  // Prepare metadata
  MetaVar meta_data;
  meta_data.attribute = attribute;
  meta_data.n = symbol.size1();
  meta_data.m = symbol.size2();
  meta_data.type = OPTI_PAR;
  meta_data.count = count_++;
  meta_data.i = count_par_++;

  symbols_.push_back(symbol);
  store_initial_[OPTI_PAR].push_back(DM::nan(symbol.sparsity()));

  set_meta(symbol, meta_data);
  return symbol;
}

MX OptiNode::parameter(const Sparsity& sp, const std::string& attribute) {
  casadi_assert_dev(attribute=="full");

  // Prepare metadata
  MetaVar meta_data;
  meta_data.attribute = attribute;
  meta_data.n = sp.size1();
  meta_data.m = sp.size2();
  meta_data.type = OPTI_PAR;
  meta_data.count = count_++;
  meta_data.i = count_par_++;

  MX symbol = MX::sym(name_prefix() + "p_" + str(count_par_), sp);
  symbols_.push_back(symbol);
  store_initial_[OPTI_PAR].push_back(DM::nan(symbol.sparsity()));

  set_meta(symbol, meta_data);
  return symbol;
}

MX OptiNode::parameter(casadi_int n, casadi_int m, const std::string& attribute) {
  casadi_assert_dev(attribute=="full");

  // Prepare metadata
  MetaVar meta_data;
  meta_data.attribute = attribute;
  meta_data.n = n;
  meta_data.m = m;
  meta_data.type = OPTI_PAR;
  meta_data.count = count_++;
  meta_data.i = count_par_++;

  MX symbol = MX::sym(name_prefix() + "p_" + str(count_par_), n, m);
  symbols_.push_back(symbol);
  store_initial_[OPTI_PAR].push_back(DM::nan(symbol.sparsity()));

  set_meta(symbol, meta_data);
  return symbol;
}

Dict OptiNode::stats() const {
  assert_solved();
  if (stats_.empty()) {
    stats_ = solver_.stats();
    is_simple_ = get_from_dict(stats_,
      "detect_simple_bounds_is_simple",
      std::vector<bool>(ng(), false));
    target_x_ = get_from_dict(stats_,
      "detect_simple_bounds_target_x",
      std::vector<casadi_int>{});
    casadi_assert_dev(is_simple_.size()==ng());

    g_index_reduce_g_.resize(ng());
    g_index_reduce_x_.resize(ng(), -1);
    g_index_unreduce_g_.resize(ng(), -1);

    casadi_int* target_x_ptr = target_x_.data();

    casadi_int k=0;
    for (casadi_int i=0;i<is_simple_.size();++i) {
      if (!is_simple_[i]) {
        g_index_reduce_g_[i] = k;
        g_index_unreduce_g_[k] = i;
        k++;
      } else {
        g_index_reduce_g_[i] = -1;
        g_index_reduce_x_[i] = *target_x_ptr;
        target_x_ptr++;
      }
    }

    reduced_ = any(is_simple_);
  }
  return stats_;
}

casadi_int OptiNode::g_index_reduce_g(casadi_int i) const {
  stats();
  return g_index_reduce_g_[i];
}
casadi_int OptiNode::g_index_unreduce_g(casadi_int i) const {
  stats();
  return g_index_unreduce_g_[i];
}
casadi_int OptiNode::g_index_reduce_x(casadi_int i) const {
  stats();
  return g_index_reduce_x_[i];
}

std::string OptiNode::return_status() const {
  Dict mystats;
  try {
    mystats = stats();
  } catch (...) {
    //
  }
  if (mystats.find("return_status")!=mystats.end()) {
    std::stringstream ss;
    ss << mystats.at("return_status");
    return ss.str();
  }
  return "unknown";
}

bool OptiNode::return_success(bool accept_limit) const {
  Dict mystats;
  try {
    mystats = stats();
  } catch (...) {
    //
  }
  bool success = false;
  if (mystats.find("success")!=mystats.end()) success = mystats.at("success");
  if (!accept_limit) return success;

  bool limited = false;
  if (mystats.find("unified_return_status")!=mystats.end())
    limited = mystats.at("unified_return_status")=="SOLVER_RET_LIMITED";
  return success || limited;
}

Function OptiNode::casadi_solver() const {
  return solver_;
}

void OptiNode::set_meta(const MX& m, const MetaVar& meta) {
  meta_[m.get()] = meta;
}

void OptiNode::set_meta_con(const MX& m, const MetaCon& meta) {
  meta_con_[m.get()] = meta;
}

void OptiNode::update_user_dict(const MX& m, const Dict& meta) {
  if (has_con(m)) {
    MetaCon m_update = get_meta_con(m);
    MetaVar m_update2 = get_meta(m_update.dual_canon);
    for (const auto & it : meta) {
        m_update.extra[it.first] = it.second;
        m_update2.extra[it.first] = it.second;
    }
    set_meta_con(m, m_update);
    set_meta(m_update.dual_canon, m_update2);
  } else {
    for (const auto & s : MX::symvar(m)) {
      MetaVar m_update = get_meta(s);
      for (const auto & it : meta)
          m_update.extra[it.first] = it.second;
      set_meta(s, m_update);
    }
  }
}

Dict OptiNode::user_dict(const MX& m) const {
  if (has_con(m)) {
    MetaCon meta = get_meta_con(m);
    return meta.extra;
  } else {
    MetaVar meta = get_meta(m);
    return meta.extra;
  }
}

MX OptiNode::dual(const MX& m) const {
  return meta_con(m).dual;
}

const MetaVar& OptiNode::meta(const MX& m) const {
  assert_has(m);
  auto find = meta_.find(m.get());
  return find->second;
}

const MetaCon& OptiNode::meta_con(const MX& m) const {
  assert_has_con(m);
  auto find = meta_con_.find(m.get());
  return find->second;
}

MetaVar& OptiNode::meta(const MX& m) {
  assert_has(m);
  auto find = meta_.find(m.get());
  return find->second;
}

MetaCon& OptiNode::meta_con(const MX& m) {
  assert_has_con(m);
  auto find = meta_con_.find(m.get());
  return find->second;
}

MetaVar OptiNode::get_meta(const MX& m) const {
  return meta(m);
}

MetaCon OptiNode::get_meta_con(const MX& m) const {
  return meta_con(m);
}

bool OptiNode::has(const MX& m) const {
  return meta_.find(m.get())!=meta_.end();
}

bool OptiNode::has_con(const MX& m) const {
  return meta_con_.find(m.get())!=meta_con_.end();
}

void OptiNode::assert_has(const MX& m) const {
  if (!has(m)) {
    VariableType vt;
    casadi_assert(m.is_symbolic(), "Symbol expected, got expression.");
    if (parse_opti_name(m.name(), vt)) {
      casadi_error("Unknown: " + describe(m));
    } else {
      casadi_error("Unknown: " + describe(m) + "\n"
        "Note: you cannot use a raw MX.sym in your Opti problem,"
        " only if you package it in a CasADi Function.");
    }
  }
}

void OptiNode::assert_has_con(const MX& m) const {
  casadi_assert(has_con(m), "Constraint not present in Opti stack.");
}

casadi_int OptiNode::instance_count_ = 0;

bool OptiNode::parse_opti_name(const std::string& name, VariableType& vt) const {
  casadi_int i = name.find("opti");
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

std::string OptiNode::variable_type_to_string(VariableType vt) const {
  auto it = VariableType2String_.find(vt);
  if (it==VariableType2String_.end()) return "unknown variable type";
  return it->second;

}
std::map<VariableType, std::string> OptiNode::VariableType2String_ =
  {{OPTI_VAR, "decision variable"},
   {OPTI_PAR, "parameter"},
   {OPTI_DUAL_G, "dual variable"}};

std::vector<MX> OptiNode::initial() const {
  std::vector<MX> ret;
  for (const auto& e : symvar()) {
    if (meta(e).type==OPTI_VAR || meta(e).type==OPTI_DUAL_G)
      ret.push_back(e==store_initial_.at(meta(e).type)[meta(e).i]);
  }
  return ret;
}

std::vector<MX> OptiNode::value_variables() const {
  std::vector<MX> ret;
  for (const auto& e : symvar()) {
    if (meta(e).type==OPTI_VAR)
      ret.push_back(e==store_latest_.at(meta(e).type)[meta(e).i]);
  }
  return ret;
}

std::vector<MX> OptiNode::value_parameters() const {
  std::vector<MX> ret;
  for (const auto& e : symvar()) {
    if (meta(e).type==OPTI_PAR)
      ret.push_back(e==store_initial_.at(meta(e).type)[meta(e).i]);
  }
  return ret;
}

void OptiNode::bake() {
  casadi_assert(!f_.is_empty() || !g_.empty(),
    "You need to specify at least an objective (y calling 'minimize'), "
    "or a constraint (by calling 'subject_to').");

  symbol_active_.clear();
  symbol_active_.resize(symbols_.size());

  // Gather all expressions
  MX total_expr = vertcat(f_, veccat(g_));

  // Categorize the symbols appearing in those expressions
  for (const auto& d : symvar(total_expr))
    symbol_active_[meta(d).count] = true;

  std::vector<MX> x = active_symvar(OPTI_VAR);
  for (casadi_int i=0;i<x.size();++i) meta(x[i]).active_i = i;

  casadi_int offset = 0;
  for (const auto& v : x) {
    meta(v).start = offset;
    offset+= v.nnz();
    meta(v).stop = offset;
  }
  std::vector<MX> p = active_symvar(OPTI_PAR);
  for (casadi_int i=0;i<p.size();++i) meta(p[i]).active_i = i;

  // Fill the nlp definition
  nlp_["x"] = veccat(x);
  nlp_["p"] = veccat(p);

  nlp_unscaled_["x"] = veccat(x);
  nlp_unscaled_["p"] = veccat(p);

  discrete_.clear();
  for (const MX& e : x) {
    discrete_.insert(discrete_.end(), e.nnz(), meta(e).domain == OPTI_DOMAIN_INTEGER);
  }

  nlp_["f"] = f_;
  nlp_unscaled_["f"] = f_;

  offset = 0;
  for (casadi_int i=0;i<g_.size();++i) {
    MetaCon& r = meta_con(g_[i]);
    MetaVar& r2 = meta(r.dual_canon);
    symbol_active_[r2.count] = true;

    // Compute offsets for this constraint:
    // location into the global constraint variable
    r.start = offset;
    offset+= r.canon.nnz();
    r.stop = offset;

    r2.start = r.start;
    r2.stop  = r.stop;

  }

  std::vector<MX> lam = active_symvar(OPTI_DUAL_G);
  for (casadi_int i=0;i<lam.size();++i) meta(lam[i]).active_i = i;

  lam_ = veccat(lam);

  // Collect bounds and canonical form of constraints
  std::vector<MX> g_all, g_unscaled_all;
  std::vector<MX> h_all, h_unscaled_all;
  std::vector<MX> lbg_all, lbg_unscaled_all;
  std::vector<MX> ubg_all, ubg_unscaled_all;

  g_linear_scale_.clear();
  h_linear_scale_.clear();
  std::vector<DM> g_linear_scale, h_linear_scale;

  equality_.clear();
  for (const auto& g : g_) {
    if (meta_con(g).type==OPTI_PSD) {
      h_all.push_back(meta_con(g).canon/meta_con(g).linear_scale);
      h_linear_scale.push_back(meta_con(g).linear_scale);
      h_unscaled_all.push_back(meta_con(g).canon);
    } else {
      g_all.push_back(meta_con(g).canon/meta_con(g).linear_scale);
      g_linear_scale.push_back(meta_con(g).linear_scale);
      g_unscaled_all.push_back(meta_con(g).canon);
      lbg_all.push_back(meta_con(g).lb/meta_con(g).linear_scale);
      lbg_unscaled_all.push_back(meta_con(g).lb);
      ubg_all.push_back(meta_con(g).ub/meta_con(g).linear_scale);
      ubg_unscaled_all.push_back(meta_con(g).ub);

      equality_.insert(equality_.end(),
        meta_con(g).canon.numel(),
        meta_con(g).type==OPTI_EQUALITY || meta_con(g).type==OPTI_GENERIC_EQUALITY);
    }
  }

  nlp_["g"] = veccat(g_all);
  g_linear_scale_ = veccat(g_linear_scale).nonzeros();
  nlp_unscaled_["g"] = veccat(g_unscaled_all);
  if (problem_type_=="conic") {
    nlp_["h"] = diagcat(h_all);
    nlp_unscaled_["h"] = diagcat(h_unscaled_all);
    h_linear_scale_ = veccat(h_linear_scale).nonzeros();
  }

  // Get scaling data
  std::vector<DM> linear_scale = active_values(OPTI_VAR, store_linear_scale_);
  std::vector<DM> linear_scale_offset = active_values(OPTI_VAR, store_linear_scale_offset_);

  linear_scale_ = veccat(linear_scale).nonzeros();
  linear_scale_offset_ = veccat(linear_scale_offset).nonzeros();

  // Unscaled version of x
  std::vector<MX> x_unscaled(x.size());
  for (casadi_int i=0;i<x.size();++i) {
    x_unscaled[i] = x[i]*linear_scale[i] + linear_scale_offset[i];
  }

  // Perform substitution
  std::vector<MX> expr = {nlp_["f"], nlp_["g"]};
  if (problem_type_=="conic") expr.push_back(nlp_["h"]);
  std::vector<MX> fgh = substitute(expr, x, x_unscaled);
  nlp_["f"] = fgh[0];
  nlp_["g"] = fgh[1];
  if (problem_type_=="conic") {
    nlp_["h"] = fgh[2];
  }

  // Scale of objective
  nlp_["f"] = nlp_["f"]/f_linear_scale_;

  // Create bounds helper function
  MXDict bounds;
  bounds["p"] = nlp_["p"];
  bounds_lbg_ = veccat(lbg_all);
  bounds_ubg_ = veccat(ubg_all);
  bounds_unscaled_lbg_ = veccat(lbg_unscaled_all);
  bounds_unscaled_ubg_ = veccat(ubg_unscaled_all);

  bounds["lbg"] = bounds_lbg_;
  bounds["ubg"] = bounds_ubg_;

  bounds_ = Function("bounds", bounds, {"p"}, {"lbg", "ubg"});
  mark_problem_dirty(false);
}

void OptiNode::solver(const std::string& solver_name, const Dict& plugin_options,
                       const Dict& solver_options) {
  solver_name_ = solver_name;
  solver_options_ = plugin_options;
  if (!solver_options.empty())
    solver_options_[solver_name] = solver_options;
  mark_solver_dirty();
}

std::vector<MX> OptiNode::sort(const std::vector<MX>& v) const {
  // We exploit the fact that std::map is ordered

  // Populate map
  std::map<casadi_int, MX> unordered;
  for (const auto& d : v)
    unordered[meta(d).count] = d;

  // Loop over map (ordered)
  std::vector<MX> ret;
  for (auto const &e : unordered)
    ret.push_back(e.second);
  return ret;
}

std::vector<MX> OptiNode::symvar() const {
  return symbols_;
}

std::vector<MX> OptiNode::symvar(const MX& expr) const {
  return sort(MX::symvar(expr));
}

std::vector<MX> OptiNode::ineq_unchain(const MX& a, bool& flipped) {
  flipped = false;
  casadi_assert_dev(a.is_op(OP_LE) || a.is_op(OP_LT));

  // Is there inequalities in the left or right leaf?
  bool left  = a.dep(0).is_op(OP_LE) || a.dep(0).is_op(OP_LT);
  bool right = a.dep(1).is_op(OP_LE) || a.dep(1).is_op(OP_LT);
  casadi_assert_dev(!left || !right);

  if (!left && !right)
    return {a.dep(0), a.dep(1)}; // Simple inequality

  // We have a chain of inequalities
  bool ineq = !left;
  std::vector<MX> ret = {a.dep(!ineq)};
  MX e = a.dep(ineq);
  while (e.is_op(OP_LE) || e.is_op(OP_LT)) {
    casadi_assert_dev(!e.is_op(OP_EQ));
    casadi_assert_dev(!e.dep(!ineq).is_op(OP_LE) && !e.dep(!ineq).is_op(OP_LT));
    ret.push_back(e.dep(!ineq));
    e = e.dep(ineq);
  }
  ret.push_back(e);
  if (left) std::reverse(ret.begin(), ret.end());
  flipped = !left;

  return ret;
}

void OptiNode::assert_only_opti_symbols(const MX& e) const {
  std::vector<MX> symbols = MX::symvar(e);
  for (const auto& s : symbols) assert_has(s);
}

void OptiNode::assert_only_opti_nondual(const MX& e) const {
  std::vector<MX> symbols = MX::symvar(e);
  for (const auto& s : symbols) {
    assert_has(s);
    casadi_assert(meta(s).type!=OPTI_DUAL_G, "Dual variables forbidden in this context.");
  }
}

bool OptiNode::is_parametric(const MX& expr) const {
  return symvar(expr, OPTI_VAR).empty();
}

MetaCon OptiNode::canon_expr(const MX& expr, const DM& linear_scale) const {
  MX c = expr;

  MetaCon con;
  con.original = expr;
  casadi_assert(linear_scale.is_scalar() || linear_scale.size()==expr.size(),
    "Linear scale must have the same size as the expression. "
    "You got linear_scale " + con.linear_scale.dim() + " while " + expr.dim() + " is expected.");
  con.linear_scale = linear_scale;

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
        casadi_assert(!parametric[0] || !parametric[1],
          "Constraint must contain decision variables.");
        if (problem_type_=="conic") {
          if (args[0].op()==OP_NORMF || args[0].op()==OP_NORM2) {
            args[0] = -soc(args[0].dep(), args[1]);
            args[1] = 0;
          }
        } else {
          con.type = OPTI_INEQUALITY;
          if (parametric[0]) {
            con.lb = args[0]*DM::ones(e.sparsity());
            con.ub = inf*DM::ones(e.sparsity());
            con.canon = args[1]*DM::ones(e.sparsity());
          } else {
            con.lb = -inf*DM::ones(e.sparsity());
            con.ub = args[1]*DM::ones(e.sparsity());
            con.canon = args[0]*DM::ones(e.sparsity());
          }
          return con;
        }
      }
      // Fall through to generic inequalities
    } else if (args.size()==3 && parametric[0] && parametric[2]) {
      // lb(p) <= g(x,p) <= ub(p)
      con.type = OPTI_DOUBLE_INEQUALITY;
      con.lb = args[0]*DM::ones(args[1].sparsity());
      con.ub = args[2]*DM::ones(args[1].sparsity());
      con.canon = args[1]*DM::ones(args[1].sparsity());
      con.flipped = flipped;
      con.n = 2;
      return con;
    }

    bool type_known = false;
    for (casadi_int j=0;j<args.size()-1;++j) {
      MX e = args[j]-args[j+1];
      if (problem_type_=="conic") {
        if (args[j].op()==OP_NORMF || args[j].op()==OP_NORM2) {
          args[j] = -soc(args[j].dep(), args[j+1]);
          args[j+1] = 0;
          e = args[j]-args[j+1];
        }
      }
      if (e.is_vector()) {
        // g1(x,p) <= g2(x,p)
        ret.push_back(e);
        casadi_assert_dev(!type_known || con.type==OPTI_GENERIC_INEQUALITY);
        type_known = true;
        con.type = OPTI_GENERIC_INEQUALITY;
        con.flipped = flipped;
      } else {
        // A(x,p) >= b(p)
        MX a = args[j+1];
        MX b = args[j];
        e = a-b;

        casadi_assert(e.size1()==e.size2(),
          "Matrix inequalities must be square. Did you mean element-wise inequality instead?");
        if (a.is_scalar()) a*= MX::eye(e.size1());
        if (b.is_scalar()) b*= MX::eye(e.size1());
        e = a-b;

        ret.push_back(e);
        casadi_assert_dev(!type_known || con.type==OPTI_PSD);
        type_known = true;
        con.type = OPTI_PSD;
      }
    }

    if (con.type==OPTI_GENERIC_INEQUALITY) {
      con.canon = veccat(ret);
      con.lb = -inf*DM::ones(con.canon.sparsity());
      con.ub = DM::zeros(con.canon.sparsity());
      con.n = ret.size();
      if (!con.linear_scale.is_scalar()) {
        con.linear_scale = repmat(con.linear_scale, ret.size(), 1);
      }
    } else {
      con.canon = diagcat(ret);
      con.n = ret.size();
    }
    return con;
  } else if (c.is_op(OP_EQ)) { // Inequalities
    casadi_assert(!is_parametric(c.dep(0)) || !is_parametric(c.dep(1)),
      "Constraint must contain decision variables.");
    MX e = c.dep(0)-c.dep(1);
    if (is_parametric(c.dep(0))) {
      con.canon = c.dep(1)*DM::ones(e.sparsity());
      con.lb = c.dep(0)*DM::ones(e.sparsity());
      con.type = OPTI_EQUALITY;
      casadi_assert(c.dep(0).size1()<=c.dep(1).size1() && c.dep(0).size2()<=c.dep(1).size2(),
        "Constraint shape mismatch.");
    } else if (is_parametric(c.dep(1))) {
      con.canon = c.dep(0)*DM::ones(e.sparsity());
      con.lb = c.dep(1)*DM::ones(e.sparsity());
      con.type = OPTI_EQUALITY;
      casadi_assert(c.dep(1).size1()<=c.dep(0).size1() && c.dep(1).size2()<=c.dep(0).size2(),
        "Constraint shape mismatch.");
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

void OptiNode::assert_solved() const {
  casadi_assert(solved(),
    "This action is forbidden since you have not solved the Opti stack yet "
    "(with calling 'solve').");
}

void OptiNode::assert_baked() const {
  casadi_assert(!problem_dirty(),
    "This action is forbidden since you have not baked the Opti stack yet "
    "(with calling 'solve').");
}

void OptiNode::assert_empty() const {
  casadi_assert_dev(g_.empty());
  casadi_assert_dev(f_.is_empty());
}

void OptiNode::minimize(const MX& f, double linear_scale) {
  assert_only_opti_nondual(f);
  mark_problem_dirty();
  casadi_assert(f.is_scalar(), "Objective must be scalar, got " + f.dim() + ".");
  f_ = f;
  f_linear_scale_ = linear_scale;
}

void OptiNode::subject_to(const MX& g, const DM& linear_scale, const Dict& options) {
  assert_only_opti_nondual(g);
  mark_problem_dirty();
  g_.push_back(g);

  casadi_assert(!g.is_empty(),    "You passed an empty expression to `subject_to`. "
                                  "Make sure the number of rows and columns is non-zero. "
                                  "Got " + g.dim(true) + ".");
  casadi_assert(g.nnz()>0,        "You passed a fully sparse expression to `subject_to`. "
                                  "Make sure the expression has at least one nonzero. "
                                  "Got " + g.dim(true) + ".");
  casadi_assert(!g.is_constant(), "You passed a constant to `subject_to`. "
                                  "You need a symbol to form a constraint.");

  // Store the meta-data
  set_meta_con(g, canon_expr(g, linear_scale));
  register_dual(meta_con(g));

  for (auto && it : options) {
    if (it.first=="stacktrace") {
      meta_con(g).extra["stacktrace"] = it.second.to_dict_vector();
      meta(meta_con(g).dual_canon).extra["stacktrace"] = it.second.to_dict_vector();
    } else if (it.first=="meta") {
      update_user_dict(g, it.second.to_dict());
    } else {
      casadi_error("Unknown option: " + it.first);
    }
  }
}

void OptiNode::subject_to() {
  mark_problem_dirty();
  g_.clear();
  store_initial_[OPTI_DUAL_G].clear();
  store_latest_[OPTI_DUAL_G].clear();
  count_dual_ = 0;
}

std::vector<MX> OptiNode::symvar(const MX& expr, VariableType type) const {
  std::vector<MX> ret;
  for (const auto& d : symvar(expr)) {
    if (meta(d).type==type) ret.push_back(d);
  }

  return ret;
}

void OptiNode::res(const DMDict& res) {
  const std::vector<double> & x_v = res.at("x").nonzeros();
  for (const auto &v : active_symvar(OPTI_VAR)) {
    casadi_int i = meta(v).i;
    std::vector<double> & data_v = store_latest_[OPTI_VAR][i].nonzeros();
    for (casadi_int i=0;i<data_v.size();++i) {
      casadi_int j = meta(v).start+i;
      data_v[i] = x_v[j]*linear_scale_[j] + linear_scale_offset_[j];
    }
  }
  if (res.find("lam_g")!=res.end() && problem_type_!="conic") {
    const std::vector<double> & lam_v = res.at("lam_g").nonzeros();
    for (const auto &v : active_symvar(OPTI_DUAL_G)) {
      casadi_int i = meta(v).i;
      std::vector<double> & data_v = store_latest_[OPTI_DUAL_G][i].nonzeros();
      for (casadi_int i=0;i<data_v.size();++i) {
        casadi_int j = meta(v).start+i;
        data_v[i] = lam_v[j]/g_linear_scale_[j]*f_linear_scale_;
      }
    }
  }
  res_ = res;
  mark_solved();
}

bool OptiNode::old_callback() const {
  if (callback_.is_null()) return false;
  InternalOptiCallback* cb = static_cast<InternalOptiCallback*>(callback_.get());
  return !cb->associated_with(this);
}

Function OptiNode::solver_construct(bool callback) {
  if (problem_dirty()) {
    bake();
  }

  // Verify the constraint types
  for (const auto& g : g_) {
    if (problem_type_!="conic") {
      if (meta_con(g).type==OPTI_PSD)
        casadi_error("Psd constraints not implemented yet. "
        "Perhaps you intended an element-wise inequality? "
        "In that case, make sure that the matrix is flattened (e.g. mat(:)).");
    }
  }

  Dict solver_options_all = solver_options_;

  if (solver_options_all.find("equality")==solver_options_all.end()) {
    solver_options_all["equality"] = equality_;
  }

  if (solver_options_all.find("discrete")==solver_options_all.end()) {
    solver_options_all["discrete"] = discrete_;
  }

  Dict opts = solver_options_all;

  // Handle callbacks
  if (callback && user_callback_) {
    callback_ = Function::create(new InternalOptiCallback(*this), Dict());
    opts["iteration_callback"] = callback_;
  }

  casadi_assert(!solver_name_.empty(),
    "You must call 'solver' on the Opti stack to select a solver. "
    "Suggestion: opti.solver('ipopt')");

  if (problem_type_=="conic") {
    return qpsol("solver", solver_name_, nlp_, opts);
  } else {
    return nlpsol("solver", solver_name_, nlp_, opts);
  }

}

// Solve the problem
OptiSol OptiNode::solve(bool accept_limit) {

  if (problem_dirty()) {
    bake();
  }

  bool solver_update =  solver_dirty() || old_callback() || (user_callback_ && callback_.is_null());

  if (solver_update) {
    solver_ = solver_construct(true);
    mark_solver_dirty(false);
  }

  solve_prepare();
  res(solve_actual(arg_));

  std::string ret = return_status();

  casadi_assert(return_success(accept_limit),
    "Solver failed. You may use opti.debug.value to investigate the latest values of variables."
    " return_status is '" + ret + "'");

  return copy();
}

// Solve the problem
void OptiNode::solve_prepare() {


  // Verify the constraint types
  for (const auto& g : g_) {
    if (meta_con(g).type==OPTI_UNKNOWN)
     casadi_error("Constraint type unknown. Use ==, >= or <= .");
  }

  if (user_callback_) {
    InternalOptiCallback* cb = static_cast<InternalOptiCallback*>(callback_.get());
    cb->reset();
  }

  // Get initial guess and parameter values
  arg_["x0"]     = (veccat(active_values(OPTI_VAR))-linear_scale_offset_)/linear_scale_;
  arg_["p"]      = veccat(active_values(OPTI_PAR));
  arg_["lam_g0"] = veccat(active_values(OPTI_DUAL_G));
  if (!arg_["p"].is_regular()) {
    std::vector<MX> s = active_symvar(OPTI_PAR);
    std::vector<DM> v = active_values(OPTI_PAR);
    for (casadi_int i=0;i<s.size();++i) {
      casadi_assert(v[i].is_regular(),
        "You have forgotten to assign a value to a parameter ('set_value'), "
        "or have set it to NaN/Inf:\n" + describe(s[i], 1));
    }
  }


  // Evaluate bounds for given parameter values
  DMDict arg;
  arg["p"] = arg_["p"];
  DMDict res = bounds_(arg);
  arg_["lbg"] = res["lbg"];
  arg_["ubg"] = res["ubg"];

  stats_.clear();
}

DMDict OptiNode::solve_actual(const DMDict& arg) {
  return solver_(arg);
}

bool override_num(const std::map<casadi_int, MX> & temp, std::vector<DM>& num, casadi_int i) {
  // Override when values are supplied
  auto it = temp.find(i);
  if (it==temp.end()) {
    return true;
  } else {
    Slice all;
    DM t = static_cast<DM>(it->second);
    num.back().set(t, false, all, all);
  }
  return false;
}

DM OptiNode::value(const MX& expr, const std::vector<MX>& values, bool scaled) const {
  std::vector<MX> x   = symvar(expr, OPTI_VAR);
  std::vector<MX> p   = symvar(expr, OPTI_PAR);
  std::vector<MX> lam = symvar(expr, OPTI_DUAL_G);

  Function helper = Function("helper", std::vector<MX>{veccat(x), veccat(p), veccat(lam)}, {expr});
  if (helper.has_free())
    casadi_error("This expression has symbols that are not defined "
      "within Opti using variable/parameter.");

  std::map<VariableType, std::map<casadi_int, MX> > temp;
  temp[OPTI_DUAL_G] = std::map<casadi_int, MX>();
  for (const auto& v : values) {
    casadi_assert_dev(v.is_op(OP_EQ));
    casadi_int i = meta(v.dep(1)).i;
    casadi_assert_dev(v.dep(0).is_constant());
    temp[meta(v.dep(1)).type][i] = v.dep(0);
  }

  bool undecided_vars = false;
  std::vector<DM> x_num;
  for (const auto& e : x) {
    casadi_int i = meta(e).i;
    x_num.push_back(store_latest_.at(OPTI_VAR).at(i));
    undecided_vars |= override_num(temp[OPTI_VAR], x_num, i);
    if (scaled) {
      x_num.back() = x_num.back()/store_linear_scale_.at(OPTI_VAR)[meta(e).i] -
                     store_linear_scale_offset_.at(OPTI_VAR)[meta(e).i];
    }
  }

  std::vector<DM> lam_num;
  for (const auto& e : lam) {
    casadi_int i = meta(e).i;
    casadi_assert(i<store_latest_.at(OPTI_DUAL_G).size(),
      "This expression has a dual for a constraint that is not given to Opti:\n" +
      describe(e, 1));
    lam_num.push_back(store_latest_.at(OPTI_DUAL_G).at(i));
    undecided_vars |= override_num(temp[OPTI_DUAL_G], lam_num, i);
  }

  std::vector<DM> p_num;
  for (const auto& e : p) {
    casadi_int i = meta(e).i;
    p_num.push_back(store_initial_.at(OPTI_PAR).at(i));
    override_num(temp[OPTI_PAR], p_num, i);
    casadi_assert(p_num.back().is_regular(),
      "This expression depends on a parameter with unset value:\n"+
      describe(e, 1));
  }

  if (undecided_vars) {
    assert_solved();
    for (const auto& e : x)
      casadi_assert(symbol_active_[meta(e).count],
        "This expression has symbols that do not appear in the constraints and objective:\n" +
        describe(e, 1));
    for (const auto& e : lam)
      casadi_assert(symbol_active_[meta(e).count],
        "This expression has a dual for a constraint that is not given to Opti:\n" +
        describe(e, 1));
  }

  std::vector<DM> arg = helper(std::vector<DM>{veccat(x_num), veccat(p_num), veccat(lam_num)});
  return arg[0];
}

void OptiNode::assert_active_symbol(const MX& m) const {
  assert_has(m);
  assert_baked();
  casadi_assert(symbol_active_[meta(m).count], "Opti symbol is not used in Solver."
    " It does not make sense to assign a value to it:\n" + describe(m, 1));
}

void OptiNode::set_initial(const std::vector<MX>& assignments) {
  for (const auto& v : assignments) {
    casadi_assert_dev(v.is_op(OP_EQ));
    casadi_assert_dev(v.dep(0).is_constant());
    if (has(v.dep(1)))
      set_initial(v.dep(1), static_cast<DM>(v.dep(0)));
  }
}

void OptiNode::set_domain(const MX& x, const std::string& domain) {
  mark_solved(false);
  mark_problem_dirty();
  casadi_assert(x.is_valid_input(), "First argument to set_domain should be a variable.");
  DomainType type;
  if (domain=="real") {
    type = OPTI_DOMAIN_REAL;
  } else if (domain=="integer") {
    type = OPTI_DOMAIN_INTEGER;
  } else {
    casadi_error("Unknown domain '" + domain + "'. Known values are 'real', 'integer'.");
  }
  for (const auto& prim : x.primitives()) {
    MetaVar& m = meta(prim);
    m.domain = type;
  }
}

void OptiNode::set_value(const std::vector<MX>& assignments) {
  for (const auto& v : assignments) {
    casadi_assert_dev(v.is_op(OP_EQ));
    casadi_assert_dev(v.dep(0).is_constant());
    if (has(v.dep(1)))
      set_value(v.dep(1), static_cast<DM>(v.dep(0)));
  }
}

void OptiNode::set_value_internal(const MX& x, const DM& v,
    std::map< VariableType, std::vector<DM> >& store) {
  mark_solved(false);
  casadi_assert_dev(v.is_regular());
  if (x.is_symbolic()) {
    DM& target = store[meta(x).type][meta(x).i];
    Slice all;
    target.set(v, false, all, all);
    return;
  }

  // Obtain symbolic primitives
  std::vector<MX> symbols = MX::symvar(x);
  MX symbols_cat = veccat(symbols);

  std::string failmessage = "You cannot set initial/value of an arbitrary expression. "
    "Use symbols or simple mappings of symbols.";

  // Assert x is linear in its symbolic primitives
  for (bool b : which_depends(x, symbols_cat, 2, false)) casadi_assert(!b, failmessage);

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
  std::vector<casadi_int> filled_rows = sum2(J).get_row();
  J = J(filled_rows, all); // NOLINT(cppcoreguidelines-slicing)

  // Get rows and columns of the mapping
  std::vector<casadi_int> row, col;
  J.sparsity().get_triplet(row, col);
  const std::vector<double>& scaling = J.nonzeros();
  const std::vector<double>& data_original = value.nonzeros();

  std::vector<double> data; data.reserve(value.nnz());
  for (casadi_int i=0;i<value.nnz();++i) {
    double v = data_original[i];
    casadi_int nz = sp_JT.colind()[i+1]-sp_JT.colind()[i];
    casadi_assert(nz<=1, failmessage);
    if (nz) {
      data.push_back(v);
    } else {
      casadi_assert(v==e[i], "In initial/value assignment: "
        "inconsistent numerical values. At nonzero " + str(i) + ", lhs has "
        + str(e[i]) + ", while rhs has " + str(v) + ".");
    }
  }

  // Contiguous workspace for nonzeros of all involved symbols
  std::vector<double> temp(symbols_cat.nnz(), casadi::nan);
  for (casadi_int k=0;k<data.size();++k) {
    double& lhs = temp[col[k]];
    double rhs = data[row[k]]/scaling[row[k]];
    if (std::isnan(lhs)) {
      // Assign in the workspace
      lhs = rhs;
    } else {
      casadi_assert(lhs==rhs, "Initial/value assignment with mapping is ambiguous.");
    }
  }

  casadi_int offset = 0;
  for (const auto & s : symbols) {
    DM& target = store[meta(s).type][meta(s).i];
    std::vector<double>& data = target.nonzeros();
    // Loop over nonzeros in each symbol
    for (casadi_int i=0;i<s.nnz();++i) {
      // Copy from the workspace (barring fields that were not set)
      double v = temp[offset+i];
      if (!std::isnan(v)) data[i] = v;
    }
    offset+=s.nnz();
  }

}

void OptiNode::set_initial(const MX& x, const DM& v) {
  for (const auto & s : MX::symvar(x))
    casadi_assert(meta(s).type!=OPTI_PAR,
      "You cannot set an initial value for a parameter. Did you mean 'set_value'?");
  set_value_internal(x, v, store_initial_);
}

void OptiNode::set_value(const MX& x, const DM& v) {
  for (const auto & s : MX::symvar(x))
    casadi_assert(meta(s).type!=OPTI_VAR,
      "You cannot set a value for a variable. Did you mean 'set_initial'?");
  set_value_internal(x, v, store_initial_);
}

void OptiNode::set_linear_scale(const MX& x, const DM& scale, const DM& offset) {
  for (const auto & s : MX::symvar(x))
    casadi_assert(meta(s).type!=OPTI_PAR,
      "You cannot set a scale value for a parameter.");
  casadi_assert(scale.is_scalar() || scale.size()==x.size(),
      "Dimension mismatch in linear_scale. Expected " + x.dim() + ", got " + scale.dim()+ ".");
  set_value_internal(x, scale, store_linear_scale_);
  casadi_assert(offset.is_scalar() || offset.size()==x.size(),
      "Dimension mismatch in linear_scale offset. Expected " + x.dim() +
      ", got " + scale.dim()+ ".");
  set_value_internal(x, offset, store_linear_scale_offset_);
}

std::vector<MX> OptiNode::active_symvar(VariableType type) const {
  if (symbol_active_.empty()) return std::vector<MX>{};
  std::vector<MX> ret;
  for (const auto& s : symbols_) {
    if (symbol_active_[meta(s).count] && meta(s).type==type)
      ret.push_back(s);
  }
  return ret;
}

std::vector<DM> OptiNode::active_values(VariableType type) const {
  return active_values(type, store_initial_);
}

std::vector<DM> OptiNode::active_values(VariableType type,
    const std::map< VariableType, std::vector<DM> >& store) const {
  if (symbol_active_.empty()) return std::vector<DM>{};
  std::vector<DM> ret;
  for (const auto& s : symbols_) {
    if (symbol_active_[meta(s).count] && meta(s).type==type) {
      ret.push_back(store.at(meta(s).type)[meta(s).i]);
    }
  }
  return ret;
}

Function OptiNode::to_function(const std::string& name,
    const std::vector<MX>& args, const std::vector<MX>& res,
    const std::vector<std::string>& name_in,
    const std::vector<std::string>& name_out,
    const Dict& opts) {
  if (problem_dirty()) return baked_copy().to_function(name, args, res, name_in, name_out, opts);

  Function solver = solver_construct(false);

  // Get initial guess and parameter values
  std::vector<MX> x0, p, lam_g;
  assign_vector(active_values(OPTI_VAR), x0);
  assign_vector(active_values(OPTI_PAR), p);
  assign_vector(active_values(OPTI_DUAL_G), lam_g);

  casadi_int k = 0;
  for (const auto& a : args) {
    casadi_assert(a.is_valid_input(), "Argument " + str(k) + " is not purely symbolic.");
    k++;
    for (const auto& prim : a.primitives()) {
      if (!symbol_active_[meta(prim).count]) continue;
      casadi_int i = meta(prim).active_i;
      if (meta(prim).type==OPTI_VAR) {
        x0.at(i) = prim;
      } else if (meta(prim).type==OPTI_PAR) {
        p.at(i) = prim;
      } else if (meta(prim).type==OPTI_DUAL_G) {
        lam_g.at(i) = prim;
      } else {
        casadi_error("Unknown");
      }
    }
  }
  MXDict arg;
  arg["p"] = veccat(p);

  // Evaluate bounds for given parameter values
  MXDict r = bounds_(arg);
  arg["x0"] = veccat(x0);
  arg["lam_g0"] = veccat(lam_g);
  arg["lbg"] = r["lbg"];
  arg["ubg"] = r["ubg"];

  r = solver(arg);

  std::vector<MX> helper_in = {veccat(active_symvar(OPTI_VAR)),
                               veccat(active_symvar(OPTI_PAR)),
                               veccat(active_symvar(OPTI_DUAL_G))};
  Function helper("helper", helper_in, {res});

  std::vector<MX> arg_in = helper(std::vector<MX>{r.at("x"), arg["p"], r.at("lam_g")});

  return Function(name, args, arg_in, name_in, name_out, opts);

}

void OptiNode::disp(std::ostream &stream, bool more) const {

}

casadi_int OptiNode::instance_number() const {
    return instance_number_;
}

void OptiNode::show_infeasibilities(double tol, const Dict& opts) const {
  stats();
  std::vector<double> lbg_ = value(lbg()).get_elements();
  std::vector<double> ubg_ = value(ubg()).get_elements();
  std::vector<double> g_scaled_ = value(nlp_.at("g"), std::vector<MX>(), true).get_elements();
  std::vector<double> lbg_scaled_ = value(bounds_lbg_, std::vector<MX>(), true).get_elements();
  std::vector<double> ubg_scaled_ = value(bounds_ubg_, std::vector<MX>(), true).get_elements();


  std::vector<double> g_ = value(g()).get_elements();
  uout() << "Violated constraints (tol " << tol << "), in order of declaration:" << std::endl;

  for (casadi_int i=0;i<g_.size();++i) {
    double err = std::max(g_[i]-ubg_[i], lbg_[i]-g_[i]);
    double err_scaled = std::max(g_scaled_[i]-ubg_scaled_[i], lbg_scaled_[i]-g_scaled_[i]);
    if (err>=tol) {
      uout() << "------- i = " << i+GlobalOptions::start_index;
      uout() << "/" << g_.size();

      if (reduced_) {
        if (is_simple_[i]) {
          if (GlobalOptions::start_index==0) {
            uout() << "  reduced to bound on x[" << g_index_reduce_x_.at(i) << "]";
          } else {
            uout() << "  reduced to bound on x(" << g_index_reduce_x_.at(i)+1 << ")";
          }
        } else {
          uout() << "  reduced to g[" << g_index_reduce_g_.at(i) << "]";
        }
      }

      uout() << " ------ " << std::endl;
      uout() << lbg_[i] << " <= " << g_[i] << " <= " << ubg_[i];
      uout() << " (viol " << err << ")" << std::endl;
      if (g_[i]!=g_scaled_[i]) {
        uout() << lbg_scaled_[i] << " <= " << g_scaled_[i] << " <= " << ubg_scaled_[i];
        uout() << " (scaled) (viol " << err_scaled << ")" << std::endl;
      }
      uout() << g_describe(i, opts) << std::endl;
    }
  }
}

} // namespace casadi
