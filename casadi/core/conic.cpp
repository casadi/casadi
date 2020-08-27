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


#include "conic_impl.hpp"
#include "nlpsol_impl.hpp"

using namespace std;
namespace casadi {

  bool has_conic(const string& name) {
    return Conic::has_plugin(name);
  }

  void load_conic(const string& name) {
    Conic::load_plugin(name);
  }

  string doc_conic(const string& name) {
    return Conic::getPlugin(name).doc;
  }

  Function conic(const string& name, const string& solver,
                const SpDict& qp, const Dict& opts) {
    return Function::create(Conic::instantiate(name, solver, qp), opts);
  }

  void conic_debug(const Function& f, const std::string &filename) {
    ofstream file;
    file.open(filename.c_str());
    conic_debug(f, file);
  }

  void conic_debug(const Function& f, std::ostream &file) {
    casadi_assert_dev(!f.is_null());
    const Conic* n = f.get<Conic>();
    return n->generateNativeCode(file);
  }

  vector<string> conic_in() {
    vector<string> ret(conic_n_in());
    for (size_t i=0; i<ret.size(); ++i) ret[i]=conic_in(i);
    return ret;
  }

  vector<string> conic_out() {
    vector<string> ret(conic_n_out());
    for (size_t i=0; i<ret.size(); ++i) ret[i]=conic_out(i);
    return ret;
  }

  string conic_in(casadi_int ind) {
    switch (static_cast<ConicInput>(ind)) {
    case CONIC_H:      return "h";
    case CONIC_G:      return "g";
    case CONIC_A:      return "a";
    case CONIC_Q:      return "q";
    case CONIC_P:      return "p";
    case CONIC_LBA:    return "lba";
    case CONIC_UBA:    return "uba";
    case CONIC_LBX:    return "lbx";
    case CONIC_UBX:    return "ubx";
    case CONIC_X0:     return "x0";
    case CONIC_LAM_X0: return "lam_x0";
    case CONIC_LAM_A0: return "lam_a0";
    case CONIC_NUM_IN: break;
    }
    return string();
  }

  string conic_out(casadi_int ind) {
    switch (static_cast<ConicOutput>(ind)) {
    case CONIC_X:     return "x";
    case CONIC_COST:  return "cost";
    case CONIC_LAM_A: return "lam_a";
    case CONIC_LAM_X: return "lam_x";
    case CONIC_NUM_OUT: break;
    }
    return string();
  }

  casadi_int conic_n_in() {
    return CONIC_NUM_IN;
  }

  casadi_int conic_n_out() {
    return CONIC_NUM_OUT;
  }

  template<typename M>
  Function qpsol_nlp(const std::string& name, const std::string& solver,
                     const std::map<std::string, M>& qp, const Dict& opts) {
    // We have: minimize    f(x) = 1/2 * x' H x + c'x
    //          subject to  lbx <= x <= ubx
    //                      lbg <= g(x) = A x + b <= ubg
    //                      h(x) >=0 (psd)

    // Extract 'expand' option
    Dict opt = opts;
    auto it = opt.find("expand");
    bool expand = false;
    if (it!=opt.end()) {
      expand = it->second;
      opt.erase(it);
    }
    if (expand && M::type_name()=="MX") {
      Function f = Function("f", qp, {"x", "p"}, {"f", "g", "h"});
      std::vector<SX> arg = f.sx_in();
      std::vector<SX> res = f(arg);
      SXDict qp_mod;
      for (casadi_int i=0;i<f.n_in();++i) qp_mod[f.name_in(i)] = arg[i];
      for (casadi_int i=0;i<f.n_out();++i) qp_mod[f.name_out(i)] = res[i];
      return qpsol_nlp(name, solver, qp_mod, opt);
    }

    M x, p, f, g, h;
    for (auto&& i : qp) {
      if (i.first=="x") {
        x = i.second;
      } else if (i.first=="p") {
        p = i.second;
      } else if (i.first=="f") {
        f = i.second;
      } else if (i.first=="g") {
        g = i.second;
      } else if (i.first=="h") {
        h = i.second;
      } else {
        casadi_error("No such field: " + i.first);
      }
    }

    if (f.is_empty()) f = 0;
    if (g.is_empty()) g = M(0, 1);

    // Dimension checks
    casadi_assert(g.is_dense() && g.is_vector(),
      "Expected a dense vector 'g', but got " + g.dim() + ".");

    casadi_assert(f.is_dense() && f.is_scalar(),
      "Expected a dense scalar 'f', but got " + f.dim() + ".");

    casadi_assert(x.is_dense() && x.is_vector(),
      "Expected a dense vector 'x', but got " + x.dim() + ".");

    casadi_assert(h.is_square(),
      "Expected a symmetric matrix 'h', but got " + h.dim() + ".");

    if (g.is_empty(true)) g = M(0, 1); // workaround

    // Gradient of the objective: gf == Hx + g
    M gf = M::gradient(f, x);

    // Identify the linear term in the objective
    M c = substitute(gf, x, M::zeros(x.sparsity()));

    // Identify the quadratic term in the objective
    M H = M::jacobian(gf, x, {{"symmetric", true}});

    // Identify constant term in the objective
    Function r("constant_qp", {x, p}, {substitute(f, x, M::zeros(x.sparsity()))});

    // Identify the constant term in the constraints
    M b = substitute(g, x, M::zeros(x.sparsity()));

    // Identify the linear term in the constraints
    M A = M::jacobian(g, x);

    // Identify the constant term in the psd constraints
    M P = substitute(h, x, M::zeros(x.sparsity()));

    // Identify the linear term in the psd constraints
    M Q = M::jacobian(h, x);

    // Create a function for calculating the required matrices vectors
    Function prob(name + "_qp", {x, p}, {H, c, A, b, Q, P});

    // Make sure that the problem is sound
    casadi_assert(!prob.has_free(), "Cannot create '" + prob.name() + "' "
                          "since " + str(prob.get_free()) + " are free.");

    // Create the QP solver
    Function conic_f = conic(name + "_qpsol", solver,
                             {{"h", H.sparsity()}, {"a", A.sparsity()},
                              {"p", P.sparsity()}, {"q", Q.sparsity()}}, opt);

    // Create an MXFunction with the right signature
    vector<MX> ret_in(NLPSOL_NUM_IN);
    ret_in[NLPSOL_X0] = MX::sym("x0", x.sparsity());
    ret_in[NLPSOL_P] = MX::sym("p", p.sparsity());
    ret_in[NLPSOL_LBX] = MX::sym("lbx", x.sparsity());
    ret_in[NLPSOL_UBX] = MX::sym("ubx", x.sparsity());
    ret_in[NLPSOL_LBG] = MX::sym("lbg", g.sparsity());
    ret_in[NLPSOL_UBG] = MX::sym("ubg", g.sparsity());
    ret_in[NLPSOL_LAM_X0] = MX::sym("lam_x0", x.sparsity());
    ret_in[NLPSOL_LAM_G0] = MX::sym("lam_g0", g.sparsity());
    vector<MX> ret_out(NLPSOL_NUM_OUT);


    vector<MX> v(NL_NUM_IN);
    v[NL_X] = ret_in[NLPSOL_X0];
    v[NL_P] = ret_in[NLPSOL_P];
    // Evaluate constant part of objective
    MX rv = r(v)[0];
    // Get expressions for the QP matrices and vectors
    v = prob(v);

    // Call the QP solver
    vector<MX> w(CONIC_NUM_IN);
    w[CONIC_H] = v.at(0);
    w[CONIC_G] = v.at(1);
    w[CONIC_A] = v.at(2);
    w[CONIC_Q] = v.at(4);
    w[CONIC_P] = v.at(5);
    w[CONIC_LBX] = ret_in[NLPSOL_LBX];
    w[CONIC_UBX] = ret_in[NLPSOL_UBX];
    w[CONIC_LBA] = ret_in[NLPSOL_LBG] - v.at(3);
    w[CONIC_UBA] = ret_in[NLPSOL_UBG] - v.at(3);
    w[CONIC_X0] = ret_in[NLPSOL_X0];
    w[CONIC_LAM_X0] = ret_in[NLPSOL_LAM_X0];
    w[CONIC_LAM_A0] = ret_in[NLPSOL_LAM_G0];
    w = conic_f(w);

    // Get expressions for the solution
    ret_out[NLPSOL_X] = reshape(w[CONIC_X], x.size());
    ret_out[NLPSOL_F] = rv + w[CONIC_COST];
    ret_out[NLPSOL_G] = reshape(mtimes(v.at(2), w[CONIC_X]), g.size()) + v.at(3);
    ret_out[NLPSOL_LAM_X] = reshape(w[CONIC_LAM_X], x.size());
    ret_out[NLPSOL_LAM_G] = reshape(w[CONIC_LAM_A], g.size());
    ret_out[NLPSOL_LAM_P] = MX::nan(p.sparsity());

    return Function(name, ret_in, ret_out, nlpsol_in(), nlpsol_out(),
                    {{"default_in", nlpsol_default_in()}});
  }

  Function qpsol(const std::string& name, const std::string& solver,
                 const SXDict& qp, const Dict& opts) {
    return qpsol_nlp(name, solver, qp, opts);
  }

  Function qpsol(const std::string& name, const std::string& solver,
                 const MXDict& qp, const Dict& opts) {
    return qpsol_nlp(name, solver, qp, opts);
  }

  // Constructor
  Conic::Conic(const std::string& name, const std::map<std::string, Sparsity> &st)
    : FunctionInternal(name) {

    P_ = Sparsity(0, 0);
    for (auto i=st.begin(); i!=st.end(); ++i) {
      if (i->first=="a") {
        A_ = i->second;
      } else if (i->first=="h") {
        H_ = i->second;
      } else if (i->first=="q") {
        Q_ = i->second;
      } else if (i->first=="p") {
        P_ = i->second;
      } else {
        casadi_error("Unrecognized field in QP structure: " + str(i->first));
      }
    }

    // We need either A or H
    casadi_assert(!A_.is_null() || !H_.is_null(),
      "Cannot determine dimension");

    // Generate A or H
    if (A_.is_null()) {
      A_ = Sparsity(0, H_.size2());
    } else if (H_.is_null()) {
      H_ = Sparsity(A_.size2(), A_.size2());
    } else {
      // Consistency check
      casadi_assert(A_.size2()==H_.size2(),
        "Got incompatible dimensions.\n"
        "min x'Hx + G'x s.t. LBA <= Ax <= UBA :\n"
        "H: " + H_.dim() + " - A: " + A_.dim() + "\n"
        "We need: H.size2()==A.size2()");
    }

    casadi_assert(H_.is_symmetric(),
      "Got incompatible dimensions. min x'Hx + G'x\n"
      "H: " + H_.dim() +
      "We need H square & symmetric");

    nx_ = A_.size2();
    na_ = A_.size1();

    // Check psd constraints (linear part)
    if (Q_.is_null()) {
      Q_ = Sparsity(0, 0);
      np_ = 0;
    } else {
      casadi_assert(Q_.size2()==nx_,
        "Got incompatible dimensions.\n"
        "Q: " + Q_.dim() +
        "We need the product Qx to exist.");
      np_ = static_cast<casadi_int>(sqrt(static_cast<double>(Q_.size1())));
      casadi_assert(np_*np_==Q_.size1(),
        "Got incompatible dimensions.\n"
        "Q: " + Q_.dim() +
        "We need Q.size1() to have an integer square root.");

      Sparsity qsum = reshape(sum2(Q_), np_, np_);

      casadi_assert(qsum.is_symmetric(),
        "Got incompatible dimensions.");
    }

    // Check psd constraints (constant part)
    if (P_.is_null()) P_ = Sparsity(np_, np_);

    casadi_assert(P_.is_symmetric(),
      "Got incompatible dimensions.\n"
      "P: " + P_.dim() +
      "We need P square & symmetric.");

    casadi_assert(P_.size1()==np_,
      "Got incompatible dimensions.\n"
      "P: " + P_.dim() +
      "We need P " + str(np_) + "-by-" + str(np_) + ".");

  }

  Sparsity Conic::get_sparsity_in(casadi_int i) {
    switch (static_cast<ConicInput>(i)) {
    case CONIC_X0:
    case CONIC_G:
    case CONIC_LBX:
    case CONIC_UBX:
    case CONIC_LAM_X0:
      return get_sparsity_out(CONIC_X);
    case CONIC_Q:
      return Q_;
    case CONIC_P:
      return P_;
    case CONIC_LBA:
    case CONIC_UBA:
    case CONIC_LAM_A0:
      return get_sparsity_out(CONIC_LAM_A);
    case CONIC_A:
      return A_;
    case CONIC_H:
      return H_;
    case CONIC_NUM_IN: break;
    }
    return Sparsity();
  }

  Sparsity Conic::get_sparsity_out(casadi_int i) {
    switch (static_cast<ConicOutput>(i)) {
    case CONIC_COST:
      return Sparsity::scalar();
    case CONIC_X:
    case CONIC_LAM_X:
      return Sparsity::dense(nx_, 1);
    case CONIC_LAM_A:
      return Sparsity::dense(na_, 1);
    case CONIC_NUM_OUT: break;
    }
    return Sparsity();
  }

  const Options Conic::options_
  = {{&FunctionInternal::options_},
     {{"discrete",
       {OT_BOOLVECTOR,
        "Indicates which of the variables are discrete, i.e. integer-valued"}},
      {"print_problem",
       {OT_BOOL,
        "Print a numeric description of the problem"}},
      {"error_on_fail",
       {OT_BOOL,
        "When the numerical process returns unsuccessfully, raise an error (default false)."}}
     }
  };

  void Conic::init(const Dict& opts) {
    // Call the init method of the base class
    FunctionInternal::init(opts);

    print_problem_ = false;
    error_on_fail_ = true;

    // Read options
    for (auto&& op : opts) {
      if (op.first=="discrete") {
        discrete_ = op.second;
      } else if (op.first=="print_problem") {
        print_problem_ = op.second;
      } else if (op.first=="error_on_fail") {
        error_on_fail_ = op.second;
      }
    }

    // Check options
    if (!discrete_.empty()) {
      casadi_assert(discrete_.size()==nx_, "\"discrete\" option has wrong length");
      if (std::find(discrete_.begin(), discrete_.end(), true)!=discrete_.end()) {
        casadi_assert(integer_support(),
                              "Discrete variables require a solver with integer support");
      }
    }

    casadi_assert(np_==0 || psd_support(),
      "Selected solver does not support psd constraints.");

  }

  /** \brief Initalize memory block */
  int Conic::init_mem(void* mem) const {
    if (ProtoFunction::init_mem(mem)) return 1;

    auto m = static_cast<ConicMemory*>(mem);

    // Problem has not been solved at this point
    m->success = false;
    m->unified_return_status = SOLVER_RET_UNKNOWN;
    m->iter_count = -1;

    return 0;
  }

  /** \brief Set the (persistent) work vectors */
  void Conic::set_work(void* mem, const double**& arg, double**& res,
                          casadi_int*& iw, double*& w) const {

    auto m = static_cast<ConicMemory*>(mem);

    // Problem has not been solved at this point
    m->success = false;
    m->unified_return_status = SOLVER_RET_UNKNOWN;
    m->iter_count = -1;
  }

  Conic::~Conic() {
  }

  void Conic::check_inputs(const double* lbx, const double* ubx,
                          const double* lba, const double* uba) const {
    for (casadi_int i=0; i<nx_; ++i) {
      double lb = lbx ? lbx[i] : 0., ub = ubx ? ubx[i] : 0.;
      casadi_assert(lb <= ub && lb!=inf && ub!=-inf,
        "Ill-posed problem detected: "
        "LBX[" + str(i) + "] <= UBX[" + str(i) + "] was violated. "
        "Got LBX[" + str(i) + "]=" + str(lb) + " and UBX[" + str(i) + "] = " + str(ub) + ".");
    }
    for (casadi_int i=0; i<na_; ++i) {
      double lb = lba ? lba[i] : 0., ub = uba ? uba[i] : 0.;
      casadi_assert(lb <= ub && lb!=inf && ub!=-inf,
        "Ill-posed problem detected: "
        "LBA[" + str(i) + "] <= UBA[" + str(i) + "] was violated. "
        "Got LBA[" + str(i) + "] = " + str(lb) + " and UBA[" + str(i) + "] = " + str(ub) + ".");
    }
  }

  void Conic::generateNativeCode(std::ostream& file) const {
    casadi_error("generateNativeCode not defined for class " + class_name());
  }

  std::map<std::string, Conic::Plugin> Conic::solvers_;

  const std::string Conic::infix_ = "conic";

  double Conic::get_default_in(casadi_int ind) const {
    switch (ind) {
    case CONIC_LBX:
    case CONIC_LBA:
      return -std::numeric_limits<double>::infinity();
    case CONIC_UBX:
    case CONIC_UBA:
      return std::numeric_limits<double>::infinity();
    default:
      return 0;
    }
  }

  int Conic::
  eval(const double** arg, double** res, casadi_int* iw, double* w, void* mem) const {
    if (print_problem_) {
      uout() << "H:";
      DM::print_dense(uout(), H_, arg[CONIC_H], false);
      uout() << std::endl;
      uout() << "G:" << std::vector<double>(arg[CONIC_G], arg[CONIC_G]+nx_) << std::endl;
      uout() << "A:";
      DM::print_dense(uout(), A_, arg[CONIC_A], false);
      uout() << std::endl;
      uout() << "lba:" << std::vector<double>(arg[CONIC_LBA], arg[CONIC_LBA]+na_) << std::endl;
      uout() << "uba:" << std::vector<double>(arg[CONIC_UBA], arg[CONIC_UBA]+na_) << std::endl;
      uout() << "lbx:" << std::vector<double>(arg[CONIC_LBX], arg[CONIC_LBX]+nx_) << std::endl;
      uout() << "ubx:" << std::vector<double>(arg[CONIC_UBX], arg[CONIC_UBX]+nx_) << std::endl;
    }
    auto m = static_cast<ConicMemory*>(mem);

    if (inputs_check_) {
      check_inputs(arg[CONIC_LBX], arg[CONIC_UBX], arg[CONIC_LBA], arg[CONIC_UBA]);
    }
    int ret = solve(arg, res, iw, w, mem);

    if (error_on_fail_ && !m->success)
      casadi_error("conic process failed. "
                   "Set 'error_on_fail' option to false to ignore this error.");
    return ret;
  }

  std::vector<std::string> conic_options(const std::string& name) {
    return Conic::plugin_options(name).all();
  }

  std::string conic_option_type(const std::string& name, const std::string& op) {
    return Conic::plugin_options(name).type(op);
  }

  std::string conic_option_info(const std::string& name, const std::string& op) {
    return Conic::plugin_options(name).info(op);
  }

  bool Conic::is_a(const std::string& type, bool recursive) const {
    return type=="Conic" || (recursive && FunctionInternal::is_a(type, recursive));
  }

  void Conic::sdp_to_socp_init(SDPToSOCPMem& mem) const {

    Sparsity qsum = reshape(sum2(Q_), np_, np_);

    // Block detection
    Sparsity aggregate = qsum+P_;

    std::vector<casadi_int> p;
    casadi_int nb = aggregate.scc(p, mem.r);

    std::string pattern_message = "Pattern not recognised";

    casadi_assert(p==range(p.size()), pattern_message);

    const casadi_int* row = aggregate.row();
    const casadi_int* colind = aggregate.colind();

    // Check fishbone-structure
    for (casadi_int i=0;i<nb;++i) {
      casadi_int block_size = mem.r[i+1]-mem.r[i];
      // number of nonzeros in column mem.r[i+1]-1
      casadi_int nz = colind[mem.r[i+1]]-colind[mem.r[i+1]-1];
      // Last column of block should be dense
      casadi_assert(nz==block_size, pattern_message);
      for (casadi_int k=0;k<block_size-1;++k) {
        casadi_assert(colind[mem.r[i]+k+1]-colind[mem.r[i]+k], pattern_message);
        casadi_assert(*(row++)==k+mem.r[i], pattern_message);
        casadi_assert(*(row++)==mem.r[i]+block_size-1, pattern_message);
      }

      for (casadi_int k=0;k<block_size;++k)
        casadi_assert(*(row++)==k+mem.r[i], pattern_message);
    }

    /**
      general soc constraints:
      ||Ax + b ||_2 <= c'x + d

      Need to represent as
      X'X <= Z'Z

      with X and Z helper variables and
      Ax  + b = X
      c'x + d = Z

      [A;c'] x - [X;Z] = [b;d]

      we look for the vertical concatenation of these constraints for all blocks:

      Q(map_Q) [x;X1;Z1;X2;Z2;...] = P.nz(map_P)

    */

    /*

    Aggregate pattern:

    x (x)
     x(x)
    xx(x)
       x (x)
        x(x)
       xx(x)

    We are interested in the parts in parenthesis (target).

    Find out which rows in Q correspond to those targets
    */

    // Lookup vector for target start and end
    std::vector<casadi_int> target_start(nb), target_stop(nb);
    for (casadi_int i=0;i<nb;++i) {
      target_start[i] = (mem.r[i+1]-1)*aggregate.size1()+mem.r[i];
      target_stop[i] = target_start[i]+mem.r[i+1]-mem.r[i];
    }

    // Collect the nonzero indices in Q that that correspond to the target area
    std::vector<casadi_int> q_nz;
    // Triplet form for map_Q sparsity
    std::vector<casadi_int> q_row, q_col;

    // Loop over Q's columns (decision variables)
    for (casadi_int j=0; j<Q_.size2(); ++j) {
      casadi_int block_index = 0;
      // Loop over Q's rows
      for (casadi_int k=Q_.colind(j); k<Q_.colind(j+1); ++k) {
        casadi_int i = Q_.row(k);

        // Increment block_index if i runs ahead
        while (i>target_stop[block_index] && block_index<nb-1) block_index++;

        if (i>=target_start[block_index] && i<target_stop[block_index]) {
          // Got a nonzero in the target region
          q_nz.push_back(k);
          q_row.push_back(mem.r[block_index]+i-target_start[block_index]);
          q_col.push_back(j);
        }
      }
    }

    mem.map_Q = IM::triplet(q_row, q_col, q_nz, mem.r[nb], nx_);

    // Add the [X1;Z1;X2;Z2;...] part
    mem.map_Q = horzcat(mem.map_Q, -IM::eye(mem.r[nb])).T();

    // Get maximum nonzero count of any column
    casadi_int max_nnz = 0;
    for (casadi_int i=0;i<mem.map_Q.size2();++i) {
      max_nnz = std::max(max_nnz, mem.map_Q.colind(i+1)-mem.map_Q.colind(i));
    }

    // ind/val size needs to cover max nonzero count
    mem.indval_size = std::max(nx_, max_nnz);

    // Collect the indices for the P target area
    mem.map_P.resize(mem.r[nb], -1);
    for (casadi_int i=0;i<nb;++i) {
      for (casadi_int k=P_.colind(mem.r[i+1]-1); k<P_.colind(mem.r[i+1]); ++k) {
        casadi_int r = P_.row(k);
        mem.map_P[r] = k;
      }
    }

    // ind/val size needs to cover blocksize
    for (casadi_int i=0;i<nb;++i)
      mem.indval_size = std::max(mem.indval_size, mem.r[i+1]-mem.r[i]);

    // Get the transpose and mapping
    mem.AT = A_.transpose(mem.A_mapping);
  }

  Dict Conic::get_stats(void* mem) const {
    Dict stats = FunctionInternal::get_stats(mem);
    auto m = static_cast<ConicMemory*>(mem);

    stats["success"] = m->success;
    stats["unified_return_status"] = string_from_UnifiedReturnStatus(m->unified_return_status);
    stats["iter_count"] = m->iter_count;
    return stats;
  }

  void Conic::serialize(SerializingStream &s, const SDPToSOCPMem& m) const {
    s.pack("Conic::SDPToSOCPMem::r", m.r);
    s.pack("Conic::SDPToSOCPMem::AT", m.AT);
    s.pack("Conic::SDPToSOCPMem::A_mapping", m.A_mapping);
    s.pack("Conic::SDPToSOCPMem::map_Q", m.map_Q);
    s.pack("Conic::SDPToSOCPMem::map_P", m.map_P);
    s.pack("Conic::SDPToSOCPMem::indval_size", m.indval_size);
  }
  void Conic::deserialize(DeserializingStream &s, SDPToSOCPMem& m) {
    s.unpack("Conic::SDPToSOCPMem::r", m.r);
    s.unpack("Conic::SDPToSOCPMem::AT", m.AT);
    s.unpack("Conic::SDPToSOCPMem::A_mapping", m.A_mapping);
    s.unpack("Conic::SDPToSOCPMem::map_Q", m.map_Q);
    s.unpack("Conic::SDPToSOCPMem::map_P", m.map_P);
    s.unpack("Conic::SDPToSOCPMem::indval_size", m.indval_size);
  }

  void Conic::serialize_body(SerializingStream &s) const {
    FunctionInternal::serialize_body(s);

    s.version("Conic", 1);
    s.pack("Conic::discrete", discrete_);
    s.pack("Conic::print_problem", print_problem_);
    s.pack("Conic::error_on_fail", error_on_fail_);
    s.pack("Conic::H", H_);
    s.pack("Conic::A", A_);
    s.pack("Conic::Q", Q_);
    s.pack("Conic::P", P_);
    s.pack("Conic::nx", nx_);
    s.pack("Conic::na", na_);
    s.pack("Conic::np", np_);
  }

  void Conic::serialize_type(SerializingStream &s) const {
    FunctionInternal::serialize_type(s);
    PluginInterface<Conic>::serialize_type(s);
  }

  ProtoFunction* Conic::deserialize(DeserializingStream& s) {
    return PluginInterface<Conic>::deserialize(s);
  }

  Conic::Conic(DeserializingStream & s) : FunctionInternal(s) {
    s.version("Conic", 1);
    s.unpack("Conic::discrete", discrete_);
    s.unpack("Conic::print_problem", print_problem_);
    s.unpack("Conic::error_on_fail", error_on_fail_);
    s.unpack("Conic::H", H_);
    s.unpack("Conic::A", A_);
    s.unpack("Conic::Q", Q_);
    s.unpack("Conic::P", P_);
    s.unpack("Conic::nx", nx_);
    s.unpack("Conic::na", na_);
    s.unpack("Conic::np", np_);
  }

} // namespace casadi
