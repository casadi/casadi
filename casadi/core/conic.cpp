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


#include "conic_impl.hpp"
#include "nlpsol_impl.hpp"
#include "filesystem_impl.hpp"

namespace casadi {

  bool has_conic(const std::string& name) {
    return Conic::has_plugin(name);
  }

  void load_conic(const std::string& name) {
    Conic::load_plugin(name);
  }

  std::string doc_conic(const std::string& name) {
    return Conic::getPlugin(name).doc;
  }

  Function conic(const std::string& name, const std::string& solver,
                const SpDict& qp, const Dict& opts) {
    return Function::create(Conic::instantiate(name, solver, qp), opts);
  }

  void conic_debug(const Function& f, const std::string &filename) {
    auto file_ptr = Filesystem::ofstream_ptr(filename);
    conic_debug(f, *file_ptr);
  }

  void conic_debug(const Function& f, std::ostream &file) {
    casadi_assert_dev(!f.is_null());
    const Conic* n = f.get<Conic>();
    n->generateNativeCode(file);
  }

  std::vector<std::string> conic_in() {
    std::vector<std::string> ret(conic_n_in());
    for (size_t i=0; i<ret.size(); ++i) ret[i]=conic_in(i);
    return ret;
  }

  std::vector<std::string> conic_out() {
    std::vector<std::string> ret(conic_n_out());
    for (size_t i=0; i<ret.size(); ++i) ret[i]=conic_out(i);
    return ret;
  }

  std::string conic_in(casadi_int ind) {
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
    return std::string();
  }

  std::string conic_out(casadi_int ind) {
    switch (static_cast<ConicOutput>(ind)) {
    case CONIC_X:     return "x";
    case CONIC_COST:  return "cost";
    case CONIC_LAM_A: return "lam_a";
    case CONIC_LAM_X: return "lam_x";
    case CONIC_NUM_OUT: break;
    }
    return std::string();
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
    bool expand = false;
    Dict opt = opts;
    extract_from_dict_inplace(opt, "expand", expand);
    bool postpone_expand = false;
    extract_from_dict_inplace(opt, "postpone_expand", postpone_expand);
    bool error_on_fail = get_from_dict(opts, "error_on_fail", true);

    if (expand && !postpone_expand && M::type_name()=="MX") {
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
    Function prob(name + "_qp", {x, p}, {H, c, A, b, Q, P},
      {"x", "p"}, {"H", "c", "A", "b", "Q", "P"});
    if (expand && postpone_expand) prob = prob.expand();

    // Make sure that the problem is sound
    casadi_assert(!prob.has_free(), "Cannot create '" + prob.name() + "' "
                          "since " + str(prob.get_free()) + " are free.");

    // Create the QP solver
    Function conic_f = conic(name + "_qpsol", solver,
                             {{"h", H.sparsity()}, {"a", A.sparsity()},
                              {"p", P.sparsity()}, {"q", Q.sparsity()}}, opt);

    // Create an MXFunction with the right signature
    std::vector<MX> ret_in(NLPSOL_NUM_IN);
    ret_in[NLPSOL_X0] = MX::sym("x0", x.sparsity());
    ret_in[NLPSOL_P] = MX::sym("p", p.sparsity());
    ret_in[NLPSOL_LBX] = MX::sym("lbx", x.sparsity());
    ret_in[NLPSOL_UBX] = MX::sym("ubx", x.sparsity());
    ret_in[NLPSOL_LBG] = MX::sym("lbg", g.sparsity());
    ret_in[NLPSOL_UBG] = MX::sym("ubg", g.sparsity());
    ret_in[NLPSOL_LAM_X0] = MX::sym("lam_x0", x.sparsity());
    ret_in[NLPSOL_LAM_G0] = MX::sym("lam_g0", g.sparsity());
    std::vector<MX> ret_out(NLPSOL_NUM_OUT);


    std::vector<MX> v(NL_NUM_IN);
    v[NL_X] = ret_in[NLPSOL_X0];
    v[NL_P] = ret_in[NLPSOL_P];
    // Evaluate constant part of objective
    MX rv = r(v)[0];
    // Get expressions for the QP matrices and vectors
    v = prob(v);

    // Call the QP solver
    std::vector<MX> w(CONIC_NUM_IN);
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

    Dict fun_opts;
    fun_opts["default_in"] = nlpsol_default_in();
    fun_opts["error_on_fail"] = error_on_fail;

    return Function(name, ret_in, ret_out, nlpsol_in(), nlpsol_out(), fun_opts);
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

    // Set default options
    error_on_fail_ = true;

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
      {"equality",
       {OT_BOOLVECTOR,
        "Indicate an upfront hint which of the constraints are equalities. "
        "Some solvers may be able to exploit this knowledge. "
        "When true, the corresponding lower and upper bounds are assumed equal. "
        "When false, the corresponding bounds may be equal or different."}},
      {"print_problem",
       {OT_BOOL,
        "Print a numeric description of the problem"}},
      {"solver_version_check",
       {OT_BOOL,
        "When the plugin loads an externally supplied solver, "
         "check that its version is compatible with the plugin [Default: true]"}},
      {"condense",
       {OT_BOOL,
        "Apply OCP partial-condensing before passing to the solver. "
        "Requires the QP to have OCP block structure (use 'structure_detection')."}},
      {"structure_detection",
       {OT_STRING,
        "OCP structure detection: 'none' (default), 'auto', or 'manual'."}},
      {"N",
       {OT_INT,
        "OCP horizon (manual structure_detection)"}},
      {"nx",
       {OT_INTVECTOR,
        "Number of states per stage, length N+1 (manual structure_detection)"}},
      {"nu",
       {OT_INTVECTOR,
        "Number of controls per stage, length N (manual structure_detection)"}},
      {"ng",
       {OT_INTVECTOR,
        "Number of path inequalities per stage, length N+1 (manual structure_detection)"}},
      {"condense_partition",
       {OT_INTVECTOR,
        "Partition M = [0, M1, ..., N] for partial condensing.  "
        "Expert override -- if set, the DP that picks M from "
        "condensed_block_count is skipped.  Mutually exclusive "
        "with condensed_block_count."}},
      {"condensed_block_count",
       {OT_INT,
        "Number of condensed blocks (n_hat).  Default 0 means the "
        "selected partition strategy picks both n_hat and M to "
        "minimize Riccati solve cost.  k > 0 fixes n_hat=k.  "
        "Mutually exclusive with condense_partition."}},
      {"condense_partition_strategy",
       {OT_STRING,
        "How to derive the partition M from condensed_block_count: "
        "'auto' (default) picks 'uniform' if nx is uniform across "
        "stages, else 'optimal' (with a warning if N>500); "
        "'uniform' assumes uniform nx and uses Frison's closed-form / "
        "equal-split partition; 'optimal' runs an O(K_max*N^2) DP."}},
      {"debug",
       {OT_BOOL,
        "Produce debug information for structure detection / condensing "
        "(default: false).  When auto-detect runs, dumps the A sparsity "
        "to debug_conic_actual.mtx and prints inferred N, nx, nu, ng."}}
     }
  };

  void Conic::init(const Dict& opts) {
    // Call the init method of the base class
    FunctionInternal::init(opts);

    print_problem_ = false;
    solver_version_check_ = true;
    debug_ = false;
    condensed_block_count_ = 0;
    condense_partition_strategy_ = "auto";
    condense_ = false;
    condense_structure_detection_ = COND_STRUCT_NONE;
    N_ = 0;
    N_hat_ = 0;
    nx_total_hat_ = 0;
    na_total_hat_ = 0;
    nx_max_ = 0;
    nu_max_block_ = 0;
    nxu_max_block_ = 0;

    casadi_int struct_cnt = 0;

    // Read options
    for (auto&& op : opts) {
      if (op.first=="discrete") {
        discrete_ = op.second;
      } else if (op.first=="equality") {
        equality_ = op.second;
      } else if (op.first=="print_problem") {
        print_problem_ = op.second;
      } else if (op.first=="solver_version_check") {
        solver_version_check_ = op.second;
      } else if (op.first=="condense") {
        condense_ = op.second;
      } else if (op.first=="structure_detection") {
        std::string v = op.second;
        if (v=="none") {
          condense_structure_detection_ = COND_STRUCT_NONE;
        } else if (v=="auto") {
          condense_structure_detection_ = COND_STRUCT_AUTO;
        } else if (v=="manual") {
          condense_structure_detection_ = COND_STRUCT_MANUAL;
        } else {
          casadi_error("Unknown 'structure_detection': '" + v + "'.");
        }
      } else if (op.first=="N") {
        N_ = op.second;
        struct_cnt++;
      } else if (op.first=="nx") {
        nxs_ = op.second.to_int_vector();
        struct_cnt++;
      } else if (op.first=="nu") {
        nus_ = op.second.to_int_vector();
        struct_cnt++;
      } else if (op.first=="ng") {
        ngs_ = op.second.to_int_vector();
        struct_cnt++;
      } else if (op.first=="condense_partition") {
        M_user_ = op.second.to_int_vector();
        M_ = M_user_;
      } else if (op.first=="condensed_block_count") {
        condensed_block_count_ = op.second;
      } else if (op.first=="condense_partition_strategy") {
        condense_partition_strategy_ = std::string(op.second);
        casadi_assert(
          condense_partition_strategy_ == "auto" ||
          condense_partition_strategy_ == "uniform" ||
          condense_partition_strategy_ == "optimal",
          "condense_partition_strategy must be 'auto', 'uniform', or "
          "'optimal'; got '" + condense_partition_strategy_ + "'.");
      } else if (op.first=="debug") {
        debug_ = op.second;
      }
    }
    casadi_assert(
      M_user_.empty() || condensed_block_count_ == 0,
      "condense_partition and condensed_block_count are mutually "
      "exclusive; set at most one.");

    if (solver_version_check_) deps_version_check("init");

    // Check options
    if (!discrete_.empty()) {
      casadi_assert(discrete_.size()==nx_, "\"discrete\" option has wrong length");
      if (std::find(discrete_.begin(), discrete_.end(), true)!=discrete_.end()) {
        casadi_assert(integer_support(),
                              "Discrete variables require a solver with integer support");
      }
    }

    if (!equality_.empty()) {
      casadi_assert(equality_.size()==na_, "\"equality\" option has wrong length. "
                                           "Expected " + str(na_) + " elements, but got " +
                                            str(equality_.size()) + " instead.");
    }

    casadi_assert(np_==0 || psd_support(),
      "Selected solver does not support psd constraints.");

    if (condense_) {
      casadi_assert(condense_structure_detection_ != COND_STRUCT_NONE,
        "condense=True requires structure_detection='auto' or 'manual'.");
      if (condense_structure_detection_ == COND_STRUCT_MANUAL) {
        casadi_assert(struct_cnt == 4,
          "structure_detection='manual' requires N, nx, nu, ng to be set.");
      } else if (condense_structure_detection_ == COND_STRUCT_AUTO) {
        casadi_assert(struct_cnt == 0,
          "structure_detection='auto' must not be combined with manual N/nx/nu/ng.");
      }
      detect_condense_structure();
      build_condense_blocks();

      // Reserve workspace -- single source of truth via casadi_condensing_work
      // (filled fields in p_cond_ tell it the totals).
      casadi_int sz_iw = 0, sz_w = 0;
      casadi_condensing_work(&p_cond_, &sz_iw, &sz_w);
      // Slack for projections (rows of target sparsities).
      sz_w += std::max({RSQsp_.size1(), ABsp_.size1(), CDsp_.size1()});
      alloc_w(sz_w, true);
      alloc_iw(sz_iw, true);
    }

    set_qp_prob();
  }

  void Conic::finalize() {
    if (solver_version_check_) deps_version_check("finalize");

    // Recursive call
    FunctionInternal::finalize();
  }

  /** \brief Initalize memory block */
  int Conic::init_mem(void* mem) const {
    if (ProtoFunction::init_mem(mem)) return 1;
    if (condense_) {
      auto *m = static_cast<ConicMemory*>(mem);
      m->add_stat("block_transform");  // user H/A -> per-stage flat layout
      m->add_stat("condensing");       // condense math (eval)
      m->add_stat("lifting");          // lift + copy-out
    }
    return 0;
  }

  /** \brief Set the (persistent) work vectors */
  void Conic::set_work(void* mem, const double**& arg, double**& res,
                          casadi_int*& iw, double*& w) const {

    auto *m = static_cast<ConicMemory*>(mem);

    casadi_qp_data<double>& d_qp = m->d_qp;

    if (condense_) {
      // ---- C++ runtime condensing wrap ----
      // casadi_condensing_set_work claims ALL workspace (per-stage flat in/
      // out, internal scratch, condensed CSC, condensed bounds, primal/
      // dual scratch).  Then we project user inputs into the per-stage
      // flat layout (block_transform), run the condense math (condensing),
      // and override d_qp to point at the condensed problem.
      casadi_condensing_data<double>& d_cond = m->d_cond;
      d_cond.prob = &p_cond_;
      casadi_condensing_set_work(&d_cond, &arg, &res, &iw, &w);

      // (1) Block transform: turn user H, A and bounds/snapshots into
      // the per-stage flat input layout that condensing_eval consumes.
      m->fstats.at("block_transform").tic();
      casadi_project(arg[CONIC_H], H_, d_cond.RSQ_val, RSQsp_, w);  /* arg[0] H   */
      d_cond.g_orig = arg[CONIC_G];                                 /* arg[1] g   */
      casadi_project(arg[CONIC_A], A_, d_cond.AB_val, ABsp_, w);    /* arg[2] A   */
      if (CDsp_.nnz() > 0) {
        casadi_project(arg[CONIC_A], A_, d_cond.CD_val, CDsp_, w);  /* arg[2] A,CD*/
      }
      d_cond.lba_orig = arg[CONIC_LBA];                             /* arg[3] lba */
      d_cond.uba_orig = arg[CONIC_UBA];                             /* arg[4] uba */
      d_cond.lbx_orig = arg[CONIC_LBX];                             /* arg[5] lbx */
      d_cond.ubx_orig = arg[CONIC_UBX];                             /* arg[6] ubx */
      m->fstats.at("block_transform").toc();

      // (2) Condense math.
      m->fstats.at("condensing").tic();
      casadi_condensing_eval(&d_cond);
      m->fstats.at("condensing").toc();

      // Override d_qp to point at the condensed problem
      d_qp.h     = d_cond.h_hat_csc;
      d_qp.g     = d_cond.qr_hat_val;
      d_qp.a     = d_cond.a_hat_csc;
      d_qp.lbx   = d_cond.lbx;
      d_qp.ubx   = d_cond.ubx;
      d_qp.lba   = d_cond.lba;
      d_qp.uba   = d_cond.uba;
      d_qp.x0    = nullptr;
      d_qp.lam_x0 = nullptr;
      d_qp.lam_a0 = nullptr;
      d_qp.x     = d_cond.x;
      d_qp.lam_x = d_cond.lam_x;
      d_qp.lam_a = d_cond.lam_a;
      d_qp.f     = res[CONIC_COST];
    } else {
      d_qp.h = arg[CONIC_H];
      d_qp.g = arg[CONIC_G];
      d_qp.a = arg[CONIC_A];
      d_qp.lbx = arg[CONIC_LBX];
      d_qp.ubx = arg[CONIC_UBX];
      d_qp.uba = arg[CONIC_UBA];
      d_qp.lba = arg[CONIC_LBA];
      d_qp.x0 = arg[CONIC_X0];
      d_qp.lam_x0 = arg[CONIC_LAM_X0];
      d_qp.lam_a0 = arg[CONIC_LAM_A0];
      d_qp.x = res[CONIC_X];
      d_qp.lam_x = res[CONIC_LAM_X];
      d_qp.lam_a = res[CONIC_LAM_A];
      d_qp.f = res[CONIC_COST];
    }

    // Problem has not been solved at this point
    d_qp.success = false;
    d_qp.unified_return_status = SOLVER_RET_UNKNOWN;
    d_qp.iter_count = -1;
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

#ifdef CASADI_WITH_THREADSAFE_SYMBOLICS
  std::mutex Conic::mutex_solvers_;
#endif // CASADI_WITH_THREADSAFE_SYMBOLICS

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
    auto *m = static_cast<ConicMemory*>(mem);

    if (inputs_check_) {
      check_inputs(arg[CONIC_LBX], arg[CONIC_UBX], arg[CONIC_LBA], arg[CONIC_UBA]);
    }

    setup(mem, arg, res, iw, w);

    int ret = solve(arg, res, iw, w, mem);

    if (condense_) {
      // Lift condensed primal/dual; copy into user-facing res[CONIC_*].
      m->fstats.at("lifting").tic();
      casadi_condensing_lift(&m->d_cond);
      casadi_copy(m->d_cond.x_lifted,     p_cond_.total_qr,
                  res[CONIC_X]);
      casadi_copy(m->d_cond.lam_x_lifted, p_cond_.total_qr,
                  res[CONIC_LAM_X]);
      casadi_copy(m->d_cond.lam_a_lifted, p_cond_.total_b + p_cond_.total_g,
                  res[CONIC_LAM_A]);
      m->fstats.at("lifting").toc();
    }

    if (m->d_qp.success) m->d_qp.unified_return_status = SOLVER_RET_SUCCESS;

    if (error_on_fail_ && !m->d_qp.success)
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
    auto *m = static_cast<ConicMemory*>(mem);

    stats["success"] = m->d_qp.success;
    stats["unified_return_status"] = string_from_UnifiedReturnStatus(m->d_qp.unified_return_status);
    stats["iter_count"] = m->d_qp.iter_count;
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

    s.version("Conic", 5);
    s.pack("Conic::discrete", discrete_);
    s.pack("Conic::equality", equality_);
    s.pack("Conic::print_problem", print_problem_);
    s.pack("Conic::solver_version_check", solver_version_check_);

    s.pack("Conic::H", H_);
    s.pack("Conic::A", A_);
    s.pack("Conic::Q", Q_);
    s.pack("Conic::P", P_);
    s.pack("Conic::nx", nx_);
    s.pack("Conic::na", na_);
    s.pack("Conic::np", np_);

    // Condense feature: pack the user-facing inputs.  Derived fields
    // (M_, *_blocks_packed_, sparsities, etc.) are recomputed in the
    // deserializing constructor by re-running the same init logic.
    s.pack("Conic::condense", condense_);
    s.pack("Conic::condense_structure_detection",
           static_cast<casadi_int>(condense_structure_detection_));
    s.pack("Conic::N", N_);
    s.pack("Conic::nxs", nxs_);
    s.pack("Conic::nus", nus_);
    s.pack("Conic::ngs", ngs_);
    s.pack("Conic::M_user", M_user_);
    s.pack("Conic::condensed_block_count", condensed_block_count_);
    s.pack("Conic::condense_partition_strategy", condense_partition_strategy_);
    s.pack("Conic::debug", debug_);
  }

  void Conic::serialize_type(SerializingStream &s) const {
    FunctionInternal::serialize_type(s);
    PluginInterface<Conic>::serialize_type(s);
  }

  ProtoFunction* Conic::deserialize(DeserializingStream& s) {
    return PluginInterface<Conic>::deserialize(s);
  }

  Conic::Conic(DeserializingStream & s) : FunctionInternal(s) {
    int version = s.version("Conic", 1, 5);
    s.unpack("Conic::discrete", discrete_);
    if (version>=3) {
      s.unpack("Conic::equality", equality_);
    }
    s.unpack("Conic::print_problem", print_problem_);
    if (version>=4) {
      s.unpack("Conic::solver_version_check", solver_version_check_);
    } else {
      solver_version_check_ = true;
    }
    if (version==1) {
      s.unpack("Conic::error_on_fail", error_on_fail_);
    }

    s.unpack("Conic::H", H_);
    s.unpack("Conic::A", A_);
    s.unpack("Conic::Q", Q_);
    s.unpack("Conic::P", P_);
    s.unpack("Conic::nx", nx_);
    s.unpack("Conic::na", na_);
    s.unpack("Conic::np", np_);

    // Condense feature: recovered for version >= 5; older streams default
    // to non-condense (the only mode that existed before).
    if (version >= 5) {
      s.unpack("Conic::condense", condense_);
      casadi_int sd_int;
      s.unpack("Conic::condense_structure_detection", sd_int);
      condense_structure_detection_ =
          static_cast<CondenseStructureDetection>(sd_int);
      s.unpack("Conic::N", N_);
      s.unpack("Conic::nxs", nxs_);
      s.unpack("Conic::nus", nus_);
      s.unpack("Conic::ngs", ngs_);
      s.unpack("Conic::M_user", M_user_);
      M_ = M_user_;
      s.unpack("Conic::condensed_block_count", condensed_block_count_);
      s.unpack("Conic::condense_partition_strategy",
               condense_partition_strategy_);
      s.unpack("Conic::debug", debug_);
    } else {
      condense_ = false;
      condense_structure_detection_ = COND_STRUCT_NONE;
      N_ = 0; N_hat_ = 0;
      condensed_block_count_ = 0;
      condense_partition_strategy_ = "auto";
      debug_ = false;
    }

    // Re-derive condense data from the unpacked state.  We skip
    // detect_condense_structure() because the auto-detect already
    // ran at original init time and its results are captured in
    // nxs_/nus_/ngs_/N_; nus_ is already padded with the trailing 0.
    // build_condense_blocks() handles both the user-supplied-M case
    // and the derive-via-strategy case via its M_.empty() check.
    if (condense_) {
      build_condense_blocks();
    }

    set_qp_prob();
  }

  void Conic::set_qp_prob() {
    if (condense_) {
      // Plugin sees the CONDENSED problem; set p_qp_ accordingly so its
      // own work() / init() / densify use condensed nx/na.
      p_qp_.sp_a = A_hat_sp_;
      p_qp_.sp_h = H_hat_sp_;
    } else {
      p_qp_.sp_a = A_;
      p_qp_.sp_h = H_;
    }
    casadi_qp_setup(&p_qp_);
  }

  // Build a sparsity from a list of (offset_r, offset_c, rows, cols) blocks,
  // optionally as identity blocks (eye=true).
  static Sparsity conic_blocksparsity(casadi_int rows, casadi_int cols,
      const std::vector<casadi_int>& packed, bool eye) {
    DM r(rows, cols);
    casadi_int n = packed.size() / 4;
    for (casadi_int i = 0; i < n; ++i) {
      casadi_int orow = packed[4*i + 0];
      casadi_int ocol = packed[4*i + 1];
      casadi_int rr = packed[4*i + 2];
      casadi_int cc = packed[4*i + 3];
      if (eye) {
        r(range(orow, orow+rr), range(ocol, ocol+cc)) = DM::eye(rr);
        casadi_assert_dev(rr == cc);
      } else {
        r(range(orow, orow+rr), range(ocol, ocol+cc)) = DM::zeros(rr, cc);
      }
    }
    return r.sparsity();
  }

  // Auto-detect OCP structure from A_ (port of fatrop_conic logic).
  void Conic::detect_condense_structure() {
    if (condense_structure_detection_ == COND_STRUCT_AUTO) {
      if (debug_) {
        A_.to_file("debug_conic_actual.mtx");
        uout() << "Conic auto-detect: A=" << A_.size1() << "x" << A_.size2()
               << " nnz=" << A_.nnz()
               << " (dumped to debug_conic_actual.mtx)" << std::endl;
      }
      // The skyline walker below indexes into A_skyline[0]/nus_[0] after
      // the loop -- precondition is at least one row in A_ (one dynamics
      // or path-constraint row).  An empty A makes those reads OOB.
      casadi_assert(na_ > 0,
        "structure_detection='auto' requires A to have at least one row "
        "(dynamics or path constraints).  Got A: " + str(A_.size1()) + "x" +
        str(A_.size2()) + ".  For QPs without OCP structure, set "
        "structure_detection='manual' and provide N/nx/nu/ng explicitly, "
        "or set condense=false.");
      Sparsity AT = A_.T();
      std::vector<casadi_int> A_skyline, A_skyline2, A_bottomline;
      for (casadi_int i = 0; i < AT.size2(); ++i) {
        casadi_int pivot = AT.colind()[i+1];
        A_bottomline.push_back(AT.row()[AT.colind()[i]]);
        if (pivot > AT.colind()[i]) {
          A_skyline.push_back(AT.row()[pivot-1]);
          if (pivot > AT.colind()[i] + 1) {
            A_skyline2.push_back(AT.row()[pivot-2]);
          } else {
            A_skyline2.push_back(-1);
          }
        } else {
          A_skyline.push_back(-1);
          A_skyline2.push_back(-1);
        }
      }
      nus_.clear();
      nxs_.clear();
      ngs_.clear();
      casadi_int pivot = 0;
      casadi_int start_pivot = pivot;
      casadi_int cg = 0;
      for (casadi_int i = 0; i < na_; ++i) {
        bool commit = false;
        if (A_skyline[i] > pivot + 1) {
          nus_.push_back(A_skyline[i] - pivot - 1);
          commit = true;
        } else if (A_skyline[i] == pivot + 1) {
          if (A_skyline2[i] < start_pivot) {
            pivot++;
          } else {
            nus_.push_back(0);
            commit = true;
          }
        } else {
          cg++;
        }
        if (commit) {
          nxs_.push_back(pivot - start_pivot + 1);
          ngs_.push_back(cg); cg = 0;
          start_pivot = A_skyline[i];
          pivot = A_skyline[i];
        }
      }
      nxs_.push_back(pivot - start_pivot + 1);
      // Correction for k==0
      nxs_[0] = A_skyline[0];
      nus_[0] = 0;
      ngs_.erase(ngs_.begin());
      casadi_int cN = 0;
      for (casadi_int i = na_ - 1; i >= 0; --i) {
        if (A_bottomline[i] < start_pivot) break;
        cN++;
      }
      ngs_.push_back(cg - cN);
      ngs_.push_back(cN);
      N_ = nus_.size();
      nus_.push_back(0);
      if (N_ > 1) {
        if (nus_[0] == 0 && nxs_[1] + nus_[1] == nxs_[0]) {
          nxs_[0] = nxs_[1];
          nus_[0] = nus_[1];
        }
      }
      if (debug_) {
        uout() << "Conic auto-detect: N=" << N_
               << " nx=" << nxs_ << " nu=" << nus_ << " ng=" << ngs_
               << std::endl;
      }
    } else {  // MANUAL
      casadi_assert((casadi_int)nxs_.size() == N_ + 1,
        "nx must have length N+1 = " + str(N_+1));
      casadi_assert((casadi_int)nus_.size() == N_,
        "nu must have length N = " + str(N_));
      casadi_assert((casadi_int)ngs_.size() == N_ + 1,
        "ng must have length N+1 = " + str(N_+1));
      nus_.push_back(0);  // pad nu[N]=0 for unified indexing
    }

  }

  void Conic::build_condense_blocks() {
    // Derive M_ if not user-supplied via condense_partition; validate;
    // set N_hat_.  This block was previously in detect_condense_structure
    // but lives here so the deserialize path (which skips detect) gets it.
    if (M_.empty()) {
      derive_condense_partition();
    }
    casadi_assert(M_.front() == 0, "Partition M must start at 0.");
    casadi_assert(M_.back() == N_, "Partition M must end at N.");
    for (size_t i = 1; i < M_.size(); ++i) {
      casadi_assert(M_[i] > M_[i-1], "Partition M must be strictly increasing.");
    }
    N_hat_ = M_.size() - 1;
    if (debug_) {
      uout() << "Conic condense partition: strategy='"
             << condense_partition_strategy_ << "', n_hat=" << N_hat_
             << ", M=" << M_ << std::endl;
    }
    const std::vector<casadi_int>& nx = nxs_;
    const std::vector<casadi_int>& ng = ngs_;
    const std::vector<casadi_int>& nu = nus_;

    // Build user-side block descriptors
    AB_blocks_packed_.clear();
    CD_blocks_packed_.clear();
    RSQ_blocks_packed_.clear();
    AB_offsets_.clear();
    CD_offsets_.clear();
    RSQ_offsets_.clear();

    casadi_int offset_r = 0, offset_c = 0;
    AB_offsets_.push_back(0);
    CD_offsets_.push_back(0);
    casadi_int off_AB = 0, off_CD = 0;
    for (casadi_int k = 0; k < N_; ++k) {
      AB_blocks_packed_.push_back(offset_r);
      AB_blocks_packed_.push_back(offset_c);
      AB_blocks_packed_.push_back(nx[k+1]);
      AB_blocks_packed_.push_back(nx[k] + nu[k]);
      off_AB += nx[k+1] * (nx[k] + nu[k]);
      AB_offsets_.push_back(off_AB);

      CD_blocks_packed_.push_back(offset_r + nx[k+1]);
      CD_blocks_packed_.push_back(offset_c);
      CD_blocks_packed_.push_back(ng[k]);
      CD_blocks_packed_.push_back(nx[k] + nu[k]);
      off_CD += ng[k] * (nx[k] + nu[k]);
      CD_offsets_.push_back(off_CD);

      offset_c += nx[k] + nu[k];
      offset_r += nx[k+1] + ng[k];
    }
    // Terminal CD block
    CD_blocks_packed_.push_back(offset_r);
    CD_blocks_packed_.push_back(offset_c);
    CD_blocks_packed_.push_back(ng[N_]);
    CD_blocks_packed_.push_back(nx[N_]);
    off_CD += ng[N_] * nx[N_];
    CD_offsets_.push_back(off_CD);  // sentinel: total nnz; CD_offsets length is N+2

    casadi_int off_RSQ = 0;
    RSQ_offsets_.push_back(0);
    casadi_int rsq_off = 0;
    for (casadi_int k = 0; k <= N_; ++k) {
      RSQ_blocks_packed_.push_back(rsq_off);
      RSQ_blocks_packed_.push_back(rsq_off);
      RSQ_blocks_packed_.push_back(nx[k] + nu[k]);
      RSQ_blocks_packed_.push_back(nx[k] + nu[k]);
      off_RSQ += (nx[k] + nu[k]) * (nx[k] + nu[k]);
      RSQ_offsets_.push_back(off_RSQ);
      rsq_off += nx[k] + nu[k];
    }

    ABsp_ = conic_blocksparsity(na_, nx_, AB_blocks_packed_, false);
    CDsp_ = conic_blocksparsity(na_, nx_, CD_blocks_packed_, false);
    RSQsp_ = conic_blocksparsity(nx_, nx_, RSQ_blocks_packed_, false);

    // Build condensed dimensions
    nx_hat_.assign(N_hat_ + 1, 0);
    nu_hat_.assign(N_hat_ + 1, 0);
    ng_hat_.assign(N_hat_ + 1, 0);
    nx_max_ = 0;
    nu_max_block_ = 0;
    nxu_max_block_ = 0;
    for (casadi_int K = 0; K <= N_hat_; ++K) {
      casadi_int k_a = M_[K];
      nx_hat_[K] = nx[k_a];
      if (nx[k_a] > nx_max_) nx_max_ = nx[k_a];
    }
    for (casadi_int K = 0; K < N_hat_; ++K) {
      casadi_int k_a = M_[K], k_b = M_[K+1];
      casadi_int sum_nu = 0, sum_ng = 0, sum_lift = 0;
      for (casadi_int j = 0; j < k_b - k_a; ++j) {
        sum_nu += nu[k_a + j];
        sum_ng += ng[k_a + j];
        if (j > 0) sum_lift += nx[k_a + j];
      }
      nu_hat_[K] = sum_nu;
      ng_hat_[K] = sum_ng + sum_lift;
      if (sum_nu > nu_max_block_) nu_max_block_ = sum_nu;
      casadi_int nxu = nx_hat_[K] + sum_nu;
      if (nxu > nxu_max_block_) nxu_max_block_ = nxu;
    }
    nu_hat_[N_hat_] = 0;
    ng_hat_[N_hat_] = ng[N_];
    casadi_int nxu_term = nx_hat_[N_hat_];
    if (nxu_term > nxu_max_block_) nxu_max_block_ = nxu_term;

    // Build condensed-side block descriptors and offsets
    AB_hat_blocks_packed_.clear();
    CD_hat_blocks_packed_.clear();
    RSQ_hat_blocks_packed_.clear();
    AB_hat_offsets_.assign(N_hat_ + 1, 0);
    CD_hat_offsets_.assign(N_hat_ + 1, 0);
    RSQ_hat_offsets_.assign(N_hat_ + 1, 0);
    casadi_int o_AB = 0, o_CD = 0, o_RSQ = 0;
    for (casadi_int K = 0; K < N_hat_; ++K) {
      casadi_int nxu = nx_hat_[K] + nu_hat_[K];
      AB_hat_blocks_packed_.push_back(0);
      AB_hat_blocks_packed_.push_back(0);
      AB_hat_blocks_packed_.push_back(nx_hat_[K+1]);
      AB_hat_blocks_packed_.push_back(nxu);
      AB_hat_offsets_[K] = o_AB;
      o_AB += nx_hat_[K+1] * nxu;
      CD_hat_blocks_packed_.push_back(0);
      CD_hat_blocks_packed_.push_back(0);
      CD_hat_blocks_packed_.push_back(ng_hat_[K]);
      CD_hat_blocks_packed_.push_back(nxu);
      CD_hat_offsets_[K] = o_CD;
      o_CD += ng_hat_[K] * nxu;
      RSQ_hat_blocks_packed_.push_back(0);
      RSQ_hat_blocks_packed_.push_back(0);
      RSQ_hat_blocks_packed_.push_back(nxu);
      RSQ_hat_blocks_packed_.push_back(nxu);
      RSQ_hat_offsets_[K] = o_RSQ;
      o_RSQ += nxu * nxu;
    }
    AB_hat_offsets_[N_hat_] = o_AB;  // sentinel
    CD_hat_blocks_packed_.push_back(0);
    CD_hat_blocks_packed_.push_back(0);
    CD_hat_blocks_packed_.push_back(ng_hat_[N_hat_]);
    CD_hat_blocks_packed_.push_back(nx_hat_[N_hat_]);
    CD_hat_offsets_[N_hat_] = o_CD;
    o_CD += ng_hat_[N_hat_] * nx_hat_[N_hat_];
    RSQ_hat_blocks_packed_.push_back(0);
    RSQ_hat_blocks_packed_.push_back(0);
    RSQ_hat_blocks_packed_.push_back(nx_hat_[N_hat_]);
    RSQ_hat_blocks_packed_.push_back(nx_hat_[N_hat_]);
    RSQ_hat_offsets_[N_hat_] = o_RSQ;

    // Total condensed dimensions
    nx_total_hat_ = 0;
    for (casadi_int K = 0; K <= N_hat_; ++K) nx_total_hat_ += nx_hat_[K] + nu_hat_[K];
    na_total_hat_ = 0;
    for (casadi_int K = 0; K < N_hat_; ++K) na_total_hat_ += nx_hat_[K+1];
    for (casadi_int K = 0; K <= N_hat_; ++K) na_total_hat_ += ng_hat_[K];

    // Build condensed H sparsity (block-diagonal of dense per-block)
    {
      DM h(nx_total_hat_, nx_total_hat_);
      casadi_int off = 0;
      for (casadi_int K = 0; K <= N_hat_; ++K) {
        casadi_int nxu = nx_hat_[K] + nu_hat_[K];
        h(range(off, off + nxu), range(off, off + nxu)) = DM::ones(nxu, nxu);
        off += nxu;
      }
      H_hat_sp_ = h.sparsity();
    }
    // Build condensed A sparsity (gap rows + path rows interleaved per K)
    {
      DM a(na_total_hat_, nx_total_hat_);
      casadi_int row_off = 0, col_off = 0;
      for (casadi_int K = 0; K <= N_hat_; ++K) {
        casadi_int nxu = nx_hat_[K] + nu_hat_[K];
        if (K < N_hat_) {
          casadi_int nxp1 = nx_hat_[K+1];
          a(range(row_off, row_off + nxp1), range(col_off, col_off + nxu))
            = DM::ones(nxp1, nxu);
          a(range(row_off, row_off + nxp1),
            range(col_off + nxu, col_off + nxu + nxp1))
            = -DM::eye(nxp1);
          row_off += nxp1;
        }
        if (ng_hat_[K] > 0) {
          a(range(row_off, row_off + ng_hat_[K]), range(col_off, col_off + nxu))
            = DM::ones(ng_hat_[K], nxu);
          row_off += ng_hat_[K];
        }
        col_off += nxu;
      }
      A_hat_sp_ = a.sparsity();
    }

    finalize_condense_prob();
  }

  // Pick M_ from condense_partition_strategy_ + condensed_block_count_.
  // Dispatches to frison_uniform_partition() or dp_optimal_partition().
  void Conic::derive_condense_partition() {
    // Detect uniformity of nx and nu (we need uniform nx for the
    // closed-form path; nu can vary but is typically uniform too).
    bool nx_uniform = true;
    for (casadi_int k = 1; k <= N_; ++k) {
      if (nxs_[k] != nxs_[0]) { nx_uniform = false; break; }
    }

    std::string strat = condense_partition_strategy_;
    if (strat == "auto") {
      if (nx_uniform) {
        strat = "uniform";
      } else {
        strat = "optimal";
        if (N_ > 500) {
          casadi_warning("condense_partition_strategy='auto' falling back "
            "to 'optimal' DP at O(K_max*N^2) for N=" + str(N_) + " with "
            "non-uniform nx; this runs once at init but may be slow.  "
            "Set condense_partition_strategy='optimal' to silence this "
            "warning, or supply condense_partition explicitly.");
        }
      }
      condense_partition_strategy_ = strat;  // record actual choice
    }

    if (strat == "uniform") {
      casadi_assert(nx_uniform || condensed_block_count_ > 0,
        "condense_partition_strategy='uniform' requires uniform nx (or "
        "an explicit condensed_block_count); got non-uniform nx and "
        "condensed_block_count=0.");
      frison_uniform_partition();
    } else if (strat == "optimal") {
      dp_optimal_partition();
    } else {
      casadi_error("Unknown condense_partition_strategy '" + strat + "'.");
    }
  }

  // Uniform-stage closed form: when condensed_block_count_ == 0, pick
  // n_hat ~= round(sqrt(N * nu / (nx + nu))) per Frison's analysis;
  // when > 0, use it directly.  Then emit M = uniform split.
  void Conic::frison_uniform_partition() {
    casadi_int n_hat = condensed_block_count_;
    if (n_hat == 0) {
      // Frison's closed-form approximate optimum for uniform stages.
      // Here nu/(nx+nu) is treated as average across stages so the
      // formula degrades gracefully if nu varies a bit.
      double nx_avg = 0, nu_avg = 0;
      for (casadi_int k = 0; k <= N_; ++k) nx_avg += nxs_[k];
      for (casadi_int k = 0; k < N_; ++k) nu_avg += nus_[k];
      nx_avg /= (N_ + 1);
      nu_avg /= std::max<casadi_int>(N_, 1);
      double ratio = nu_avg / std::max(1e-30, nx_avg + nu_avg);
      n_hat = static_cast<casadi_int>(
          std::round(std::sqrt(static_cast<double>(N_) * ratio)));
      n_hat = std::max<casadi_int>(1, std::min<casadi_int>(n_hat, N_));
    }
    casadi_assert(n_hat >= 1 && n_hat <= N_,
      "condensed_block_count must satisfy 1 <= n_hat <= N=" + str(N_)
      + "; got " + str(n_hat) + ".");
    M_.assign(n_hat + 1, 0);
    for (casadi_int K = 0; K <= n_hat; ++K) {
      // floor((K * N) / n_hat) -- equal split rounding down
      M_[K] = (K * N_) / n_hat;
    }
    M_[n_hat] = N_;  // ensure exact endpoint
  }

  // O(K_max * N^2) DP: find M minimizing sum_K block_cost(M_K, M_{K+1})
  // where block_cost(a, b) = (nx[a] + sum_{j=a..b-1} nu[j])^3.
  // If condensed_block_count_ == 0, scans K=1..K_max and picks argmin;
  // else fixes K = condensed_block_count_.
  void Conic::dp_optimal_partition() {
    // Prefix sum of nu for O(1) block-cost lookup.
    std::vector<casadi_int> nu_pref(N_ + 1, 0);
    for (casadi_int k = 0; k < N_; ++k) nu_pref[k+1] = nu_pref[k] + nus_[k];
    auto blk_cost = [&](casadi_int a, casadi_int b) -> double {
      double nxu = static_cast<double>(nxs_[a])
                 + static_cast<double>(nu_pref[b] - nu_pref[a]);
      return nxu * nxu * nxu;
    };

    casadi_int K_max;
    if (condensed_block_count_ > 0) {
      K_max = condensed_block_count_;
    } else {
      // Heuristic cap: ~3x the Frison optimum, clamped to [2, 100, N].
      double nx_avg = 0, nu_avg = 0;
      for (casadi_int k = 0; k <= N_; ++k) nx_avg += nxs_[k];
      for (casadi_int k = 0; k < N_; ++k) nu_avg += nus_[k];
      nx_avg /= (N_ + 1);
      nu_avg /= std::max<casadi_int>(N_, 1);
      double ratio = nu_avg / std::max(1e-30, nx_avg + nu_avg);
      casadi_int est = static_cast<casadi_int>(
          std::ceil(3.0 * std::sqrt(static_cast<double>(N_) * ratio)));
      K_max = std::max<casadi_int>(2,
              std::min<casadi_int>(N_,
              std::min<casadi_int>(100, est)));
    }
    casadi_assert(K_max >= 1 && K_max <= N_,
      "DP K_max out of range [1, N]; got " + str(K_max));

    const double INF = std::numeric_limits<double>::infinity();
    std::vector<std::vector<double>> f(K_max + 1,
        std::vector<double>(N_ + 1, INF));
    std::vector<std::vector<casadi_int>> prev(K_max + 1,
        std::vector<casadi_int>(N_ + 1, -1));
    f[0][0] = 0.0;

    for (casadi_int K = 0; K < K_max; ++K) {
      // Need to keep room for (K_max-K-1) more blocks after this one;
      // last block must end at N.
      casadi_int b_hi_global = N_ - (K_max - K - 1);
      for (casadi_int a = K; a <= b_hi_global - 1; ++a) {
        if (f[K][a] == INF) continue;
        for (casadi_int b = a + 1; b <= b_hi_global; ++b) {
          double c = f[K][a] + blk_cost(a, b);
          if (c < f[K+1][b]) {
            f[K+1][b] = c;
            prev[K+1][b] = a;
          }
        }
      }
    }

    casadi_int n_hat;
    if (condensed_block_count_ > 0) {
      n_hat = condensed_block_count_;
      casadi_assert(f[n_hat][N_] != INF,
        "DP could not place exactly " + str(n_hat) + " blocks over N=" +
        str(N_) + " (infeasible).");
    } else {
      n_hat = 1;
      double best = f[1][N_];
      for (casadi_int K = 2; K <= K_max; ++K) {
        if (f[K][N_] < best) { best = f[K][N_]; n_hat = K; }
      }
    }

    // Backtrack
    M_.assign(n_hat + 1, 0);
    M_[n_hat] = N_;
    for (casadi_int K = n_hat; K > 0; --K) {
      M_[K-1] = prev[K][M_[K]];
    }

    if (debug_) {
      uout() << "Conic DP: K_max=" << K_max << ", chose n_hat=" << n_hat
             << ", cost=" << f[n_hat][N_] << std::endl;
    }
  }

  void Conic::unpack_ocp_blocks(const std::vector<casadi_int>& packed,
                                std::vector<casadi_ocp_block>& dst) {
    dst.resize(packed.size() / 4);
    for (size_t i = 0; i < dst.size(); ++i) {
      dst[i].offset_r = packed[4*i + 0];
      dst[i].offset_c = packed[4*i + 1];
      dst[i].rows     = packed[4*i + 2];
      dst[i].cols     = packed[4*i + 3];
    }
  }

  std::vector<casadi_int> Conic::len_prefixed(casadi_int n,
                                              const std::vector<casadi_int>& v) {
    std::vector<casadi_int> r;
    r.reserve(1 + v.size());
    r.push_back(n);
    r.insert(r.end(), v.begin(), v.end());
    return r;
  }

  // Unpack the *_blocks_packed_ vectors into native casadi_ocp_block
  // vectors and wire p_cond_ to point at our member arrays.  Called
  // from build_condense_blocks() and from the deserializing constructor.
  void Conic::finalize_condense_prob() {
    unpack_ocp_blocks(AB_blocks_packed_, AB_blocks_);
    unpack_ocp_blocks(CD_blocks_packed_, CD_blocks_);
    unpack_ocp_blocks(RSQ_blocks_packed_, RSQ_blocks_);

    p_cond_.nx = get_ptr(nxs_);
    p_cond_.nu = get_ptr(nus_);
    p_cond_.ng = get_ptr(ngs_);
    p_cond_.N = N_;
    p_cond_.AB = get_ptr(AB_blocks_);
    p_cond_.CD = get_ptr(CD_blocks_);
    p_cond_.RSQ = get_ptr(RSQ_blocks_);
    p_cond_.AB_offsets = get_ptr(AB_offsets_);
    p_cond_.CD_offsets = get_ptr(CD_offsets_);
    p_cond_.RSQ_offsets = get_ptr(RSQ_offsets_);
    p_cond_.M = get_ptr(M_);
    p_cond_.N_hat = N_hat_;
    // Per-condensed-stage *_hat dim/offset/block-descriptor arrays now
    // live on casadi_condensing_data and are populated per call by
    // casadi_condensing_set_work.  Setup only fills prob scalars here.
    casadi_condensing_setup(&p_cond_);
  }

  void Conic::qp_codegen_body(CodeGenerator& g) const {
    g.add_auxiliary(CodeGenerator::AUX_QP);
    g.local("d_qp", "struct casadi_qp_data");
    g.local("p_qp", "struct casadi_qp_prob");

    if (condense_) {
      // ---- Condensing pre-block ----
      g.add_auxiliary(CodeGenerator::AUX_CONDENSING);
      g.add_auxiliary(CodeGenerator::AUX_PROJECT);
      g.add_auxiliary(CodeGenerator::AUX_COPY);

      // Read-only block descriptors / dim / offset arrays go to file-scope
      // deduplicated const pool via g.constant() (emits as casadi_sNN[]).
      // Block-pack consts are length-prefixed ([N, off_r, off_c, rows, cols, ...])
      // so casadi_unpack_ocp_blocks can read them directly -- mirrors the
      // fatrop/hpipm pattern (no runtime copy through a static scratch).
      const std::string s_nx          = g.constant(nxs_);
      const std::string s_nu          = g.constant(nus_);
      const std::string s_ng          = g.constant(ngs_);
      const std::string s_M           = g.constant(M_);
      const std::string s_AB_off      = g.constant(AB_offsets_);
      const std::string s_CD_off      = g.constant(CD_offsets_);
      const std::string s_RSQ_off     = g.constant(RSQ_offsets_);
      const std::string s_AB_lp       = g.constant(len_prefixed(N_,     AB_blocks_packed_));
      const std::string s_CD_lp       = g.constant(len_prefixed(N_+1,   CD_blocks_packed_));
      const std::string s_RSQ_lp      = g.constant(len_prefixed(N_+1,   RSQ_blocks_packed_));

      // Prob-side scratch: original (non-hat) block descriptor structs
      // -- target of the unpack from length-prefixed const.  Mirrors
      // fatrop/hpipm codegen_unpack_block.  Does NOT consume iw / w.
      // The *_hat block descriptors and *_hat dim/offset arrays now live
      // on casadi_condensing_data and are populated by
      // casadi_condensing_set_work from the prob recurrence -- no static
      // scratch needed for them.
      g.local("cond_AB_blk[" + str(std::max<casadi_int>(N_, 1)) + "]",
              "static struct casadi_ocp_block");
      g.local("cond_CD_blk[" + str(N_+1) + "]",
              "static struct casadi_ocp_block");
      g.local("cond_RSQ_blk[" + str(N_+1) + "]",
              "static struct casadi_ocp_block");

      g << "casadi_unpack_ocp_blocks(cond_AB_blk, "  << s_AB_lp  << ");\n";
      g << "casadi_unpack_ocp_blocks(cond_CD_blk, "  << s_CD_lp  << ");\n";
      g << "casadi_unpack_ocp_blocks(cond_RSQ_blk, " << s_RSQ_lp << ");\n";

      // Set up condensing prob (scalars; per-stage *_hat lives on data)
      g.local("p_cond", "struct casadi_condensing_prob");
      g << "p_cond.nx = " << s_nx << ";\n";
      g << "p_cond.nu = " << s_nu << ";\n";
      g << "p_cond.ng = " << s_ng << ";\n";
      g << "p_cond.N = " << N_ << ";\n";
      g << "p_cond.AB = cond_AB_blk;\n";
      g << "p_cond.CD = cond_CD_blk;\n";
      g << "p_cond.RSQ = cond_RSQ_blk;\n";
      g << "p_cond.AB_offsets = " << s_AB_off << ";\n";
      g << "p_cond.CD_offsets = " << s_CD_off << ";\n";
      g << "p_cond.RSQ_offsets = " << s_RSQ_off << ";\n";
      g << "p_cond.M = " << s_M << ";\n";
      g << "p_cond.N_hat = " << N_hat_ << ";\n";
      g << "casadi_condensing_setup(&p_cond);\n";

      // d_cond + ALL workspace (per-condensed-stage *_hat dim/offset/
      // block descriptors via iw, per-stage flat in/out via w, internal
      // scratch, condensed CSC, bounds, primal/dual) is claimed and
      // populated by casadi_condensing_set_work in this single call.
      g.local("d_cond", "struct casadi_condensing_data");
      g << "d_cond.prob = &p_cond;\n";
      g << "casadi_condensing_set_work(&d_cond, &arg, &res, &iw, &w);\n";

      // Wire user inputs in arg[] order (arg[0] H, [1] g, [2] A,
      // [3] lba, [4] uba, [5] lbx, [6] ubx).  H, A go through
      // casadi_project; the rest are pointer snapshots read by eval.
      g << g.project("arg[" + str(CONIC_H) + "]", H_,
                     "d_cond.RSQ_val", RSQsp_, "w") << "\n";
      g << "d_cond.g_orig   = arg[" << CONIC_G << "];\n";
      g << g.project("arg[" + str(CONIC_A) + "]", A_,
                     "d_cond.AB_val", ABsp_, "w") << "\n";
      if (CDsp_.nnz() > 0) {
        g << g.project("arg[" + str(CONIC_A) + "]", A_,
                       "d_cond.CD_val", CDsp_, "w") << "\n";
      }
      g << "d_cond.lba_orig = arg[" << CONIC_LBA << "];\n";
      g << "d_cond.uba_orig = arg[" << CONIC_UBA << "];\n";
      g << "d_cond.lbx_orig = arg[" << CONIC_LBX << "];\n";
      g << "d_cond.ubx_orig = arg[" << CONIC_UBX << "];\n";
      g << "casadi_condensing_eval(&d_cond);\n";

      // Set up d_qp/p_qp pointing at the CONDENSED problem
      g << "d_qp.prob = &p_qp;\n";
      g << "p_qp.sp_a = " << g.sparsity(A_hat_sp_) << ";\n";
      g << "p_qp.sp_h = " << g.sparsity(H_hat_sp_) << ";\n";
      g << "casadi_qp_setup(&p_qp);\n";
      g << "casadi_qp_set_work(&d_qp, &arg, &res, &iw, &w);\n";

      g << "d_qp.h = d_cond.h_hat_csc;\n";
      g << "d_qp.g = d_cond.qr_hat_val;\n";
      g << "d_qp.a = d_cond.a_hat_csc;\n";
      g << "d_qp.lbx = d_cond.lbx;\n";
      g << "d_qp.ubx = d_cond.ubx;\n";
      g << "d_qp.lba = d_cond.lba;\n";
      g << "d_qp.uba = d_cond.uba;\n";
      g << "d_qp.x0 = 0;\n";
      g << "d_qp.lam_x0 = 0;\n";
      g << "d_qp.lam_a0 = 0;\n";

      g << "d_qp.f = res[" << CONIC_COST << "];\n";
      g << "d_qp.x = d_cond.x;\n";
      g << "d_qp.lam_x = d_cond.lam_x;\n";
      g << "d_qp.lam_a = d_cond.lam_a;\n";
      return;
    }

    // ---- Standard (non-condensing) path ----
    g << "d_qp.prob = &p_qp;\n";
    g << "p_qp.sp_a = " << g.sparsity(A_) << ";\n";
    g << "p_qp.sp_h = " << g.sparsity(H_) << ";\n";
    g << "casadi_qp_setup(&p_qp);\n";
    g << "casadi_qp_set_work(&d_qp, &arg, &res, &iw, &w);\n";


    g << "d_qp.h = arg[" << CONIC_H << "];\n";
    g << "d_qp.g = arg[" << CONIC_G << "];\n";
    g << "d_qp.a = arg[" << CONIC_A << "];\n";
    g << "d_qp.lbx = arg[" << CONIC_LBX << "];\n";
    g << "d_qp.ubx = arg[" << CONIC_UBX << "];\n";
    g << "d_qp.lba = arg[" << CONIC_LBA << "];\n";
    g << "d_qp.uba = arg[" << CONIC_UBA << "];\n";
    g << "d_qp.x0 = arg[" << CONIC_X0 << "];\n";
    g << "d_qp.lam_x0 = arg[" << CONIC_LAM_X0 << "];\n";
    g << "d_qp.lam_a0 = arg[" << CONIC_LAM_A0 << "];\n";

    g << "d_qp.f = res[" << CONIC_COST << "];\n";
    g << "d_qp.x = res[" << CONIC_X << "];\n";
    g << "d_qp.lam_x = res[" << CONIC_LAM_X << "];\n";
    g << "d_qp.lam_a = res[" << CONIC_LAM_A << "];\n";
  }

  void Conic::qp_codegen_post(CodeGenerator& g) const {
    if (!condense_) return;
    // Lift condensed primal/dual into d_cond.{x,lam_x,lam_a}_lifted,
    // then copy into the user-facing res[CONIC_*] slots.
    g.add_auxiliary(CodeGenerator::AUX_COPY);
    g << "casadi_condensing_lift(&d_cond);\n";
    g << "casadi_copy(d_cond.x_lifted, "     << p_cond_.total_qr
      << ", res[" << CONIC_X << "]);\n";
    g << "casadi_copy(d_cond.lam_x_lifted, " << p_cond_.total_qr
      << ", res[" << CONIC_LAM_X << "]);\n";
    g << "casadi_copy(d_cond.lam_a_lifted, " << (p_cond_.total_b + p_cond_.total_g)
      << ", res[" << CONIC_LAM_A << "]);\n";
  }

} // namespace casadi
