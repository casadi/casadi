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


#include "nlpsol_impl.hpp"
#include "external.hpp"
#include "casadi/core/timing.hpp"
#include "nlp_builder.hpp"

using namespace std;
namespace casadi {

  bool has_nlpsol(const string& name) {
    return Nlpsol::has_plugin(name);
  }

  void load_nlpsol(const string& name) {
    Nlpsol::load_plugin(name);
  }

  string doc_nlpsol(const string& name) {
    return Nlpsol::getPlugin(name).doc;
  }

  Function nlpsol(const string& name, const string& solver,
                  const SXDict& nlp, const Dict& opts) {
    return nlpsol(name, solver, Nlpsol::create_oracle(nlp, opts), opts);
  }

  Function nlpsol(const string& name, const string& solver,
                  const MXDict& nlp, const Dict& opts) {
    return nlpsol(name, solver, Nlpsol::create_oracle(nlp, opts), opts);
  }

  template<typename XType>
  Function Nlpsol::create_oracle(const std::map<std::string, XType>& d,
                                 const Dict& opts) {
    std::vector<XType> nl_in(NL_NUM_IN), nl_out(NL_NUM_OUT);
    for (auto&& i : d) {
      if (i.first=="x") {
        nl_in[NL_X]=i.second;
      } else if (i.first=="p") {
        nl_in[NL_P]=i.second;
      } else if (i.first=="f") {
        nl_out[NL_F]=i.second;
      } else if (i.first=="g") {
        nl_out[NL_G]=i.second;
      } else {
        casadi_error("No such field: " + i.first);
      }
    }
    if (nl_out[NL_F].is_empty()) nl_out[NL_F] = 0;
    if (nl_out[NL_G].is_empty()) nl_out[NL_G] = XType(0, 1);

    // Options for the oracle
    Dict oracle_options;
    Dict::const_iterator it = opts.find("oracle_options");
    if (it!=opts.end()) {
      // "oracle_options" has been set
      oracle_options = it->second;
    } else if ((it=opts.find("verbose")) != opts.end()) {
      // "oracle_options" has not been set, but "verbose" has
      oracle_options["verbose"] = it->second;
    }

    // Create oracle
    return Function("nlp", nl_in, nl_out, NL_INPUTS, NL_OUTPUTS, oracle_options);
  }

  Function nlpsol(const std::string& name, const std::string& solver,
                  const NlpBuilder& nl, const Dict& opts) {
     MXDict nlp;
     nlp["x"] = vertcat(nl.x);
     nlp["f"] = nl.f;
     nlp["g"] = vertcat(nl.g);
     return nlpsol(name, solver, nlp, opts);
  }

  Function nlpsol(const std::string& name, const std::string& solver,
                  const std::string& fname, const Dict& opts) {
    // If fname ends with .c, JIT
    if (fname.size()>2 && fname.compare(fname.size()-2, fname.size(), ".c")==0) {
      Importer compiler(fname, "clang");
      return nlpsol(name, solver, compiler, opts);
    } else {
      return nlpsol(name, solver, external("nlp", fname), opts);
    }
  }

  Function nlpsol(const string& name, const string& solver,
                  const Importer& compiler, const Dict& opts) {
    return nlpsol(name, solver, external("nlp", compiler), opts);
  }

  Function nlpsol(const string& name, const string& solver,
                  const Function& nlp, const Dict& opts) {
    // Make sure that nlp is sound
    if (nlp.has_free()) {
      casadi_error("Cannot create '" + name + "' since " + str(nlp.get_free()) + " are free.");
    }
    return Function::create(Nlpsol::instantiate(name, solver, nlp), opts);
  }

  vector<string> nlpsol_in() {
    vector<string> ret(nlpsol_n_in());
    for (size_t i=0; i<ret.size(); ++i) ret[i]=nlpsol_in(i);
    return ret;
  }

  vector<string> nlpsol_out() {
    vector<string> ret(nlpsol_n_out());
    for (size_t i=0; i<ret.size(); ++i) ret[i]=nlpsol_out(i);
    return ret;
  }

  double nlpsol_default_in(casadi_int ind) {
    switch (ind) {
    case NLPSOL_LBX:
    case NLPSOL_LBG:
      return -std::numeric_limits<double>::infinity();
    case NLPSOL_UBX:
    case NLPSOL_UBG:
      return std::numeric_limits<double>::infinity();
    default:
      return 0;
    }
  }

  std::vector<double> nlpsol_default_in() {
    vector<double> ret(nlpsol_n_in());
    for (size_t i=0; i<ret.size(); ++i) ret[i]=nlpsol_default_in(i);
    return ret;
  }

  string nlpsol_in(casadi_int ind) {
    switch (static_cast<NlpsolInput>(ind)) {
    case NLPSOL_X0:     return "x0";
    case NLPSOL_P:      return "p";
    case NLPSOL_LBX:    return "lbx";
    case NLPSOL_UBX:    return "ubx";
    case NLPSOL_LBG:    return "lbg";
    case NLPSOL_UBG:    return "ubg";
    case NLPSOL_LAM_X0: return "lam_x0";
    case NLPSOL_LAM_G0: return "lam_g0";
    case NLPSOL_NUM_IN: break;
    }
    return string();
  }

  string nlpsol_out(casadi_int ind) {
    switch (static_cast<NlpsolOutput>(ind)) {
    case NLPSOL_X:     return "x";
    case NLPSOL_F:     return "f";
    case NLPSOL_G:     return "g";
    case NLPSOL_LAM_X: return "lam_x";
    case NLPSOL_LAM_G: return "lam_g";
    case NLPSOL_LAM_P: return "lam_p";
    case NLPSOL_NUM_OUT: break;
    }
    return string();
  }

  casadi_int nlpsol_n_in() {
    return NLPSOL_NUM_IN;
  }

  casadi_int nlpsol_n_out() {
    return NLPSOL_NUM_OUT;
  }

  Nlpsol::Nlpsol(const std::string& name, const Function& oracle)
    : OracleFunction(name, oracle) {

    // Set default options
    callback_step_ = 1;
    eval_errors_fatal_ = false;
    warn_initial_bounds_ = false;
    iteration_callback_ignore_errors_ = false;
    print_time_ = true;
    calc_multipliers_ = false;
    bound_consistency_ = true;
    min_lam_ = 0;
    calc_lam_x_ = calc_f_ = calc_g_ = false;
    calc_lam_p_ = true;
    no_nlp_grad_ = false;
    error_on_fail_ = false;
  }

  Nlpsol::~Nlpsol() {
    clear_mem();
  }

  Sparsity Nlpsol::get_sparsity_in(casadi_int i) {
    switch (static_cast<NlpsolInput>(i)) {
    case NLPSOL_X0:
    case NLPSOL_LBX:
    case NLPSOL_UBX:
    case NLPSOL_LAM_X0:
      return get_sparsity_out(NLPSOL_X);
    case NLPSOL_LBG:
    case NLPSOL_UBG:
    case NLPSOL_LAM_G0:
      return get_sparsity_out(NLPSOL_G);
    case NLPSOL_P:
      return oracle_.sparsity_in(NL_P);
    case NLPSOL_NUM_IN: break;
    }
    return Sparsity();
  }

  Sparsity Nlpsol::get_sparsity_out(casadi_int i) {
    switch (static_cast<NlpsolOutput>(i)) {
    case NLPSOL_F:
      return oracle_.sparsity_out(NL_F);
    case NLPSOL_X:
    case NLPSOL_LAM_X:
      return oracle_.sparsity_in(NL_X);
    case NLPSOL_LAM_G:
    case NLPSOL_G:
      return oracle_.sparsity_out(NL_G);
    case NLPSOL_LAM_P:
      return get_sparsity_in(NLPSOL_P);
    case NLPSOL_NUM_OUT: break;
    }
    return Sparsity();
  }

  const Options Nlpsol::options_
  = {{&OracleFunction::options_},
     {{"iteration_callback",
       {OT_FUNCTION,
        "A function that will be called at each iteration with the solver as input. "
        "Check documentation of Callback."}},
      {"iteration_callback_step",
       {OT_INT,
        "Only call the callback function every few iterations."}},
      {"iteration_callback_ignore_errors",
       {OT_BOOL,
        "If set to true, errors thrown by iteration_callback will be ignored."}},
      {"ignore_check_vec",
       {OT_BOOL,
        "If set to true, the input shape of F will not be checked."}},
      {"warn_initial_bounds",
       {OT_BOOL,
        "Warn if the initial guess does not satisfy LBX and UBX"}},
      {"eval_errors_fatal",
       {OT_BOOL,
        "When errors occur during evaluation of f,g,...,"
        "stop the iterations"}},
      {"verbose_init",
       {OT_BOOL,
        "Print out timing information about "
        "the different stages of initialization"}},
      {"discrete",
       {OT_BOOLVECTOR,
        "Indicates which of the variables are discrete, i.e. integer-valued"}},
      {"calc_multipliers",
      {OT_BOOL,
       "Calculate Lagrange multipliers in the Nlpsol base class"}},
      {"calc_lam_x",
       {OT_BOOL,
        "Calculate 'lam_x' in the Nlpsol base class"}},
      {"calc_lam_p",
       {OT_BOOL,
        "Calculate 'lam_p' in the Nlpsol base class"}},
      {"calc_f",
       {OT_BOOL,
        "Calculate 'f' in the Nlpsol base class"}},
      {"calc_g",
       {OT_BOOL,
        "Calculate 'g' in the Nlpsol base class"}},
      {"no_nlp_grad",
       {OT_BOOL,
        "Prevent the creation of the 'nlp_grad' function"}},
      {"bound_consistency",
       {OT_BOOL,
        "Ensure that primal-dual solution is consistent with the bounds"}},
      {"min_lam",
       {OT_DOUBLE,
        "Minimum allowed multiplier value"}},
      {"oracle_options",
       {OT_DICT,
        "Options to be passed to the oracle function"}},
      {"error_on_fail",
       {OT_BOOL,
        "When the numerical process returns unsuccessfully, raise an error (default false)."}}
     }
  };

  void Nlpsol::init(const Dict& opts) {
    // Call the initialization method of the base class
    OracleFunction::init(opts);

    // Read options
    for (auto&& op : opts) {
      if (op.first=="iteration_callback") {
        fcallback_ = op.second;
      } else if (op.first=="iteration_callback_step") {
        callback_step_ = op.second;
      } else if (op.first=="eval_errors_fatal") {
        eval_errors_fatal_ = op.second;
      } else if (op.first=="warn_initial_bounds") {
        warn_initial_bounds_ = op.second;
      } else if (op.first=="iteration_callback_ignore_errors") {
        iteration_callback_ignore_errors_ = op.second;
      } else if (op.first=="discrete") {
        discrete_ = op.second;
      } else if (op.first=="calc_multipliers") {
        calc_multipliers_ = op.second;
      } else if (op.first=="calc_lam_x") {
        calc_lam_x_ = op.second;
      } else if (op.first=="calc_lam_p") {
        calc_lam_p_ = op.second;
      } else if (op.first=="calc_f") {
        calc_f_ = op.second;
      } else if (op.first=="calc_g") {
        calc_g_ = op.second;
      } else if (op.first=="no_nlp_grad") {
        no_nlp_grad_ = op.second;
      } else if (op.first=="bound_consistency") {
        bound_consistency_ = op.second;
      } else if (op.first=="min_lam") {
        min_lam_ = op.second;
      } else if (op.first=="error_on_fail") {
        error_on_fail_ = op.second;
      }
    }

    // Deprecated option
    if (calc_multipliers_) {
      calc_lam_x_ = true;
      calc_lam_p_ = true;
    }

    // Get dimensions
    nx_ = nnz_out(NLPSOL_X);
    np_ = nnz_in(NLPSOL_P);
    ng_ = nnz_out(NLPSOL_G);

    // No need to calculate non-existant quantities
    if (np_==0) calc_lam_p_ = false;
    if (ng_==0) calc_g_ = false;

    // Consistency check
    if (no_nlp_grad_) {
      casadi_assert(!calc_lam_p_, "Options 'no_nlp_grad' and 'calc_lam_p' inconsistent");
      casadi_assert(!calc_lam_x_, "Options 'no_nlp_grad' and 'calc_lam_x' inconsistent");
      casadi_assert(!calc_f_, "Options 'no_nlp_grad' and 'calc_f' inconsistent");
      casadi_assert(!calc_g_, "Options 'no_nlp_grad' and 'calc_g' inconsistent");
    }

    // Dimension checks
    casadi_assert(sparsity_out_.at(NLPSOL_G).is_dense()
                          && sparsity_out_.at(NLPSOL_G).is_vector(),
        "Expected a dense vector 'g', but got " + sparsity_out_.at(NLPSOL_G).dim(true) + ".");

    casadi_assert(sparsity_out_.at(NLPSOL_F).is_dense(),
        "Expected a dense 'f', but got " + sparsity_out_.at(NLPSOL_F).dim(true) + ".");

    casadi_assert(sparsity_out_.at(NLPSOL_X).is_dense()
                          && sparsity_out_.at(NLPSOL_X).is_vector(),
      "Expected a dense vector 'x', but got " + sparsity_out_.at(NLPSOL_X).dim(true) + ".");

    // Discrete marker
    mi_ = false;
    if (!discrete_.empty()) {
      casadi_assert(discrete_.size()==nx_, "\"discrete\" option has wrong length");
      if (std::find(discrete_.begin(), discrete_.end(), true)!=discrete_.end()) {
        casadi_assert(integer_support(),
                              "Discrete variables require a solver with integer support");
        mi_ = true;
      }
    }

    set_nlpsol_prob();

    // Allocate memory
    casadi_int sz_w, sz_iw;
    casadi_nlpsol_work(&p_nlp_, &sz_iw, &sz_w);
    alloc_iw(sz_iw, true);
    alloc_w(sz_w, true);

    if (!fcallback_.is_null()) {
      // Consistency checks
      casadi_assert_dev(!fcallback_.is_null());
      casadi_assert(fcallback_.n_out()==1 && fcallback_.numel_out()==1,
        "Callback function must return a scalar.");
      casadi_assert(fcallback_.n_in()==n_out_,
        "Callback input signature must match the NLP solver output signature");
      for (casadi_int i=0; i<n_out_; ++i) {
        casadi_assert(fcallback_.size_in(i)==size_out(i),
          "Callback function input size mismatch. For argument '" + nlpsol_out(i) + "', "
          "callback has shape " + fcallback_.sparsity_in(i).dim() + " while NLP has " +
          sparsity_out_.at(i).dim() + ".");
        // TODO(@jaeandersson): Wrap fcallback_ in a function with correct sparsity
        casadi_assert(fcallback_.sparsity_in(i)==sparsity_out_.at(i),
          "Callback function input size mismatch. "
          "For argument " + nlpsol_out(i) + "', callback has shape " +
          fcallback_.sparsity_in(i).dim() + " while NLP has " +
          sparsity_out_.at(i).dim() + ".");
      }

      // Allocate temporary memory
      alloc(fcallback_);
    }

    // Function calculating f, g and the gradient of the Lagrangian w.r.t. x and p
    if (!no_nlp_grad_) {
      create_function("nlp_grad", {"x", "p", "lam:f", "lam:g"},
                      {"f", "g", "grad:gamma:x", "grad:gamma:p"},
                      {{"gamma", {"f", "g"}}});
    }
  }

  void Nlpsol::set_nlpsol_prob() {
    p_nlp_.nx = nx_;
    p_nlp_.ng = ng_;
    p_nlp_.np = np_;
  }

  int Nlpsol::init_mem(void* mem) const {
    if (OracleFunction::init_mem(mem)) return 1;
    auto m = static_cast<NlpsolMemory*>(mem);
    m->add_stat("callback_fun");
    m->success = false;
    m->unified_return_status = SOLVER_RET_UNKNOWN;
    return 0;
  }

  void Nlpsol::check_inputs(void* mem) const {
    auto m = static_cast<NlpsolMemory*>(mem);
    auto d_nlp = &m->d_nlp;

    // Skip check?
    if (!inputs_check_) return;

    const double inf = std::numeric_limits<double>::infinity();

    // Number of equality constraints
    casadi_int n_eq = 0;

    // Detect ill-posed problems (simple bounds)
    for (casadi_int i=0; i<nx_; ++i) {
      double lb = d_nlp->lbz[i];
      double ub = d_nlp->ubz[i];
      double x0 = d_nlp->z[i];
      casadi_assert(lb <= ub && lb!=inf && ub!=-inf,
          "Ill-posed problem detected: "
          "LBX[" + str(i) + "] <= UBX[" + str(i) + "] was violated. "
          "Got LBX[" + str(i) + "]=" + str(lb) + " and UBX[" + str(i) + "] = " + str(ub) + ".");
      if (warn_initial_bounds_ && (x0>ub || x0<lb)) {
        casadi_warning("Nlpsol: The initial guess does not satisfy LBX and UBX. "
          "Option 'warn_initial_bounds' controls this warning.");
        break;
      }
      if (lb==ub) n_eq++;
    }

    // Detect ill-posed problems (nonlinear bounds)
    for (casadi_int i=0; i<ng_; ++i) {
      double lb = d_nlp->lbz[nx_ + i];
      double ub = d_nlp->ubz[nx_ + i];
      casadi_assert(lb <= ub && lb!=inf && ub!=-inf,
        "Ill-posed problem detected: "
        "LBG[" + str(i) + "] <= UBG[" + str(i) + "] was violated. "
        "Got LBG[" + str(i) + "] = " + str(lb) + " and UBG[" + str(i) + "] = " + str(ub) + ".");
      if (lb==ub) n_eq++;
    }

    // Make sure enough degrees of freedom
    using casadi::str; // Workaround, MingGW bug, cf. CasADi issue #890
    if (n_eq> nx_) {
      casadi_warning("NLP is overconstrained: There are " + str(n_eq) +
      " equality constraints but only " + str(nx_) + " variables.");
    }
  }

  std::map<std::string, Nlpsol::Plugin> Nlpsol::solvers_;

  const std::string Nlpsol::infix_ = "nlpsol";

  DM Nlpsol::getReducedHessian() {
    casadi_error("getReducedHessian not defined for class " + class_name());
    return DM();
  }

  void Nlpsol::setOptionsFromFile(const std::string & file) {
    casadi_error("setOptionsFromFile not defined for class " + class_name());
  }

  void Nlpsol::bound_consistency(casadi_int n, double* z, double* lam,
                                 const double* lbz, const double* ubz) {
    casadi_assert_dev(z!=nullptr);
    casadi_assert_dev(lam!=nullptr);
    casadi_assert_dev(lbz!=nullptr);
    casadi_assert_dev(ubz!=nullptr);
    // Local variables
    casadi_int i;
    // Loop over variables
    for (i=0; i<n; ++i) {
      // Make sure bounds are respected
      z[i] = std::fmin(std::fmax(z[i], lbz[i]), ubz[i]);
      // Adjust multipliers
      if (std::isinf(lbz[i]) && std::isinf(ubz[i])) {
        // Both multipliers are infinite
        lam[i] = 0.;
      } else if (std::isinf(lbz[i]) || z[i] - lbz[i] > ubz[i] - z[i]) {
        // Infinite lower bound or closer to upper bound than lower bound
        lam[i] = std::fmax(0., lam[i]);
      } else if (std::isinf(ubz[i]) || z[i] - lbz[i] < ubz[i] - z[i]) {
        // Infinite upper bound or closer to lower bound than upper bound
        lam[i] = std::fmin(0., lam[i]);
      }
    }
  }

  int Nlpsol::eval(const double** arg, double** res, casadi_int* iw, double* w, void* mem) const {
    auto m = static_cast<NlpsolMemory*>(mem);

    auto d_nlp = &m->d_nlp;

    // Bounds, given parameter values
    d_nlp->p = arg[NLPSOL_P];
    const double *lbx = arg[NLPSOL_LBX];
    const double *ubx = arg[NLPSOL_UBX];
    const double *lbg = arg[NLPSOL_LBG];
    const double *ubg = arg[NLPSOL_UBG];

    // Get input pointers
    const double *x0 = arg[NLPSOL_X0];
    const double *lam_x0 = arg[NLPSOL_LAM_X0];
    const double *lam_g0 = arg[NLPSOL_LAM_G0];
    arg += NLPSOL_NUM_IN;

    // Get output pointers
    double *x = res[NLPSOL_X];
    double *f = res[NLPSOL_F];
    double *g = res[NLPSOL_G];
    double *lam_x = res[NLPSOL_LAM_X];
    double *lam_g = res[NLPSOL_LAM_G];
    double *lam_p = res[NLPSOL_LAM_P];
    res += NLPSOL_NUM_OUT;

    // Reset the solver, prepare for solution
    setup(m, arg, res, iw, w);

    // Set initial guess
    casadi_copy(x0, nx_, d_nlp->z);
    casadi_copy(lam_x0, nx_, d_nlp->lam);
    casadi_copy(lam_g0, ng_, d_nlp->lam + nx_);

    // Set multipliers to nan
    casadi_fill(d_nlp->lam_p, np_, nan);

    // Reset f, g
    d_nlp->f = nan;
    casadi_fill(d_nlp->z + nx_, ng_, nan);

    // Get bounds
    casadi_copy(lbx, nx_, d_nlp->lbz);
    casadi_copy(lbg, ng_, d_nlp->lbz + nx_);
    casadi_copy(ubx, nx_, d_nlp->ubz);
    casadi_copy(ubg, ng_, d_nlp->ubz + nx_);

    // Check the provided inputs
    check_inputs(m);

    // Solve the NLP
    int flag = solve(m);

    // Calculate multiplers
    if ((calc_f_ || calc_g_ || calc_lam_x_ || calc_lam_p_) && !flag) {
      const double lam_f = 1.;
      m->arg[0] = d_nlp->z;
      m->arg[1] = d_nlp->p;
      m->arg[2] = &lam_f;
      m->arg[3] = d_nlp->lam + nx_;
      m->res[0] = calc_f_ ? &d_nlp->f : nullptr;
      m->res[1] = calc_g_ ? d_nlp->z + nx_ : nullptr;
      m->res[2] = calc_lam_x_ ? d_nlp->lam : nullptr;
      m->res[3] = calc_lam_p_ ? d_nlp->lam_p : nullptr;
      if (calc_function(m, "nlp_grad")) {
        casadi_warning("Failed to calculate multipliers");
      }
      if (calc_lam_x_) casadi_scal(nx_, -1., d_nlp->lam);
      if (calc_lam_p_) casadi_scal(np_, -1., d_nlp->lam_p);
    }

    // Make sure that an optimal solution is consistant with bounds
    if (bound_consistency_ && !flag) {
      bound_consistency(nx_+ng_, d_nlp->z, d_nlp->lam, d_nlp->lbz, d_nlp->ubz);
    }

    // Get optimal solution
    casadi_copy(d_nlp->z, nx_, x);
    casadi_copy(d_nlp->z + nx_, ng_, g);
    casadi_copy(d_nlp->lam, nx_, lam_x);
    casadi_copy(d_nlp->lam + nx_, ng_, lam_g);
    casadi_copy(d_nlp->lam_p, np_, lam_p);
    casadi_copy(&d_nlp->f, 1, f);

    if (error_on_fail_ && !m->success)
      casadi_error("nlpsol process failed. "
                   "Set 'error_on_fail' option to false to ignore this error.");
    return flag;
  }

  void Nlpsol::set_work(void* mem, const double**& arg, double**& res,
                        casadi_int*& iw, double*& w) const {
    auto m = static_cast<NlpsolMemory*>(mem);

    // Problem has not been solved at this point
    m->success = false;
    m->unified_return_status = SOLVER_RET_UNKNOWN;

    m->d_nlp.prob = &p_nlp_;

    casadi_nlpsol_init(&m->d_nlp, &iw, &w);
  }

  std::vector<std::string> nlpsol_options(const std::string& name) {
    return Nlpsol::plugin_options(name).all();
  }

  std::string nlpsol_option_type(const std::string& name, const std::string& op) {
    return Nlpsol::plugin_options(name).type(op);
  }

  std::string nlpsol_option_info(const std::string& name, const std::string& op) {
    return Nlpsol::plugin_options(name).info(op);
  }

  void Nlpsol::disp_more(std::ostream& stream) const {
    stream << "minimize f(x;p) subject to lbx<=x<=ubx, lbg<=g(x;p)<=ubg defined by:\n";
    oracle_.disp(stream, true);
  }

  Function Nlpsol::kkt() const {
    // Quick return if cached
    if (kkt_.alive()) {
      return shared_cast<Function>(kkt_.shared());
    }

    // Generate KKT function
    Function ret = oracle_.factory("kkt", {"x", "p", "lam:f", "lam:g"},
      {"jac:g:x", "sym:hess:gamma:x:x"}, {{"gamma", {"f", "g"}}});

    // Cache and return
    kkt_ = ret;
    return ret;
  }


  Function Nlpsol::
  get_forward(casadi_int nfwd, const std::string& name,
              const std::vector<std::string>& inames,
              const std::vector<std::string>& onames,
              const Dict& opts) const {

    // Symbolic expression for the input
    vector<MX> arg = mx_in(), res = mx_out();

    // Initial guesses not used for derivative calculations
    for (NlpsolInput i : {NLPSOL_X0, NLPSOL_LAM_X0, NLPSOL_LAM_G0}) {
      std::string name = arg[i].is_symbolic() ? arg[i].name() : "temp";
      arg[i] = MX::sym(name, Sparsity(arg[i].size()));
    }

    // Optimal solution
    MX x = res[NLPSOL_X];
    MX lam_g = res[NLPSOL_LAM_G];
    MX lam_x = res[NLPSOL_LAM_X];
    MX lam_p = res[NLPSOL_LAM_P];
    MX f = res[NLPSOL_F];
    MX g = res[NLPSOL_G];

    // Inputs used
    MX lbx = arg[NLPSOL_LBX];
    MX ubx = arg[NLPSOL_UBX];
    MX lbg = arg[NLPSOL_LBG];
    MX ubg = arg[NLPSOL_UBG];
    MX p = arg[NLPSOL_P];

    // Get KKT function
    Function kkt = this->kkt();

    // Hessian of the Lagrangian, Jacobian of the constraints
    vector<MX> HJ_res = kkt({x, p, 1, lam_g});
    MX JG = HJ_res.at(0);
    MX HL = HJ_res.at(1);

    // Active set (assumed known and given by the multiplier signs)
    MX ubIx = lam_x > min_lam_;
    MX lbIx = lam_x < -min_lam_;
    MX bIx = ubIx + lbIx;
    MX iIx = 1-bIx;
    MX ubIg = lam_g > min_lam_;
    MX lbIg = lam_g < -min_lam_;
    MX bIg = ubIg + lbIg;
    MX iIg = 1-bIg;

    // KKT matrix
    MX H_11 = mtimes(diag(iIx), HL) + diag(bIx);
    MX H_12 = mtimes(diag(iIx), JG.T());
    MX H_21 = mtimes(diag(bIg), JG);
    MX H_22 = diag(-iIg);
    MX H = MX::blockcat({{H_11, H_12}, {H_21, H_22}});

    // Sensitivity inputs
    vector<MX> fseed(NLPSOL_NUM_IN);
    MX fwd_lbx = fseed[NLPSOL_LBX] = MX::sym("fwd_lbx", repmat(x.sparsity(), 1, nfwd));
    MX fwd_ubx = fseed[NLPSOL_UBX] = MX::sym("fwd_ubx", repmat(x.sparsity(), 1, nfwd));
    MX fwd_lbg = fseed[NLPSOL_LBG] = MX::sym("fwd_lbg", repmat(g.sparsity(), 1, nfwd));
    MX fwd_ubg = fseed[NLPSOL_UBG] = MX::sym("fwd_ubg", repmat(g.sparsity(), 1, nfwd));
    MX fwd_p = fseed[NLPSOL_P] = MX::sym("fwd_p", repmat(p.sparsity(), 1, nfwd));

    // Guesses are unused
    for (NlpsolInput i : {NLPSOL_X0, NLPSOL_LAM_X0, NLPSOL_LAM_G0}) {
      fseed[i] = MX(repmat(Sparsity(arg[i].size()), 1, nfwd));
    }

    // nlp_grad has the signature
    // (x, p, lam_f, lam_g) -> (f, g, grad_x, grad_p)
    // with lam_f=1 and lam_g=lam_g, grad_x = -lam_x, grad_p=-lam_p
    Function nlp_grad = get_function("nlp_grad");

    // fwd_nlp_grad has the signature
    // (x, p, lam_f, lam_g, f, g, grad_x, grad_p,
    //  fwd_x, fwd_p, fwd_lam_f, fwd_lam_g)
    // -> (fwd_f, fwd_g, fwd_grad_x, fwd_grad_p)
    Function fwd_nlp_grad = nlp_grad.forward(nfwd);

    // Calculate sensitivities from fwd_p
    vector<MX> vv = {x, p, 1, lam_g, f, g, -lam_x, -lam_p, 0., fwd_p, 0., 0.};
    vv = fwd_nlp_grad(vv);
    MX fwd_g_p = vv.at(1);
    MX fwd_gL_p = vv.at(2);

    // Propagate forward seeds
    MX fwd_alpha_x = (if_else(lbIx, fwd_lbx, 0) + if_else(ubIx, fwd_ubx, 0))
                   - if_else(iIx, fwd_gL_p, 0);
    MX fwd_alpha_g = (if_else(ubIg, fwd_ubg, 0) + if_else(lbIg, fwd_lbg, 0))
                   - if_else(bIg, fwd_g_p, 0);
    MX v = MX::vertcat({fwd_alpha_x, fwd_alpha_g});

    // Solve
    v = MX::solve(H, v, "qr");

    // Extract sensitivities in x, lam_x and lam_g
    vector<MX> v_split = vertsplit(v, {0, nx_, nx_+ng_});
    MX fwd_x = v_split.at(0);
    MX fwd_lam_g = v_split.at(1);

    // Calculate sensitivities in lam_x, lam_g
    vv = {x, p, 1, lam_g, f, g, -lam_x, -lam_p,
          fwd_x, fwd_p, 0, fwd_lam_g};
    vv = fwd_nlp_grad(vv);
    MX fwd_f = vv.at(0);
    MX fwd_g = vv.at(1);
    MX fwd_lam_x = -vv.at(2);
    MX fwd_lam_p = -vv.at(3);

    // Forward sensitivities
    vector<MX> fsens(NLPSOL_NUM_OUT);
    fsens[NLPSOL_X] = fwd_x;
    fsens[NLPSOL_F] = fwd_f;
    fsens[NLPSOL_G] = fwd_g;
    fsens[NLPSOL_LAM_X] = fwd_lam_x;
    fsens[NLPSOL_LAM_G] = fwd_lam_g;
    fsens[NLPSOL_LAM_P] = fwd_lam_p;

    // Gather return values
    arg.insert(arg.end(), res.begin(), res.end());
    arg.insert(arg.end(), fseed.begin(), fseed.end());
    res = fsens;

    return Function(name, arg, res, inames, onames, opts);
  }

  Function Nlpsol::
  get_reverse(casadi_int nadj, const std::string& name,
              const std::vector<std::string>& inames,
              const std::vector<std::string>& onames,
              const Dict& opts) const {
    // Symbolic expression for the input
    vector<MX> arg = mx_in(), res = mx_out();

    // Initial guesses not used for derivative calculations
    for (NlpsolInput i : {NLPSOL_X0, NLPSOL_LAM_X0, NLPSOL_LAM_G0}) {
      std::string name = arg[i].is_symbolic() ? arg[i].name() : "temp";
      arg[i] = MX::sym(name, Sparsity(arg[i].size()));
    }

    // Optimal solution
    MX x = res[NLPSOL_X];
    MX lam_g = res[NLPSOL_LAM_G];
    MX lam_x = res[NLPSOL_LAM_X];
    MX lam_p = res[NLPSOL_LAM_P];
    MX f = res[NLPSOL_F];
    MX g = res[NLPSOL_G];

    // Inputs used
    MX lbx = arg[NLPSOL_LBX];
    MX ubx = arg[NLPSOL_UBX];
    MX lbg = arg[NLPSOL_LBG];
    MX ubg = arg[NLPSOL_UBG];
    MX p = arg[NLPSOL_P];

    // Get KKT function
    Function kkt = this->kkt();

    // Hessian of the Lagrangian, Jacobian of the constraints
    vector<MX> HJ_res = kkt({x, p, 1, lam_g});
    MX JG = HJ_res.at(0);
    MX HL = HJ_res.at(1);

    // Active set (assumed known and given by the multiplier signs)
    MX ubIx = lam_x > min_lam_;
    MX lbIx = lam_x < -min_lam_;
    MX bIx = ubIx + lbIx;
    MX iIx = 1-bIx;
    MX ubIg = lam_g > min_lam_;
    MX lbIg = lam_g < -min_lam_;
    MX bIg = ubIg + lbIg;
    MX iIg = 1-bIg;

    // KKT matrix
    MX H_11 = mtimes(diag(iIx), HL) + diag(bIx);
    MX H_12 = mtimes(diag(iIx), JG.T());
    MX H_21 = mtimes(diag(bIg), JG);
    MX H_22 = diag(-iIg);
    MX H = MX::blockcat({{H_11, H_12}, {H_21, H_22}});

    // Sensitivity inputs
    vector<MX> aseed(NLPSOL_NUM_OUT);
    MX adj_x = aseed[NLPSOL_X] = MX::sym("adj_x", repmat(x.sparsity(), 1, nadj));
    MX adj_lam_g = aseed[NLPSOL_LAM_G] = MX::sym("adj_lam_g", repmat(g.sparsity(), 1, nadj));
    MX adj_lam_x = aseed[NLPSOL_LAM_X] = MX::sym("adj_lam_x", repmat(x.sparsity(), 1, nadj));
    MX adj_lam_p = aseed[NLPSOL_LAM_P] = MX::sym("adj_lam_p", repmat(p.sparsity(), 1, nadj));
    MX adj_f = aseed[NLPSOL_F] = MX::sym("adj_f", Sparsity::dense(1, nadj));
    MX adj_g = aseed[NLPSOL_G] = MX::sym("adj_g", repmat(g.sparsity(), 1, nadj));

    // nlp_grad has the signature
    // (x, p, lam_f, lam_g) -> (f, g, grad_x, grad_p)
    // with lam_f=1 and lam_g=lam_g, grad_x = -lam_x, grad_p=-lam_p
    Function nlp_grad = get_function("nlp_grad");

    // rev_nlp_grad has the signature
    // (x, p, lam_f, lam_g, f, g, grad_x, grad_p,
    //  adj_f, adj_g, adj_grad_x, adj_grad_p)
    // -> (adj_x, adj_p, adj_lam_f, adj_lam_g)
    Function rev_nlp_grad = nlp_grad.reverse(nadj);

    // Calculate sensitivities from f, g and lam_x
    vector<MX> vv = {x, p, 1, lam_g, f, g, -lam_x, -lam_p,
                     adj_f, adj_g, -adj_lam_x, -adj_lam_p};
    vv = rev_nlp_grad(vv);
    MX adj_x0 = vv.at(0);
    MX adj_p0 = vv.at(1);
    MX adj_lam_g0 = vv.at(3);

    // Solve to get beta_x_bar, beta_g_bar
    MX v = MX::vertcat({adj_x + adj_x0, adj_lam_g + adj_lam_g0});
    v = MX::solve(H.T(), v, "qr");
    vector<MX> v_split = vertsplit(v, {0, nx_, nx_+ng_});
    MX beta_x_bar = v_split.at(0);
    MX beta_g_bar = v_split.at(1);

    // Calculate sensitivities in p
    vv = {x, p, 1, lam_g, f, g, -lam_x, -lam_p,
          0, bIg*beta_g_bar, iIx*beta_x_bar, 0};
    vv = rev_nlp_grad(vv);
    MX adj_p = vv.at(1);

    // Reverse sensitivities
    vector<MX> asens(NLPSOL_NUM_IN);
    asens[NLPSOL_UBX] = if_else(ubIx, beta_x_bar, 0);
    asens[NLPSOL_LBX] = if_else(lbIx, beta_x_bar, 0);
    asens[NLPSOL_UBG] = if_else(ubIg, beta_g_bar, 0);
    asens[NLPSOL_LBG] = if_else(lbIg, beta_g_bar, 0);
    asens[NLPSOL_P] = adj_p0 - adj_p;

    // Guesses are unused
    for (NlpsolInput i : {NLPSOL_X0, NLPSOL_LAM_X0, NLPSOL_LAM_G0}) {
      asens[i] = MX(repmat(Sparsity(arg[i].size()), 1, nadj));
    }

    // Gather return values
    arg.insert(arg.end(), res.begin(), res.end());
    arg.insert(arg.end(), aseed.begin(), aseed.end());
    res = asens;

    return Function(name, arg, res, inames, onames, opts);
  }

  int Nlpsol::callback(NlpsolMemory* m) const {
    // Quick return if no callback function
    if (fcallback_.is_null()) return 0;
    // Callback inputs
    fill_n(m->arg, fcallback_.n_in(), nullptr);

    auto d_nlp = &m->d_nlp;

    m->arg[NLPSOL_X] = d_nlp->z;
    m->arg[NLPSOL_F] = &d_nlp->f;
    m->arg[NLPSOL_G] = d_nlp->z + nx_;
    m->arg[NLPSOL_LAM_G] = d_nlp->lam + nx_;
    m->arg[NLPSOL_LAM_X] = d_nlp->lam;

    // Callback outputs
    fill_n(m->res, fcallback_.n_out(), nullptr);
    double ret = 0;
    m->res[0] = &ret;

    // Start timer
    m->fstats.at("callback_fun").tic();
    try {
      // Evaluate
      fcallback_(m->arg, m->res, m->iw, m->w, 0);
    } catch(KeyboardInterruptException& ex) {
      throw;
    } catch(exception& ex) {
      print("WARNING: intermediate_callback error: %s\n", ex.what());
      if (!iteration_callback_ignore_errors_) ret=1;
    }

    // User user interruption?
    if (static_cast<casadi_int>(ret)) return 1;

    // Stop timer
    m->fstats.at("callback_fun").toc();

    return 0;
  }

  Dict Nlpsol::get_stats(void* mem) const {
    Dict stats = OracleFunction::get_stats(mem);
    auto m = static_cast<NlpsolMemory*>(mem);
    stats["success"] = m->success;
    stats["unified_return_status"] = string_from_UnifiedReturnStatus(m->unified_return_status);
    return stats;
  }

  void Nlpsol::codegen_body(CodeGenerator& g) const {
    g.local("d_nlp", "struct casadi_nlpsol_data");
    g.local("p_nlp", "struct casadi_nlpsol_prob");

    g << "d_nlp.prob = &p_nlp;\n";
    g << "p_nlp.nx = " << nx_ << ";\n";
    g << "p_nlp.ng = " << ng_ << ";\n";
    g << "p_nlp.np = " << np_ << ";\n";
    g << "casadi_nlpsol_init(&d_nlp, &iw, &w);\n";
  }

  void Nlpsol::serialize_body(SerializingStream &s) const {
    OracleFunction::serialize_body(s);

    s.version("Nlpsol", 1);
    s.pack("Nlpsol::nx", nx_);
    s.pack("Nlpsol::ng", ng_);
    s.pack("Nlpsol::np", np_);
    s.pack("Nlpsol::fcallback", fcallback_);
    s.pack("Nlpsol::callback_step", callback_step_);
    s.pack("Nlpsol::error_on_fail", error_on_fail_);
    s.pack("Nlpsol::eval_errors_fatal", eval_errors_fatal_);
    s.pack("Nlpsol::warn_initial_bounds", warn_initial_bounds_);
    s.pack("Nlpsol::iteration_callback_ignore_errors", iteration_callback_ignore_errors_);
    s.pack("Nlpsol::calc_multipliers", calc_multipliers_);
    s.pack("Nlpsol::calc_lam_x", calc_lam_x_);
    s.pack("Nlpsol::calc_lam_p", calc_lam_p_);
    s.pack("Nlpsol::calc_f", calc_f_);
    s.pack("Nlpsol::calc_g", calc_g_);
    s.pack("Nlpsol::min_lam", min_lam_);
    s.pack("Nlpsol::bound_consistency", bound_consistency_);
    s.pack("Nlpsol::no_nlp_grad", no_nlp_grad_);
    s.pack("Nlpsol::discrete", discrete_);
    s.pack("Nlpsol::mi", mi_);
  }

  void Nlpsol::serialize_type(SerializingStream &s) const {
    OracleFunction::serialize_type(s);
    PluginInterface<Nlpsol>::serialize_type(s);
  }

  ProtoFunction* Nlpsol::deserialize(DeserializingStream& s) {
    return PluginInterface<Nlpsol>::deserialize(s);
  }

  Nlpsol::Nlpsol(DeserializingStream & s) : OracleFunction(s) {
    s.version("Nlpsol", 1);
    s.unpack("Nlpsol::nx", nx_);
    s.unpack("Nlpsol::ng", ng_);
    s.unpack("Nlpsol::np", np_);
    s.unpack("Nlpsol::fcallback", fcallback_);
    s.unpack("Nlpsol::callback_step", callback_step_);
    s.unpack("Nlpsol::error_on_fail", error_on_fail_);
    s.unpack("Nlpsol::eval_errors_fatal", eval_errors_fatal_);
    s.unpack("Nlpsol::warn_initial_bounds", warn_initial_bounds_);
    s.unpack("Nlpsol::iteration_callback_ignore_errors", iteration_callback_ignore_errors_);
    s.unpack("Nlpsol::calc_multipliers", calc_multipliers_);
    s.unpack("Nlpsol::calc_lam_x", calc_lam_x_);
    s.unpack("Nlpsol::calc_lam_p", calc_lam_p_);
    s.unpack("Nlpsol::calc_f", calc_f_);
    s.unpack("Nlpsol::calc_g", calc_g_);
    s.unpack("Nlpsol::min_lam", min_lam_);
    s.unpack("Nlpsol::bound_consistency", bound_consistency_);
    s.unpack("Nlpsol::no_nlp_grad", no_nlp_grad_);
    s.unpack("Nlpsol::discrete", discrete_);
    s.unpack("Nlpsol::mi", mi_);
    set_nlpsol_prob();
  }

} // namespace casadi
