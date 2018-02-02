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

  Options Nlpsol::options_
  = {{&OracleFunction::options_},
     {{"expand",
       {OT_BOOL,
        "Replace MX with SX expressions in problem formulation [false]"}},
      {"iteration_callback",
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
      {"oracle_options",
      {OT_DICT,
       "Options to be passed to the oracle function"}}
     }
  };

  void Nlpsol::init(const Dict& opts) {
    // Call the initialization method of the base class
    OracleFunction::init(opts);

    // Default options
    bool expand = false;

    // Read options
    for (auto&& op : opts) {
      if (op.first=="expand") {
        expand = op.second;
      } else if (op.first=="iteration_callback") {
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
      }
    }

    // Replace MX oracle with SX oracle?
    if (expand) oracle_ = oracle_.expand();

    // Get dimensions
    nx_ = nnz_out(NLPSOL_X);
    np_ = nnz_in(NLPSOL_P);
    ng_ = nnz_out(NLPSOL_G);

    // Dimension checks
    casadi_assert(sparsity_out_.at(NLPSOL_G).is_dense()
                          && sparsity_out_.at(NLPSOL_G).is_vector(),
        "Expected a dense vector 'g', but got " + sparsity_out_.at(NLPSOL_G).dim() + ".");

    casadi_assert(sparsity_out_.at(NLPSOL_F).is_dense(),
        "Expected a dense 'f', but got " + sparsity_out_.at(NLPSOL_F).dim() + ".");

    casadi_assert(sparsity_out_.at(NLPSOL_X).is_dense()
                          && sparsity_out_.at(NLPSOL_X).is_vector(),
      "Expected a dense vector 'x', but got " + sparsity_out_.at(NLPSOL_X).dim() + ".");

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

    // Function to calculate multiplers
    if (calc_multipliers_) {
      create_function("nlp_mult", {"x", "p", "lam:f", "lam:g"},
                      {"grad:gamma:x", "grad:gamma:p"}, {{"gamma", {"f", "g"}}});
    }
  }

  int Nlpsol::init_mem(void* mem) const {
    if (OracleFunction::init_mem(mem)) return 1;
    auto m = static_cast<NlpsolMemory*>(mem);
    m->add_stat(name_);
    m->add_stat("callback_fun");
    return 0;
  }

  void Nlpsol::check_inputs(void* mem) const {
    auto m = static_cast<NlpsolMemory*>(mem);

    // Skip check?
    if (!inputs_check_) return;

    const double inf = std::numeric_limits<double>::infinity();

    // Number of equality constraints
    casadi_int n_eq = 0;

    // Detect ill-posed problems (simple bounds)
    for (casadi_int i=0; i<nx_; ++i) {
      double lb = m->lbx ? m->lbx[i] : 0;
      double ub = m->ubx ? m->ubx[i] : 0;
      double x0 = m->x0 ? m->x0[i] : 0;
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
      double lb = m->lbg ? m->lbg[i] : 0;
      double ub = m->ubg ? m->ubg[i] : 0;
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

  int Nlpsol::eval(const double** arg, double** res, casadi_int* iw, double* w, void* mem) const {
    auto m = static_cast<NlpsolMemory*>(mem);

    // Reset statistics
    for (auto&& s : m->fstats) s.second.reset();
    m->fstats.at(name_).tic();

    // Reset the solver, prepare for solution
    setup(m, arg, res, iw, w);

    // Set multipliers to nan
    casadi_fill(m->lam_x, nx_, nan);
    casadi_fill(m->lam_g, ng_, nan);
    casadi_fill(m->lam_p, np_, nan);

    // Solve the NLP
    solve(m);

    // Calculate multiplers
    if (calc_multipliers_) {
      casadi_assert(m->x!=0, "Not implemented");
      casadi_assert(ng_==0 || m->lam_g!=0, "Not implemented");
      double lam_f = 1.;
      m->arg[0] = m->x;
      m->arg[1] = m->p;
      m->arg[2] = &lam_f;
      m->arg[3] = m->lam_g;
      m->res[0] = m->lam_x;
      m->res[1] = m->lam_p;
      if (calc_function(m, "nlp_mult")) {
        casadi_warning("Failed to calculate multipliers");
      }

      if (m->lam_x) {
        casadi_scal(nx_, -1., m->lam_x);
        for (casadi_int i=0; i<nx_; ++i) {
          if (m->lam_x[i]>0) {
            // If upper bound isn't active, multiplier is zero
            if (m->x[i] < (m->ubx ? m->ubx[i] : 0)) m->lam_x[i] = 0;
          } else if (m->lam_x[i]<0) {
            // If lower bound isn't active, multiplier is zero
            if (m->x[i] > (m->lbx ? m->lbx[i] : 0)) m->lam_x[i] = 0;
          }
        }
      }
      if (m->lam_p) {
        casadi_scal(np_, -1., m->lam_p);
      }
    }

    // Finalize/print statistics
    m->fstats.at(name_).toc();
    if (print_time_)  print_fstats(m);
    return 0;
  }

  void Nlpsol::set_work(void* mem, const double**& arg, double**& res,
                        casadi_int*& iw, double*& w) const {
    auto m = static_cast<NlpsolMemory*>(mem);

    // Get input pointers
    m->x0 = arg[NLPSOL_X0];
    m->p = arg[NLPSOL_P];
    m->lbx = arg[NLPSOL_LBX];
    m->ubx = arg[NLPSOL_UBX];
    m->lbg = arg[NLPSOL_LBG];
    m->ubg = arg[NLPSOL_UBG];
    m->lam_x0 = arg[NLPSOL_LAM_X0];
    m->lam_g0 = arg[NLPSOL_LAM_G0];
    arg += NLPSOL_NUM_IN;

    // Get output pointers
    m->x = res[NLPSOL_X];
    m->f = res[NLPSOL_F];
    m->g = res[NLPSOL_G];
    m->lam_x = res[NLPSOL_LAM_X];
    m->lam_g = res[NLPSOL_LAM_G];
    m->lam_p = res[NLPSOL_LAM_P];
    res += NLPSOL_NUM_OUT;
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

  Function Nlpsol::
  get_forward(casadi_int nfwd, const std::string& name,
              const std::vector<std::string>& inames,
              const std::vector<std::string>& onames,
              const Dict& opts) const {
    // Symbolic expression for the input
    vector<MX> arg = mx_in(), res = mx_out();

    // Initial guesses not used for derivative calculations
    for (NlpsolInput i : {NLPSOL_X0, NLPSOL_LAM_X0, NLPSOL_LAM_G0}) {
      arg[i] = MX::sym(arg[i].name(), Sparsity(arg[i].size()));
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

    // Calculates Hessian of the Lagrangian and Jacobian of the constraints
    Function HJ_fun = oracle_.factory("HJ", {"x", "p", "lam:f", "lam:g"},
      {"jac:g:x", "sym:hess:gamma:x:x"}, {{"gamma", {"f", "g"}}});

    // Hessian of the Lagrangian, Jacobian of the constraints
    vector<MX> HJ_res = HJ_fun({x, p, 1, lam_g});
    MX JG = HJ_res[0];
    MX HL = HJ_res[1];

    // Active bounds
    MX lam_x_pos = lam_x>0;
    MX lam_x_neg = lam_x<0;
    MX lam_g_pos = lam_g>0;
    MX lam_g_neg = lam_g<0;

    // Common
    MX alpha_x = if_else(lam_x_pos, ubx, 0) + if_else(lam_x_neg, lbx, 0);
    MX alpha_g = if_else(lam_g_pos, ubg, 0) + if_else(lam_g_neg, lbg, 0);
    MX a = alpha_x-x;

    // KKT matrix
    MX H_11 = mtimes(diag(a), HL) + diag(lam_x);
    MX H_12 = mtimes(diag(a), JG.T());
    MX H_21 = -mtimes(diag(lam_g), JG);
    MX H_22 = diag(alpha_g - g);
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

    // Propagate forward seeds
    MX fwd_alpha_x = if_else(lam_x_pos, fwd_ubx, 0) + if_else(lam_x_neg, fwd_lbx, 0);
    MX fwd_alpha_g = if_else(lam_g_pos, fwd_ubg, 0) + if_else(lam_g_neg, fwd_lbg, 0);
    MX v_x = -lam_x * fwd_alpha_x;
    MX v_lam_g = lam_g * fwd_alpha_g;
    MX v = MX::vertcat({v_x, v_lam_g});

    // Solve
    v = -MX::solve(H, v, "qr");

    // Extract sensitivities in x, lam_x and lam_g
    vector<MX> v_split = vertsplit(v, {0, nx_, nx_+ng_});
    MX fwd_x = v_split[0];
    MX fwd_lam_g = v_split[1];

    // Calculate sensitivities in f and g
    Function fwd_oracle = oracle_.forward(nfwd);
    vector<MX> vv = {x, p, f, g, fwd_x, fwd_p};
    vv = fwd_oracle(vv);
    MX fwd_f = vv[NL_F];
    MX fwd_g = vv[NL_G];

    // Calculate sensitivities in lam_x, lam_g
    Function rev_oracle = oracle_.reverse(1);
    // rev_reverse has the signature
    // (x, p, out_f, out_g, adj_f, adj_g) -> (adj_x, adj_p)
    // with adj_f=1, adj_g=lam_g, adj_x = -lam_x, adj_p = -lam_p
    Function fwd_rev_oracle = rev_oracle.forward(nfwd);
    // fwd_rev_oracle has the signature
    // (x, p, out_f, out_g, adj_f, adj_g, out_adj_x, out_adj_p,
    //  fwd_x, fwd_p, fwd_out_f, fwd_out_g, fwd_adj_f, fwd_adj_g)
    // -> (fwd_adj_x, fwd_adj_p)
    vv = {x, p, f, g, 1, lam_g, -lam_x, -lam_p,
          fwd_x, fwd_p, fwd_f, fwd_g, 0, fwd_lam_g};
    vv = fwd_rev_oracle(vv);
    MX fwd_lam_x = -vv[NL_X];
    MX fwd_lam_p = -vv[NL_P];

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

} // namespace casadi
