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
#include "casadi/core/timing.hpp"
#include <chrono>

using namespace std;
namespace casadi {

  bool has_nlpsol(const string& name) {
    return Nlpsol::hasPlugin(name);
  }

  void load_nlpsol(const string& name) {
    Nlpsol::loadPlugin(name);
  }

  string doc_nlpsol(const string& name) {
    return Nlpsol::getPlugin(name).doc;
  }

  Function nlpsol(const string& name, const string& solver,
                  const SXDict& nlp, const Dict& opts) {
    return nlpsol(name, solver, Nlpsol::map2problem(nlp), opts);
  }

  Function nlpsol(const string& name, const string& solver,
                  const MXDict& nlp, const Dict& opts) {
    return nlpsol(name, solver, Nlpsol::map2problem(nlp), opts);
  }

  Function nlpsol(const std::string& name, const std::string& solver,
                  const std::string& fname, const Dict& opts) {
    return nlpsol(name, solver, Oracle::construct(fname, "nlp"), opts);
  }

  Function nlpsol(const string& name, const string& solver,
                  const Compiler& compiler, const Dict& opts) {
    return nlpsol(name, solver, Oracle::construct(compiler, "nlp"), opts);
  }

  Function nlpsol(const string& name, const string& solver,
                  const Function& nlp, const Dict& opts) {
    if (nlp.is_a("sxfunction")) {
      return nlpsol(name, solver, Nlpsol::fun2problem<SX>(nlp), opts);
    } else {
      return nlpsol(name, solver, Nlpsol::fun2problem<MX>(nlp), opts);
    }
  }

  Function nlpsol(const string& name, const string& solver,
                  Oracle* nlp, const Dict& opts) {
    Function ret;
    ret.assignNode(Nlpsol::instantiatePlugin(name, solver, nlp));
    ret->construct(opts);
    return ret;
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

  string nlpsol_in(int ind) {
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

  string nlpsol_out(int ind) {
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

  int nlpsol_n_in() {
    return NLPSOL_NUM_IN;
  }

  int nlpsol_n_out() {
    return NLPSOL_NUM_OUT;
  }

  Nlpsol::Nlpsol(const std::string& name, Oracle* nlp)
    : FunctionInternal(name), nlp_(nlp) {

    // Set default options
    callback_step_ = 1;
    eval_errors_fatal_ = false;
    warn_initial_bounds_ = false;
    iteration_callback_ignore_errors_ = false;
  }

  Nlpsol::~Nlpsol() {
    if (nlp_) delete nlp_;
  }

  Sparsity Nlpsol::get_sparsity_in(int ind) const {
    switch (static_cast<NlpsolInput>(ind)) {
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
      return nlp_->sparsity_in(NL_P);
    case NLPSOL_NUM_IN: break;
    }
    return Sparsity();
  }

  Sparsity Nlpsol::get_sparsity_out(int ind) const {
    switch (static_cast<NlpsolOutput>(ind)) {
    case NLPSOL_F:
      return Sparsity::scalar();
    case NLPSOL_X:
    case NLPSOL_LAM_X:
      return nlp_->sparsity_in(NL_X);
    case NLPSOL_LAM_G:
    case NLPSOL_G:
      return nlp_->sparsity_out(NL_G);
    case NLPSOL_LAM_P:
      return get_sparsity_in(NLPSOL_P);
    case NLPSOL_NUM_OUT: break;
    }
    return Sparsity();
  }

  Options Nlpsol::options_
  = {{&FunctionInternal::options_},
     {{"expand",
       {OT_BOOL,
        "Expand the NLP function in terms of scalar operations, i.e. MX->SX"}},
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
        "the different stages of initialization"}}
     }
  };

  void Nlpsol::init(const Dict& opts) {
    // Call the initialization method of the base class
    FunctionInternal::init(opts);

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
      }
    }

    // Get dimensions
    nx_ = nnz_out(NLPSOL_X);
    np_ = nnz_in(NLPSOL_P);
    ng_ = nnz_out(NLPSOL_G);

    if (!fcallback_.is_null()) {
      // Consistency checks
      casadi_assert(!fcallback_.is_null());
      casadi_assert_message(fcallback_.n_out()==1 && fcallback_.numel_out()==1,
                            "Callback function must return a scalar");
      casadi_assert_message(fcallback_.n_in()==n_out(),
                            "Callback input signature must match the NLP solver output signature");
      for (int i=0; i<n_out(); ++i) {
        casadi_assert_message(fcallback_.size_in(i)==size_out(i),
                              "Callback function input size mismatch");
        // TODO(@jaeandersson): Wrap fcallback_ in a function with correct sparsity
        casadi_assert_message(fcallback_.sparsity_in(i)==sparsity_out(i),
                              "Not implemented");
      }

      // Allocate temporary memory
      alloc(fcallback_);
    }
  }

  void Nlpsol::init_memory(Memory& mem) const {
    NlpsolMemory& m = dynamic_cast<NlpsolMemory&>(mem);

    // Create statistics
    for (const Function& f : all_functions_) {
      m.fstats[f.name()] = NlpsolMemory::FStats();
    }
  }

  void Nlpsol::checkInputs(Memory& mem) const {
    NlpsolMemory& m = dynamic_cast<NlpsolMemory&>(mem);

    // Skip check?
    if (!inputs_check_) return;

    const double inf = std::numeric_limits<double>::infinity();

    // Detect ill-posed problems (simple bounds)
    for (int i=0; i<nx_; ++i) {
      double lbx = m.lbx ? m.lbx[i] : 0;
      double ubx = m.ubx ? m.ubx[i] : 0;
      double x0 = m.x0 ? m.x0[i] : 0;
      casadi_assert_message(!(lbx==inf || lbx>ubx || ubx==-inf),
                            "Ill-posed problem detected (x bounds)");
      if (warn_initial_bounds_ && (x0>ubx || x0<lbx)) {
        casadi_warning("Nlpsol: The initial guess does not satisfy LBX and UBX. "
                       "Option 'warn_initial_bounds' controls this warning.");
        break;
      }
    }

    // Detect ill-posed problems (nonlinear bounds)
    for (int i=0; i<ng_; ++i) {
      double lbg = m.lbg ? m.lbg[i] : 0;
      double ubg = m.ubg ? m.ubg[i] : 0;
      casadi_assert_message(!(lbg==inf || lbg>ubg || ubg==-inf),
                            "Ill-posed problem detected (g bounds)");
    }
  }

  std::map<std::string, Nlpsol::Plugin> Nlpsol::solvers_;

  const std::string Nlpsol::infix_ = "nlpsol";

  DM Nlpsol::getReducedHessian() {
    casadi_error("Nlpsol::getReducedHessian not defined for class "
                 << typeid(*this).name());
    return DM();
  }

  void Nlpsol::setOptionsFromFile(const std::string & file) {
    casadi_error("Nlpsol::setOptionsFromFile not defined for class "
                 << typeid(*this).name());
  }

  double Nlpsol::default_in(int ind) const {
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

  void Nlpsol::eval(Memory& mem, const double** arg, double** res, int* iw, double* w) const {
    // Reset the solver, prepare for solution
    setup(mem, arg, res, iw, w);

    // Solve the NLP
    solve(mem);
  }

  void Nlpsol::set_work(Memory& mem, const double**& arg, double**& res,
                        int*& iw, double*& w) const {
    NlpsolMemory& m = dynamic_cast<NlpsolMemory&>(mem);

    // Get input pointers
    m.x0 = arg[NLPSOL_X0];
    m.p = arg[NLPSOL_P];
    m.lbx = arg[NLPSOL_LBX];
    m.ubx = arg[NLPSOL_UBX];
    m.lbg = arg[NLPSOL_LBG];
    m.ubg = arg[NLPSOL_UBG];
    m.lam_x0 = arg[NLPSOL_LAM_X0];
    m.lam_g0 = arg[NLPSOL_LAM_G0];
    arg += NLPSOL_NUM_IN;

    // Get output pointers
    m.x = res[NLPSOL_X];
    m.f = res[NLPSOL_F];
    m.g = res[NLPSOL_G];
    m.lam_x = res[NLPSOL_LAM_X];
    m.lam_g = res[NLPSOL_LAM_G];
    m.lam_p = res[NLPSOL_LAM_P];
    res += NLPSOL_NUM_OUT;
  }

  void Nlpsol::set_temp(Memory& mem, const double** arg, double** res,
                        int* iw, double* w) const {
    NlpsolMemory& m = dynamic_cast<NlpsolMemory&>(mem);
    m.arg = arg;
    m.res = res;
    m.iw = iw;
    m.w = w;
  }

  int Nlpsol::calc_function(NlpsolMemory& m, const Function& fcn,
                            std::initializer_list<const double*> arg,
                            std::initializer_list<double*> res) const {
    // Respond to a possible Crl+C signals
    InterruptHandler::check();

    // Get statistics structure
    NlpsolMemory::FStats& fstats = m.fstats.at(fcn.name());

    // Number of inputs and outputs
    int n_in = fcn.n_in(), n_out = fcn.n_out();

    // Input buffers
    fill_n(m.arg, n_in, nullptr);
    auto arg_it = arg.begin();
    for (int i=0; i<n_in; ++i) m.arg[i] = *arg_it++;
    casadi_assert(arg_it==arg.end());

    // Output buffers
    fill_n(m.res, n_out, nullptr);
    auto res_it = res.begin();
    for (int i=0; i<n_out; ++i) m.res[i] = *res_it++;
    casadi_assert(res_it==res.end());

    // Prepare stats, start timer
    fstats.n_calc += 1;
    auto t_start = chrono::system_clock::now();

    // Evaluate memory-less
    try {
      fcn(m.arg, m.res, m.iw, m.w, 0);
    } catch(exception& ex) {
      // Fatal error
      userOut<true, PL_WARN>()
        << name() << ":" << fcn.name() << " failed:" << ex.what() << endl;
      return 1;
    }

    // Stop timer
    auto t_stop = chrono::system_clock::now();

    // Make sure not NaN or Inf
    for (int i=0; i<n_out; ++i) {
      if (!m.res[i]) continue;
      if (!all_of(m.res[i], m.res[i]+fcn.nnz_out(i), [](double v) { return isfinite(v);})) {
        userOut<true, PL_WARN>()
          << name() << ":" << fcn.name() << " failed: NaN or Inf detected for output "
          << fcn.name_out(i) << endl;
        return -1;
      }
    }

    // Update stats
    fstats.n_calc += 1;
    fstats.t_calc += chrono::duration<double>(t_stop - t_start).count();

    // Success
    return 0;
  }

  int Nlpsol::
  calc_fg(NlpsolMemory& m, const Function& fcn,
          const double* x, const double* p,
          double* f, double* g) const {
    fill_n(m.arg, fcn.n_in(), nullptr);
    m.arg[0] = x;
    m.arg[1] = p;
    fill_n(m.res, fcn.n_out(), nullptr);
    m.res[0] = f;
    m.res[1] = g;
    fcn(m.arg, m.res, m.iw, m.w, 0);

    // Success
    return 0;
  }

  int Nlpsol::
  calc_gf_jg(NlpsolMemory& m, const Function& fcn,
             const double* x, const double* p,
             double* gf, double* jg) const {
    fill_n(m.arg, fcn.n_in(), nullptr);
    m.arg[0] = x;
    m.arg[1] = p;
    fill_n(m.res, fcn.n_out(), nullptr);
    m.res[0] = gf;
    m.res[1] = jg;
    fcn(m.arg, m.res, m.iw, m.w, 0);

    // Success
    return 0;
  }

  int Nlpsol::calc_jac_g(NlpsolMemory& m, const Function& fcn,
                         const double* x, const double* p,
                         double* g, double* jac_g) const {
    // Respond to a possible Crl+C signals
    InterruptHandler::check();
    casadi_assert(jac_g!=0);

    // Evaluate User function
    fill_n(m.arg, fcn.n_in(), nullptr);
    m.arg[0] = x;
    m.arg[1] = p;
    fill_n(m.res, fcn.n_out(), nullptr);
    m.res[0] = g;
    m.res[1] = jac_g;
    fcn(m.arg, m.res, m.iw, m.w, 0);

    // Success
    return 0;
  }

  int Nlpsol::calc_jac_f(NlpsolMemory& m, const Function& fcn,
                         const double* x, const double* p,
                         double* f, double* jac_f) const {
    // Respond to a possible Crl+C signals
    InterruptHandler::check();
    casadi_assert(jac_f!=0);

    // Evaluate User function
    fill_n(m.arg, fcn.n_in(), nullptr);
    m.arg[0] = x;
    m.arg[1] = p;
    fill_n(m.res, fcn.n_out(), nullptr);
    m.res[0] = f;
    m.res[1] = jac_f;
    fcn(m.arg, m.res, m.iw, m.w, 0);

    // Success
    return 0;
  }

  int Nlpsol::calc_hess_l(NlpsolMemory& m, const Function& fcn,
                          const double* x, const double* p,
                          const double* sigma, const double* lambda,
                          double* hl) const {
    // Respond to a possible Crl+C signals
    InterruptHandler::check();

    // Evaluate User function
    fill_n(m.arg, fcn.n_in(), nullptr);
    m.arg[HL_X] = x;
    m.arg[HL_P] = p;
    m.arg[HL_LAM_F] = sigma;
    m.arg[HL_LAM_G] = lambda;
    fill_n(m.res, fcn.n_out(), nullptr);
    m.res[HL_HL] = hl;
    auto t_start = chrono::system_clock::now(); // start timer
    try {
      fcn(m.arg, m.res, m.iw, m.w, 0);
    } catch(exception& ex) {
      // Fatal error
      userOut<true, PL_WARN>() << name() << ":calc_hess_l failed:" << ex.what() << endl;
      return 1;
    }
    auto t_stop = chrono::system_clock::now(); // stop timer

    // Make sure not NaN or Inf
    if (!all_of(hl, hl+hesslag_sp_.nnz(), [](double v) { return isfinite(v);})) {
      userOut<true, PL_WARN>() << name() << ":calc_hess_l failed: NaN or Inf detected" << endl;
      return -1;
    }

    // Update stats
    m.n_calc_hess_l += 1;
    m.t_calc_hess_l += chrono::duration<double>(t_stop - t_start).count();

    // Success
    return 0;
  }

  void Nlpsol::generate_dependencies(const std::string& fname, const Dict& opts) {
    CodeGenerator gen(opts);
    gen.add(nlp_->all_io("nlp"));
    for (const Function& f : all_functions_) gen.add(f);
    gen.generate(fname);
  }

  Function Nlpsol::create_function(const std::string& fname,
                                   const std::vector<std::string>& s_in,
                                   const std::vector<std::string>& s_out,
                                   const std::vector<LinComb>& lincomb,
                                   const Dict& opts, bool reg) {
    // Generate the function
    casadi_assert(nlp_!=0);
    Function ret = nlp_->create(fname, s_in, s_out, lincomb, opts);
    if (reg) register_function(ret);
    return ret;
  }

  void Nlpsol::register_function(const Function& fcn) {
    all_functions_.push_back(fcn);
    alloc(fcn);
  }

  std::vector<bool> Nlpsol::nl_var(const std::string& s_in,
                                   const std::vector<std::string>& s_out) const {
    return nlp_->nl_var(s_in, s_out);
  }

} // namespace casadi
