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
#include "../misc/nlp_builder.hpp"

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
    return nlpsol(name, solver, Nlpsol::map2problem(nlp), opts);
  }

  Function nlpsol(const string& name, const string& solver,
                  const MXDict& nlp, const Dict& opts) {
    return nlpsol(name, solver, Nlpsol::map2problem(nlp), opts);
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

  double nlpsol_default_in(int ind) {
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

  Nlpsol::Nlpsol(const std::string& name, const Function& oracle)
    : OracleFunction(name, oracle) {

    // Set default options
    callback_step_ = 1;
    eval_errors_fatal_ = false;
    warn_initial_bounds_ = false;
    iteration_callback_ignore_errors_ = false;
    print_time_ = true;
  }

  Nlpsol::~Nlpsol() {
    clear_memory();
  }

  Sparsity Nlpsol::get_sparsity_in(int i) {
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

  Sparsity Nlpsol::get_sparsity_out(int i) {
    switch (static_cast<NlpsolOutput>(i)) {
    case NLPSOL_F:
      return Sparsity::scalar();
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
      {"print_time",
         {OT_BOOL,
          "print information about execution time"}},
      {"verbose_init",
       {OT_BOOL,
        "Print out timing information about "
        "the different stages of initialization"}},
      {"discrete",
       {OT_BOOLVECTOR,
        "Indicates which of the variables are discrete, i.e. integer-valued"}}
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
      } else if (op.first=="print_time") {
        print_time_ = op.second;
      } else if (op.first=="discrete") {
        discrete_ = op.second;
      }
    }

    // Replace MX oracle with SX oracle?
    if (expand) oracle_ = oracle_.expand();

    // Get dimensions
    nx_ = nnz_out(NLPSOL_X);
    np_ = nnz_in(NLPSOL_P);
    ng_ = nnz_out(NLPSOL_G);

    // Discrete marker
    mi_ = false;
    if (!discrete_.empty()) {
      casadi_assert_message(discrete_.size()==nx_, "\"discrete\" option has wrong length");
      if (std::find(discrete_.begin(), discrete_.end(), true)!=discrete_.end()) {
        casadi_assert_message(integer_support(),
                              "Discrete variables require a solver with integer support");
        mi_ = true;
      }
    }

    if (!fcallback_.is_null()) {
      // Consistency checks
      casadi_assert(!fcallback_.is_null());
      casadi_assert_message(fcallback_.n_out()==1 && fcallback_.numel_out()==1,
                            "Callback function must return a scalar.");
      casadi_assert_message(fcallback_.n_in()==n_out(),
                            "Callback input signature must match the NLP solver output signature");
      for (int i=0; i<n_out(); ++i) {
        casadi_assert_message(fcallback_.size_in(i)==size_out(i),
                              "Callback function input size mismatch. " <<
                              "For argument '" << nlpsol_out(i) << "', callback has shape "
                              << fcallback_.sparsity_in(i).dim() << " while NLP has " <<
                              sparsity_out(i).dim() << ".");
        // TODO(@jaeandersson): Wrap fcallback_ in a function with correct sparsity
        casadi_assert_message(fcallback_.sparsity_in(i)==sparsity_out(i),
                              "Callback function input size mismatch. " <<
                              "For argument " << nlpsol_out(i) << "', callback has shape "
                              << fcallback_.sparsity_in(i).dim() << " while NLP has " <<
                              sparsity_out(i).dim() << ".");
      }

      // Allocate temporary memory
      alloc(fcallback_);
    }
  }

  void Nlpsol::init_memory(void* mem) const {
    OracleFunction::init_memory(mem);
    auto m = static_cast<NlpsolMemory*>(mem);

    m->fstats["mainloop"] = FStats();
    m->fstats["callback_fun"] = FStats();
    m->fstats["callback_prep"] = FStats();
  }

  void Nlpsol::checkInputs(void* mem) const {
    auto m = static_cast<NlpsolMemory*>(mem);

    // Skip check?
    if (!inputs_check_) return;

    const double inf = std::numeric_limits<double>::infinity();

    // Number of equality constraints
    int n_eq = 0;

    // Detect ill-posed problems (simple bounds)
    for (int i=0; i<nx_; ++i) {
      double lbx = m->lbx ? m->lbx[i] : 0;
      double ubx = m->ubx ? m->ubx[i] : 0;
      double x0 = m->x0 ? m->x0[i] : 0;
      casadi_assert_message(!(lbx==inf || lbx>ubx || ubx==-inf),
                            "Ill-posed problem detected (x bounds)");
      if (warn_initial_bounds_ && (x0>ubx || x0<lbx)) {
        casadi_warning("Nlpsol: The initial guess does not satisfy LBX and UBX. "
                       "Option 'warn_initial_bounds' controls this warning.");
        break;
      }
      if (lbx==ubx) n_eq++;
    }

    // Detect ill-posed problems (nonlinear bounds)
    for (int i=0; i<ng_; ++i) {
      double lbg = m->lbg ? m->lbg[i] : 0;
      double ubg = m->ubg ? m->ubg[i] : 0;
      casadi_assert_message(!(lbg==inf || lbg>ubg || ubg==-inf),
                            "Ill-posed problem detected (g bounds)");
      if (lbg==ubg) n_eq++;
    }

    // Make sure enough degrees of freedom
    using casadi::to_string; // Workaround, MingGW bug, cf. CasADi issue #890
    casadi_assert_message(n_eq <= nx_, "NLP is overconstrained: There are " + to_string(n_eq)
                         + " equality constraints but only " + to_string(nx_) + " variables.");
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

  void Nlpsol::eval(void* mem, const double** arg, double** res, int* iw, double* w) const {
    // Reset the solver, prepare for solution
    setup(mem, arg, res, iw, w);

    // Solve the NLP
    solve(mem);

    // Show statistics
    if (print_time_)  print_fstats(static_cast<OracleMemory*>(mem));
  }

  void Nlpsol::set_work(void* mem, const double**& arg, double**& res,
                        int*& iw, double*& w) const {
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

} // namespace casadi
