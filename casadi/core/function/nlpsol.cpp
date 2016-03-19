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

#include <iostream>
#include <iomanip>

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

  Nlpsol::Nlpsol(const std::string& name, Oracle* nlp)
    : FunctionInternal(name), nlp_(nlp) {

    // Set default options
    callback_step_ = 1;
    eval_errors_fatal_ = false;
    warn_initial_bounds_ = false;
    iteration_callback_ignore_errors_ = false;
    print_time_ = true;
  }

  Nlpsol::~Nlpsol() {
    if (nlp_) delete nlp_;
    clear_memory();
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
        "the different stages of initialization"}}
     }
  };

  void Nlpsol::init(const Dict& opts) {
    // Call the initialization method of the base class
    FunctionInternal::init(opts);

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
      }
    }

    // Replace MX oracle with SX oracle?
    if (expand && nlp_) {
      Oracle* nlp_new = nlp_->expand();
      delete nlp_;
      nlp_ = nlp_new;
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

  void Nlpsol::init_memory(void* mem) const {
    auto m = static_cast<NlpsolMemory*>(mem);

    // Create statistics
    for (const Function& f : all_functions_) {
      m->fstats[f.name()] = FStats();
    }

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

  void Nlpsol::set_temp(void* mem, const double** arg, double** res,
                        int* iw, double* w) const {
    auto m = static_cast<NlpsolMemory*>(mem);
    m->arg = arg;
    m->res = res;
    m->iw = iw;
    m->w = w;
  }

  // Convert a float to a string of an exact length.
  // First it tries fixed precision, then falls back to exponential notation.
  //
  // todo(jaeandersson,jgillis): needs either review or unit tests
  // because it throws exceptions if it fail.
  std::string formatFloat(double x, int totalWidth, int maxPrecision, int fallbackPrecision) {
    std::ostringstream out0;
    out0 << fixed << setw(totalWidth) << setprecision(maxPrecision) << x;
    std::string ret0 = out0.str();
    if (ret0.length() == totalWidth) {
      return ret0;
    } else if (ret0.length() > totalWidth) {
      std::ostringstream out1;
      out1 << setw(totalWidth) << setprecision(fallbackPrecision) << x;
      std::string ret1 = out1.str();
      if (ret1.length() != totalWidth)
        casadi_error(
          "ipopt timing formatting fallback is bugged, sorry about that."
          << "expected " << totalWidth <<  " digits, but got " << ret1.length()
          << ", string: \"" << ret1 << "\", number: " << x);
      return ret1;
    } else {
      casadi_error("ipopt timing formatting is bugged, sorry about that.");
    }
  }

  void print_stats_line(int maxNameLen, std::string label,
      double n_call, double t_proc, double t_wall) {
    // Skip when not called
    if (n_call == 0) return;

    std::stringstream s;

    s
      << setw(maxNameLen) << label << " "
      << formatFloat(t_proc, 9, 3, 3) << " [s]  "
      << formatFloat(t_wall, 9, 3, 3) << " [s]";
    if (n_call == -1) {
      // things like main loop don't have # evals
      s << endl;
    } else {
      s
        << " "
        << setw(5) << n_call;
      if (n_call < 2) {
        s << endl;
      } else {
        // only print averages if there is more than 1 eval
        s
          << " "
          << formatFloat(1000.0*t_proc/n_call, 10, 2, 3) << " [ms]  "
          << formatFloat(1000.0*t_wall/n_call, 10, 2, 3) << " [ms]"
          << endl;
      }
    }
    userOut() << s.str();
  }

  void Nlpsol::print_fstats(const NlpsolMemory* m) const {

    size_t maxNameLen=0;

    // Retrieve all nlp keys
    std::vector<std::string> keys;
    std::vector<std::string> keys_other;
    for (auto &&s : m->fstats) {
      maxNameLen = max(s.first.size(), maxNameLen);
      if (s.first.find("nlp")!=std::string::npos) {
        keys.push_back(s.first);
      } else if (s.first.find("mainloop")==std::string::npos) {
        keys_other.push_back(s.first);
      } else {
        continue;
      }
    }

    maxNameLen = max(std::string("all previous").size(), maxNameLen);
    maxNameLen = max(std::string("solver").size(), maxNameLen);

    // Print header
    std::stringstream s;
    std::string blankName(maxNameLen, ' ');
    s
      << blankName
      << "      proc           wall      num           mean             mean"
      << endl << blankName
      << "      time           time     evals       proc time        wall time";
    userOut() << s.str() << endl;

    // Sort the keys according to order
    std::vector<std::string> keys_order0;
    std::vector<std::string> keys_order1;
    std::vector<std::string> keys_order2;
    for (auto k : keys) {
      if (k.find("hess")!=std::string::npos) {
        keys_order2.push_back(k);
        continue;
      }
      if (k.find("grad")!=std::string::npos ||
          k.find("jac")!=std::string::npos) {
        keys_order1.push_back(k);
        continue;
      }
      keys_order0.push_back(k);
    }

    // Print all NLP stats
    for (auto keys : {&keys_order0, &keys_order1, &keys_order2}) {
        std::sort(keys->begin(), keys->end());
        for (auto k : *keys) {
          const FStats& fs = m->fstats.at(k);
          print_stats_line(maxNameLen, k, fs.n_call, fs.t_proc, fs.t_wall);
        }
    }

    // Sum the previously printed stats
    double t_wall_all_previous = 0;
    double t_proc_all_previous = 0;
    for (auto k : keys) {
      const FStats& fs = m->fstats.at(k);
      t_proc_all_previous += fs.t_proc;
      t_wall_all_previous += fs.t_wall;
    }
    print_stats_line(maxNameLen, "all previous", -1, t_proc_all_previous, t_wall_all_previous);

    // Sort and show the remainder of keys
    std::sort(keys_other.begin(), keys_other.end());
    for (std::string k : keys_other) {
      const FStats& fs = m->fstats.at(k);
      print_stats_line(maxNameLen, k, fs.n_call, fs.t_proc, fs.t_wall);
      t_proc_all_previous += fs.t_proc;
      t_wall_all_previous += fs.t_wall;
    }

    // Show the mainloop stats
    const FStats& fs_mainloop = m->fstats.at("mainloop");
    if (fs_mainloop.n_call>0) {
      print_stats_line(maxNameLen, "solver", -1,
        fs_mainloop.t_proc-t_proc_all_previous, fs_mainloop.t_wall-t_wall_all_previous);
      print_stats_line(maxNameLen, "mainloop", -1, fs_mainloop.t_proc, fs_mainloop.t_wall);
    }

  }

  Dict Nlpsol::get_stats(void *mem) const {
    auto m = static_cast<NlpsolMemory*>(mem);

    // Add timing statistics
    Dict stats;
    for (auto&& s : m->fstats) {
      stats["n_call_" +s.first] = s.second.n_call;
      stats["t_wall_" +s.first] = s.second.t_wall;
      stats["t_proc_" +s.first] = s.second.t_proc;
    }
    return stats;
  }

  int Nlpsol::calc_function(NlpsolMemory* m, const Function& fcn,
                            std::initializer_list<const double*> arg,
                            std::initializer_list<double*> res) const {
    // Respond to a possible Crl+C signals
    InterruptHandler::check();

    // Get statistics structure
    FStats& fstats = m->fstats.at(fcn.name());

    // Prepare stats, start timer
    fstats.tic();

    // Number of inputs and outputs
    int n_in = fcn.n_in(), n_out = fcn.n_out();

    // Input buffers
    fill_n(m->arg, n_in, nullptr);
    auto arg_it = arg.begin();
    for (int i=0; i<n_in; ++i) m->arg[i] = *arg_it++;
    casadi_assert(arg_it==arg.end());

    // Output buffers
    fill_n(m->res, n_out, nullptr);
    auto res_it = res.begin();
    for (int i=0; i<n_out; ++i) m->res[i] = *res_it++;
    casadi_assert(res_it==res.end());

    // Evaluate memory-less
    try {
      fcn(m->arg, m->res, m->iw, m->w, 0);
    } catch(exception& ex) {
      // Fatal error
      userOut<true, PL_WARN>()
        << name() << ":" << fcn.name() << " failed:" << ex.what() << endl;
      return 1;
    }

    // Make sure not NaN or Inf
    for (int i=0; i<n_out; ++i) {
      if (!m->res[i]) continue;
      if (!all_of(m->res[i], m->res[i]+fcn.nnz_out(i), [](double v) { return isfinite(v);})) {
        userOut<true, PL_WARN>()
          << name() << ":" << fcn.name() << " failed: NaN or Inf detected for output "
          << fcn.name_out(i) << endl;
        return -1;
      }
    }

    // Update stats
    fstats.toc();

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
    // Make sure that the function names are unique
    for (const Function& f : all_functions_) casadi_assert(fcn.name()!=f.name());

    // Add to the list
    all_functions_.push_back(fcn);
    alloc(fcn);
  }

  std::vector<bool> Nlpsol::nl_var(const std::string& s_in,
                                   const std::vector<std::string>& s_out) const {
    return nlp_->nl_var(s_in, s_out);
  }

} // namespace casadi
