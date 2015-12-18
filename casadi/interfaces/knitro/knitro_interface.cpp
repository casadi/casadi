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


#include "knitro_interface.hpp"
#include "casadi/core/std_vector_tools.hpp"
#include <ctime>
#include <cstdio>
#include <cstdlib>

using namespace std;
namespace casadi {

  extern "C"
  int CASADI_NLPSOL_KNITRO_EXPORT
  casadi_register_nlpsol_knitro(Nlpsol::Plugin* plugin) {
    plugin->creator = KnitroInterface::creator;
    plugin->name = "knitro";
    plugin->doc = KnitroInterface::meta_doc.c_str();
    plugin->version = 23;
    return 0;
  }

  extern "C"
  void CASADI_NLPSOL_KNITRO_EXPORT casadi_load_nlpsol_knitro() {
    Nlpsol::registerPlugin(casadi_register_nlpsol_knitro);
  }

  KnitroInterface::KnitroInterface(const std::string& name, const XProblem& nlp)
    : Nlpsol(name, nlp) {

    // Monitors
    addOption("monitor",      OT_STRINGVECTOR, GenericType(),  "",
              "eval_f|eval_g|eval_jac_g|eval_grad_f|eval_h", true);

    // Not yet ready
    //addOption("algorithm",                OT_STRING, GenericType(),
    // "Which algorithm to use. See KNITRO documentation.", "auto|direct|cg|active");
    //addOption("bar_directinterval",       OT_INTEGER, GenericType(),
    //  "When using the Interior/Direct algorithm, this parameter controls the maximum number of "
    //  "consecutive CG steps before trying to force the algorithm to take a direct step again. "
    //  "See KNITRO documentation.");
    //addOption("bar_feasible",             OT_STRING, GenericType(),
    //  "Whether feasibility is given special emphasis. See KNITRO documentation.",
    //  "no|stay|get|get_stay");
    //addOption("bar_feasmodetol",          OT_REAL, GenericType(),
    //  "Specifies the tolerance for entering the stay feasible mode See KNITRO documentation.");
    //addOption("bar_initmu",               OT_INTEGER, GenericType(),
    //  "Initial value for the barrier parameter. See KNITRO documentation.");
    //addOption("bar_initpt",               OT_STRING, GenericType(),
    //  "Whether to use the initial point strategy with barrier algorithms. "
    //  "See KNITRO documentation.", "auto|yes|no");
    //addOption("bar_maxbacktrack",         OT_INTEGER, GenericType(),
    //  "Maximum allowable number of backtracks during the linesearch of the Interior Direct "
    //  "algorithm before reverting to a CG step. See KNITRO documentation.");
    //addOption("bar_maxrefactor",          OT_INTEGER, GenericType(),
    //  "Maximum number of refactorizations of the KKT system per iteration of the Interior "
    //  "Direct algorithm before reverting to a CG step. See KNITRO documentation.");

    //addOption("Alg", OT_INTEGER,0, "Algorithm");
    addOption("BarRule", OT_INTEGER, 0, "Barrier Rule");
    addOption("NewPoint", OT_BOOLEAN, 0, "Select new-point feature");
    addOption("GradOpt", OT_INTEGER, 1, "Gradient calculation method");
    addOption("HessOpt", OT_INTEGER, 1, "Hessian calculation method");
    addOption("Feasible", OT_BOOLEAN, 1, "Allow infeasible iterations");
    addOption("HonorBnds", OT_BOOLEAN, 0, "Enforce bounds");
    addOption("LpSolver", OT_BOOLEAN, 0, "Use LpSolver");
    addOption("Multistart", OT_BOOLEAN, 0, "Use multistart");
    //addOption("MsMaxSolves", OT_INTEGER, 1, "Maximum multistart points");
    addOption("MaxCgIt", OT_INTEGER, 0, "Maximum conjugate gradient iterations");
    //addOption("MaxCrossTt", OT_INTEGER, 0, "Maximum crossover iterations");
    addOption("MaxIt", OT_INTEGER, 10000, "Iteration limit");
    //addOption("MaxTimeCPU", OT_REAL, 1e8, "CPU Time limit");
    //addOption("MaxTimeReal", OT_REAL, 1e8, "Time limit");
    addOption("LmSize", OT_INTEGER, 10, "Memory pairsize limit");
    addOption("Scale", OT_BOOLEAN, 1, "Perform scaling");
    addOption("ShiftInit", OT_BOOLEAN, 1, "Interior-point shifting initial point");
    addOption("Soc", OT_INTEGER, 1, "Second order correction");
    addOption("InitPt", OT_BOOLEAN, 0, "Use initial point strategy");
    addOption("Delta", OT_REAL, 1.0, "Initial region scaling factor");
    addOption("FeasModeTol", OT_REAL, 1e-4, "Feasible mode tolerance");
    addOption("FeasTol", OT_REAL, 1e-6, "Feasible tolerance");
    addOption("FeasTolAbs", OT_REAL, 0, "Absolute feasible tolerance");
    addOption("OptTol", OT_REAL, 1e-6, "Relative optimality tolerance");
    addOption("OptTolAbs", OT_REAL, 0, "Absolute optimality tolerance");
    addOption("Pivot", OT_REAL, 1e-8, "Initial pivot threshold");
    addOption("XTol", OT_REAL, 1e-15, "Relative solution change tolerance");
    addOption("Mu", OT_REAL, 0.1, "Initial barrier parameter");
    addOption("ObjRange", OT_REAL, 1e20, "Maximum objective value");
    addOption("OutLev", OT_INTEGER, 2, "Log output level");
    addOption("Debug", OT_INTEGER, 0, "Debug level");



    kc_handle_ = 0;

    addOption("contype", OT_INTEGERVECTOR);
  }


  KnitroInterface::~KnitroInterface() {
    // Free KNITRO memory
    if (kc_handle_) {
      /*    KTR_free(&kc_handle_);
            kc_handle_ = 0;*/
    }
  }

  void KnitroInterface::init() {
    // Call the init method of the base class
    Nlpsol::init();

    //if (hasSetOption("Alg")) int_param_["alg"] = option("Alg");
    if (hasSetOption("BarRule")) int_param_["barrule"] = option("BarRule");
    if (hasSetOption("NewPoint")) int_param_["newpoint"] = option("NewPoint");
    if (hasSetOption("GradOpt")) int_param_["gradopt"] = option("GradOpt");
    if (hasSetOption("HessOpt")) int_param_["hessopt"] = option("HessOpt");
    if (hasSetOption("Feasible")) int_param_["feasible"] = option("Feasible");
    if (hasSetOption("HonorBnds")) int_param_["honorbnds"] = option("HonorBnds");
    if (hasSetOption("LpSolver")) int_param_["lpsolver"] = option("LpSolver");
    if (hasSetOption("Multistart")) int_param_["multistart"] = option("Multistart");
    //if (hasSetOption("MsMaxSolves")) int_param_["msmaxsolves"] = option("MsMaxSolves");
    if (hasSetOption("MaxCgIt")) int_param_["maxcgit"] = option("MaxCgIt");
    //if (hasSetOption("MaxCrossTt")) int_param_["maxcrosstt"] = option("MaxCrossTt");
    if (hasSetOption("MaxIt")) int_param_["maxit"] = option("MaxIt");
    //if (hasSetOption("MaxTimeCPU")) double_param_["maxtimecpu"] = option("MaxTimeCPU");
    //if (hasSetOption("MaxTimeReal")) double_param_["maxtimereal"] = option("MaxTimeReal");
    if (hasSetOption("LmSize")) int_param_["lmsize"] = option("LmSize");
    if (hasSetOption("Scale")) int_param_["scale"] = option("Scale");
    if (hasSetOption("ShiftInit")) int_param_["shiftinit"] = option("ShiftInit");
    if (hasSetOption("Soc")) int_param_["soc"] = option("Soc");
    if (hasSetOption("InitPt")) int_param_["initpt"] = option("InitPt");
    if (hasSetOption("Delta")) double_param_["delta"] = option("Delta");
    if (hasSetOption("FeasModeTol")) double_param_["feasmodetol"] = option("FeasModeTol");
    if (hasSetOption("FeasTol")) double_param_["feastol"] = option("FeasTol");
    if (hasSetOption("FeasTolAbs")) double_param_["feastolabs"] = option("FeasTolAbs");
    if (hasSetOption("OptTol")) double_param_["opttol"] = option("OptTol");
    if (hasSetOption("OptTolAbs")) double_param_["opttolabs"] = option("OptTolAbs");
    if (hasSetOption("Pivot")) double_param_["pivot"] = option("Pivot");
    if (hasSetOption("XTol")) double_param_["xtol"] = option("XTol");
    if (hasSetOption("Mu")) double_param_["mu"] = option("Mu");
    if (hasSetOption("ObjRange")) double_param_["objrange"] = option("ObjRange");
    if (hasSetOption("OutLev")) int_param_["outlev"] = option("OutLev");
    if (hasSetOption("Debug")) int_param_["debug"] = option("Debug");

    // Get/generate required functions
    gradF();
    jacG();
    if (true) { // NOTE: should be only if HessOpt
      hessLag();
    }

    // Commented out since I have not found out how to change the bounds
    // Allocate KNITRO memory block
    /*  casadi_assert(kc_handle_==0);
        kc_handle_ = KTR_new();*/

    // Allocate memory
    alloc_w(nx_, true); // wx_
    alloc_w(nx_, true); // wlbx_
    alloc_w(nx_, true); // wubx_
    alloc_w(ng_, true); // wlbg_
    alloc_w(ng_, true); // wubg_
  }

  void KnitroInterface::reset(void* mem, const double**& arg, double**& res, int*& iw, double*& w) {
    // Reset the base classes
    Nlpsol::reset(mem, arg, res, iw, w);

    // Copy inputs to temporary arrays
    wx_ = w; w += nx_;
    wlbx_ = w; w += nx_;
    wubx_ = w; w += nx_;
    wlbg_ = w; w += ng_;
    wubg_ = w; w += ng_;
  }

  void KnitroInterface::solve(void* mem) {
    // Allocate KNITRO memory block (move back to init!)
    casadi_assert(kc_handle_==0);
    kc_handle_ = KTR_new();
    casadi_assert(kc_handle_!=0);
    int status;

    // Jacobian sparsity
    vector<int> Jcol, Jrow;
    if (!jacG_.isNull()) {
      Jcol = jacG_.sparsity_out(0).get_col();
      int sz = jacG_.nnz_out(0);
      const int* row = jacG_.sparsity_out(0).row();
      Jrow = vector<int>(row, row+sz);
    }

    // Hessian sparsity
    int nnzH = hessLag_.isNull() ? 0 : hessLag_.output().nnz_lower();
    vector<int> Hcol(nnzH), Hrow(nnzH);
    if (nnzH>0) {
      const int* colind = hessLag_.sparsity_out(0).colind();
      int ncol = hessLag_.size2_out(0);
      const int* row = hessLag_.sparsity_out(0).row();
      int nz=0;
      for (int cc=0; cc<ncol; ++cc) {
        for (int el=colind[cc]; el<colind[cc+1] && row[el]<=cc; ++el) {
          Hcol[nz] = cc;
          Hrow[nz] = row[el];
          nz++;
        }
      }
      casadi_assert(nz==nnzH);

      status = KTR_set_int_param_by_name(kc_handle_, "hessopt", KTR_HESSOPT_EXACT);
      casadi_assert_message(status==0, "KTR_set_int_param failed");
    } else {
      status = KTR_set_int_param_by_name(kc_handle_, "hessopt", KTR_HESSOPT_LBFGS);
      casadi_assert_message(status==0, "KTR_set_int_param failed");
    }

    // Set user set options
    for (auto it=double_param_.begin(); it!=double_param_.end(); ++it) {
      status = KTR_set_double_param_by_name(kc_handle_, it->first.c_str(), it->second);
      if (status!=0) {
        throw CasadiException("KnitroInterface::evaluate: cannot set " + it->first);
      }
    }

    for (std::map<std::string, int>::iterator it=int_param_.begin(); it!=int_param_.end(); ++it) {
      status = KTR_set_int_param_by_name(kc_handle_, it->first.c_str(), it->second);
      if (status!=0) {
        throw CasadiException("KnitroInterface::evaluate: cannot set " + it->first);
      }
    }

    for (std::map<std::string, std::string>::iterator it=string_param_.begin();
        it!=string_param_.end(); ++it) {
      status = KTR_set_char_param_by_name(kc_handle_, it->first.c_str(), it->second.c_str());
      if (status!=0) {
        throw CasadiException("KnitroInterface::evaluate: cannot set " + it->first);
      }
    }

    // Type of constraints
    vector<int> cType(ng_, KTR_CONTYPE_GENERAL);
    if (hasSetOption("contype")) {
      vector<int> contype = option("contype");
      casadi_assert(contype.size()==cType.size());
      copy(contype.begin(), contype.end(), cType.begin());
    }

    // "Correct" upper and lower bounds
    casadi_copy(x0_, nx_, wx_);
    casadi_copy(lbx_, nx_, wlbx_);
    casadi_copy(ubx_, nx_, wubx_);
    casadi_copy(lbg_, ng_, wlbg_);
    casadi_copy(ubg_, ng_, wubg_);
    for (int i=0; i<nx_; ++i) if (isinf(wlbx_[i])) wlbx_[i] = -KTR_INFBOUND;
    for (int i=0; i<nx_; ++i) if (isinf(wubx_[i])) wubx_[i] =  KTR_INFBOUND;
    for (int i=0; i<ng_; ++i) if (isinf(wlbg_[i])) wlbg_[i] = -KTR_INFBOUND;
    for (int i=0; i<ng_; ++i) if (isinf(wubg_[i])) wubg_[i] =  KTR_INFBOUND;

    // Initialize KNITRO
    status = KTR_init_problem(kc_handle_, nx_, KTR_OBJGOAL_MINIMIZE,
                              KTR_OBJTYPE_GENERAL, wlbx_, wubx_, ng_, getPtr(cType),
                              wlbg_, wubg_, Jcol.size(), getPtr(Jcol), getPtr(Jrow),
                              nnzH, getPtr(Hrow), getPtr(Hcol), wx_, 0); // initial lambda
    casadi_assert_message(status==0, "KTR_init_problem failed");

    // Register callback functions
    status = KTR_set_func_callback(kc_handle_, &callback);
    casadi_assert_message(status==0, "KTR_set_func_callback failed");

    status = KTR_set_grad_callback(kc_handle_, &callback);
    casadi_assert_message(status==0, "KTR_set_grad_callbackfailed");

    if (nnzH>0) {
      status = KTR_set_hess_callback(kc_handle_, &callback);
      casadi_assert_message(status==0, "KTR_set_hess_callbackfailed");
    }

    // Lagrange multipliers
    vector<double> lambda(nx_+ng_);

    // Solve NLP
    double f;
    status = KTR_solve(kc_handle_, wx_, getPtr(lambda), 0, &f,
                       0, 0, 0, 0, 0, static_cast<void*>(this));

    // Output primal solution
    casadi_copy(wx_, nx_, x_);

    // Output dual solution
    casadi_copy(getPtr(lambda), ng_, lam_g_);
    casadi_copy(getPtr(lambda)+ng_, nx_, lam_x_);

    // Output optimal cost
    if (f_) *f_ = f;

    stats_["return_status"] = status;

    // Copy constraints
    casadi_copy(nlp_.output(NL_G).ptr(), ng_, g_);

    // Free memory (move to destructor!)
    KTR_free(&kc_handle_);
    kc_handle_ = 0;
  }


  int KnitroInterface::callback(const int evalRequestCode, const int n, const int m, const int nnzJ,
                               const int nnzH, const double* const x, const double* const lambda,
                               double* const obj, double* const c, double* const objGrad,
                               double* const jac, double* const hessian, double* const hessVector,
                               void *userParams) {
    try {
      // Get a pointer to the calling object
      KnitroInterface* this_ = static_cast<KnitroInterface*>(userParams);

      // Direct to the correct function
      switch (evalRequestCode) {
      case KTR_RC_EVALFC: this_->evalfc(x, *obj, c); break;
      case KTR_RC_EVALGA: this_->evalga(x, objGrad, jac); break;
      case KTR_RC_EVALH:  this_->evalh(x, lambda, hessian); break;
      default: casadi_assert_message(0, "KnitroInterface::callback: unknown method");
      }

      return 0;
    } catch(exception& ex) {
      userOut<true, PL_WARN>() << "KnitroInterface::callback caugth exception: "
                               << ex.what() << endl;
      return -1;
    }
  }

  void KnitroInterface::evalfc(const double* x, double& obj, double *c) {
    // Pass the argument to the function
    nlp_.setInputNZ(x, NL_X);
    if (p_) {
      nlp_.setInputNZ(p_, NL_P);
    } else {
      nlp_.setInput(0., NL_P);
    }

    // Evaluate the function
    nlp_.evaluate();

    // Get the result
    nlp_.output(NL_F).get(obj);
    nlp_.output(NL_G).get(c);

    // Printing
    if (monitored("eval_f")) {
      userOut() << "x = " << nlp_.input(NL_X) << endl;
      userOut() << "f = " << nlp_.output(NL_F) << endl;
    }
    if (monitored("eval_g")) {
      userOut() << "x = " << nlp_.input(NL_X) << endl;
      userOut() << "g = " << nlp_.output(NL_G) << endl;
    }
  }

  void KnitroInterface::evalga(const double* x, double* objGrad, double* jac) {
    // Pass the argument to the function
    gradF_.setInputNZ(x, NL_X);
    if (p_) {
      gradF_.setInputNZ(p_, NL_P);
    } else {
      gradF_.setInput(0., NL_P);
    }

    // Evaluate the function using adjoint mode AD
    gradF_.evaluate();

    // Get the result
    gradF_.output().get(objGrad);

    // Printing
    if (monitored("eval_grad_f")) {
      userOut() << "x = " << gradF_.input(NL_X) << endl;
      userOut() << "grad_f = " << gradF_.output() << endl;
    }

    if (!jacG_.isNull()) {
      // Pass the argument to the Jacobian function
      jacG_.setInputNZ(x, NL_X);
      if (p_) {
        jacG_.setInputNZ(p_, NL_P);
      } else {
        jacG_.setInput(0., NL_P);
      }

      // Evaluate the Jacobian function
      jacG_.evaluate();

      // Get the result
      jacG_.output().getNZ(jac);

      // Printing
      if (monitored("eval_jac_g")) {
        userOut() << "x = " << jacG_.input(NL_X) << endl;
        userOut() << "jac_g = " << jacG_.output() << endl;
      }
    }
  }

  void KnitroInterface::evalh(const double* x, const double* lambda, double* hessian) {
    // Pass the argument to the function
    hessLag_.setInputNZ(x, NL_X);
    if (p_) {
      hessLag_.setInputNZ(p_, NL_P);
    } else {
      hessLag_.setInput(0., NL_P);
    }
    hessLag_.setInput(1.0, NL_NUM_IN+NL_F);
    hessLag_.setInputNZ(lambda, NL_NUM_IN+NL_G);

    // Evaluate
    hessLag_.evaluate();

    // Get results
    hessLag_.output().getSym(hessian);

    // Printing
    if (monitored("eval_h")) {
      userOut() << "eval_h" << endl;
      userOut() << "x = " << hessLag_.input(0) << endl;
      userOut() << "lambda = " << hessLag_.input(1) << endl;
      userOut() << "scale = " << hessLag_.input(2) << endl;
      userOut() << "H = " << hessLag_ << endl;
    }
  }

} // namespace casadi
