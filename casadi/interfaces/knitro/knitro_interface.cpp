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
#include "casadi/core/matrix/matrix_tools.hpp"
#include <ctime>
#include <cstdio>
#include <cstdlib>

using namespace std;
namespace casadi {

  extern "C"
  int CASADI_NLP_KNITRO_EXPORT
  casadi_register_nlp_knitro(NlpSolverInternal::Plugin* plugin) {
    plugin->creator = KnitroInterface::creator;
    plugin->name = "knitro";
    plugin->doc = KnitroInterface::meta_doc.c_str();
    plugin->version = 21;
    return 0;
  }

  extern "C"
  void CASADI_NLP_KNITRO_EXPORT casadi_load_nlp_knitro() {
    NlpSolverInternal::registerPlugin(casadi_register_nlp_knitro);
  }

  KnitroInterface::KnitroInterface(const Function& nlp) : NlpSolverInternal(nlp) {
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
    NlpSolverInternal::init();

    //if (hasSetOption("Alg")) int_param_["alg"] = getOption("Alg");
    if (hasSetOption("BarRule")) int_param_["barrule"] = getOption("BarRule");
    if (hasSetOption("NewPoint")) int_param_["newpoint"] = getOption("NewPoint");
    if (hasSetOption("GradOpt")) int_param_["gradopt"] = getOption("GradOpt");
    if (hasSetOption("HessOpt")) int_param_["hessopt"] = getOption("HessOpt");
    if (hasSetOption("Feasible")) int_param_["feasible"] = getOption("Feasible");
    if (hasSetOption("HonorBnds")) int_param_["honorbnds"] = getOption("HonorBnds");
    if (hasSetOption("LpSolver")) int_param_["lpsolver"] = getOption("LpSolver");
    if (hasSetOption("Multistart")) int_param_["multistart"] = getOption("Multistart");
    //if (hasSetOption("MsMaxSolves")) int_param_["msmaxsolves"] = getOption("MsMaxSolves");
    if (hasSetOption("MaxCgIt")) int_param_["maxcgit"] = getOption("MaxCgIt");
    //if (hasSetOption("MaxCrossTt")) int_param_["maxcrosstt"] = getOption("MaxCrossTt");
    if (hasSetOption("MaxIt")) int_param_["maxit"] = getOption("MaxIt");
    //if (hasSetOption("MaxTimeCPU")) double_param_["maxtimecpu"] = getOption("MaxTimeCPU");
    //if (hasSetOption("MaxTimeReal")) double_param_["maxtimereal"] = getOption("MaxTimeReal");
    if (hasSetOption("LmSize")) int_param_["lmsize"] = getOption("LmSize");
    if (hasSetOption("Scale")) int_param_["scale"] = getOption("Scale");
    if (hasSetOption("ShiftInit")) int_param_["shiftinit"] = getOption("ShiftInit");
    if (hasSetOption("Soc")) int_param_["soc"] = getOption("Soc");
    if (hasSetOption("InitPt")) int_param_["initpt"] = getOption("InitPt");
    if (hasSetOption("Delta")) double_param_["delta"] = getOption("Delta");
    if (hasSetOption("FeasModeTol")) double_param_["feasmodetol"] = getOption("FeasModeTol");
    if (hasSetOption("FeasTol")) double_param_["feastol"] = getOption("FeasTol");
    if (hasSetOption("FeasTolAbs")) double_param_["feastolabs"] = getOption("FeasTolAbs");
    if (hasSetOption("OptTol")) double_param_["opttol"] = getOption("OptTol");
    if (hasSetOption("OptTolAbs")) double_param_["opttolabs"] = getOption("OptTolAbs");
    if (hasSetOption("Pivot")) double_param_["pivot"] = getOption("Pivot");
    if (hasSetOption("XTol")) double_param_["xtol"] = getOption("XTol");
    if (hasSetOption("Mu")) double_param_["mu"] = getOption("Mu");
    if (hasSetOption("ObjRange")) double_param_["objrange"] = getOption("ObjRange");
    if (hasSetOption("OutLev")) int_param_["outlev"] = getOption("OutLev");
    if (hasSetOption("Debug")) int_param_["debug"] = getOption("Debug");

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

  }

  void KnitroInterface::evaluate() {
    // Allocate KNITRO memory block (move back to init!)
    casadi_assert(kc_handle_==0);
    kc_handle_ = KTR_new();
    casadi_assert(kc_handle_!=0);
    int status;

    // Jacobian sparsity
    vector<int> Jcol, Jrow;
    if (!jacG_.isNull()) {
      Jcol = jacG_.output().sparsity().getCol();
      Jrow = jacG_.output().row();
    }

    // Hessian sparsity
    int nnzH = hessLag_.isNull() ? 0 : hessLag_.output().sizeL();
    vector<int> Hcol(nnzH), Hrow(nnzH);
    if (nnzH>0) {
      const vector<int> &colind = hessLag_.output().colind();
      const vector<int> &row = hessLag_.output().row();
      int nz=0;
      for (int cc=0; cc<colind.size()-1; ++cc) {
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
    for (std::map<std::string, double>::iterator it=double_param_.begin();
        it!=double_param_.end(); ++it) {
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
      vector<int> contype = getOption("contype");
      casadi_assert(contype.size()==cType.size());
      copy(contype.begin(), contype.end(), cType.begin());
    }

    // "Correct" upper and lower bounds
    for (vector<double>::iterator it=input(NLP_SOLVER_LBX).begin();
        it!=input(NLP_SOLVER_LBX).end(); ++it)
      if (isinf(*it)) *it = -KTR_INFBOUND;
    for (vector<double>::iterator it=input(NLP_SOLVER_UBX).begin();
        it!=input(NLP_SOLVER_UBX).end(); ++it)
      if (isinf(*it)) *it =  KTR_INFBOUND;
    for (vector<double>::iterator it=input(NLP_SOLVER_LBG).begin();
        it!=input(NLP_SOLVER_LBG).end(); ++it)
      if (isinf(*it)) *it = -KTR_INFBOUND;
    for (vector<double>::iterator it=input(NLP_SOLVER_UBG).begin();
        it!=input(NLP_SOLVER_UBG).end(); ++it)
      if (isinf(*it)) *it =  KTR_INFBOUND;

    // Initialize KNITRO
    status = KTR_init_problem(kc_handle_,
                              nx_,
                              KTR_OBJGOAL_MINIMIZE,
                              KTR_OBJTYPE_GENERAL,
                              input(NLP_SOLVER_LBX).ptr(),
                              input(NLP_SOLVER_UBX).ptr(),
                              ng_,
                              getPtr(cType),
                              input(NLP_SOLVER_LBG).ptr(),
                              input(NLP_SOLVER_UBG).ptr(),
                              Jcol.size(),
                              getPtr(Jcol),
                              getPtr(Jrow),
                              nnzH,
                              getPtr(Hrow),
                              getPtr(Hcol),
                              input(NLP_SOLVER_X0).ptr(),
                              0); // initial lambda
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
    status = KTR_solve(kc_handle_,
                       output(NLP_SOLVER_X).ptr(),
                       getPtr(lambda),
                       0,  // not used
                       output(NLP_SOLVER_F).ptr(),
                       0,  // not used
                       0,  // not used
                       0,  // not used
                       0,  // not used
                       0,  // not used
                       static_cast<void*>(this)); // to be retrieved in the callback function
    stats_["return_status"] = status;

    // Copy constraints
    nlp_.output(NL_G).get(output(NLP_SOLVER_G));

    // Copy lagrange multipliers
    output(NLP_SOLVER_LAM_G).set(getPtr(lambda));
    output(NLP_SOLVER_LAM_X).set(getPtr(lambda)+ng_);

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
      cerr << "KnitroInterface::callback caugth exception: " << ex.what() << endl;
      return -1;
    }
  }

  void KnitroInterface::evalfc(const double* x, double& obj, double *c) {
    // Pass the argument to the function
    nlp_.setInput(x, NL_X);
    nlp_.setInput(input(NLP_SOLVER_P), NL_P);

    // Evaluate the function
    nlp_.evaluate();

    // Get the result
    nlp_.output(NL_F).get(obj);
    nlp_.output(NL_G).get(c, DENSE);

    // Printing
    if (monitored("eval_f")) {
      cout << "x = " << nlp_.input(NL_X) << endl;
      cout << "f = " << nlp_.output(NL_F) << endl;
    }
    if (monitored("eval_g")) {
      cout << "x = " << nlp_.input(NL_X) << endl;
      cout << "g = " << nlp_.output(NL_G) << endl;
    }
  }

  void KnitroInterface::evalga(const double* x, double* objGrad, double* jac) {
    // Pass the argument to the function
    gradF_.setInput(x, NL_X);
    gradF_.setInput(input(NLP_SOLVER_P), NL_P);

    // Evaluate the function using adjoint mode AD
    gradF_.evaluate();

    // Get the result
    gradF_.output().get(objGrad, DENSE);

    // Printing
    if (monitored("eval_grad_f")) {
      cout << "x = " << gradF_.input(NL_X) << endl;
      cout << "grad_f = " << gradF_.output() << endl;
    }

    if (!jacG_.isNull()) {
      // Pass the argument to the Jacobian function
      jacG_.setInput(x, NL_X);
      jacG_.setInput(input(NLP_SOLVER_P), NL_P);

      // Evaluate the Jacobian function
      jacG_.evaluate();

      // Get the result
      jacG_.output().get(jac);

      // Printing
      if (monitored("eval_jac_g")) {
        cout << "x = " << jacG_.input(NL_X) << endl;
        cout << "jac_g = " << jacG_.output() << endl;
      }
    }
  }

  void KnitroInterface::evalh(const double* x, const double* lambda, double* hessian) {
    // Pass the argument to the function
    hessLag_.setInput(x, NL_X);
    hessLag_.setInput(input(NLP_SOLVER_P), NL_P);
    hessLag_.setInput(1.0, NL_NUM_IN+NL_F);
    hessLag_.setInput(lambda, NL_NUM_IN+NL_G);

    // Evaluate
    hessLag_.evaluate();

    // Get results
    hessLag_.output().get(hessian, SPARSESYM);

    // Printing
    if (monitored("eval_h")) {
      cout << "eval_h" << endl;
      cout << "x = " << hessLag_.input(0) << endl;
      cout << "lambda = " << hessLag_.input(1) << endl;
      cout << "scale = " << hessLag_.input(2) << endl;
      cout << "H = " << hessLag_ << endl;
    }
  }

} // namespace casadi
