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


#include "idas_interface.hpp"
#include "casadi/core/std_vector_tools.hpp"
#include "casadi/core/function/linear_solver_internal.hpp"
#include "casadi/core/function/mx_function.hpp"
#include "casadi/core/sx/sx_tools.hpp"
#include "casadi/core/mx/mx_tools.hpp"

#ifdef WITH_SYSTEM_SUNDIALS
#include <external_packages/sundials-2.5mod/idas/idas_spils_impl.h>
#endif

using namespace std;
namespace casadi {

  extern "C"
  int CASADI_INTEGRATOR_IDAS_EXPORT
      casadi_register_integrator_idas(IntegratorInternal::Plugin* plugin) {
    plugin->creator = IdasInterface::creator;
    plugin->name = "idas";
    plugin->doc = IdasInterface::meta_doc.c_str();
    plugin->version = 21;
    return 0;
  }

  extern "C"
  void CASADI_INTEGRATOR_IDAS_EXPORT casadi_load_integrator_idas() {
    IntegratorInternal::registerPlugin(casadi_register_integrator_idas);
  }

  IdasInterface* IdasInterface::clone() const {
    // Return a deep copy
    IdasInterface* node = new IdasInterface(f_, g_);
    node->setOption(dictionary());
    node->jac_ = jac_;
    node->linsol_ = linsol_;
    return node;
  }

  void IdasInterface::deepCopyMembers(std::map<SharedObjectNode*, SharedObject>& already_copied) {
    SundialsInterface::deepCopyMembers(already_copied);
  }

  IdasInterface::IdasInterface(const Function& f, const Function& g) : SundialsInterface(f, g) {
    addOption("suppress_algebraic",          OT_BOOLEAN,          false,
              "Suppress algebraic variables in the error testing");
    addOption("calc_ic",                     OT_BOOLEAN,          true,
              "Use IDACalcIC to get consistent initial conditions.");
    addOption("calc_icB",                    OT_BOOLEAN,          GenericType(),
              "Use IDACalcIC to get consistent initial conditions for "
              "backwards system [default: equal to calc_ic].");
    addOption("abstolv",                     OT_REALVECTOR);
    addOption("fsens_abstolv",               OT_REALVECTOR);
    addOption("max_step_size",               OT_REAL,             0,
              "Maximim step size");
    addOption("first_time",                  OT_REAL,             GenericType(),
              "First requested time as a fraction of the time interval");
    addOption("cj_scaling",                  OT_BOOLEAN,          false,
              "IDAS scaling on cj for the user-defined linear solver module");
    addOption("extra_fsens_calc_ic",         OT_BOOLEAN,          false,
              "Call calc ic an extra time, with fsens=0");
    addOption("disable_internal_warnings",   OT_BOOLEAN,          false,
              "Disable IDAS internal warning messages");
    addOption("monitor",                     OT_STRINGVECTOR,     GenericType(), "",
              "correctInitialConditions|res|resS|resB|rhsQB|bjacB|jtimesB|psetupB|psolveB|psetup",
              true);
    addOption("init_xdot",                   OT_REALVECTOR,       GenericType(),
              "Initial values for the state derivatives");

    mem_ = 0;

    xz_  = 0;
    xzdot_ = 0,
      q_ = 0;

    rxz_ = 0;
    rxzdot_ = 0;
    rq_ = 0;

    isInitAdj_ = false;
    isInitTaping_ = false;
    disable_internal_warnings_ = false;
  }

  IdasInterface::~IdasInterface() {
    freeIDAS();
  }

  void IdasInterface::freeIDAS() {
    if (mem_) { IDAFree(&mem_); mem_ = 0; }

    // Forward integration
    if (xz_) { N_VDestroy_Serial(xz_); xz_ = 0; }
    if (xzdot_) { N_VDestroy_Serial(xzdot_); xzdot_ = 0; }
    if (q_) { N_VDestroy_Serial(q_); q_ = 0; }

    // Backward integration
    if (rxz_) { N_VDestroy_Serial(rxz_); rxz_ = 0; }
    if (rxzdot_) { N_VDestroy_Serial(rxzdot_); rxzdot_ = 0; }
    if (rq_) { N_VDestroy_Serial(rq_); rq_ = 0; }

    // Forward problem
    for (vector<N_Vector>::iterator it=xzF_.begin(); it != xzF_.end(); ++it)
        if (*it) { N_VDestroy_Serial(*it); *it = 0; }
    for (vector<N_Vector>::iterator it=xzdotF_.begin(); it != xzdotF_.end(); ++it)
        if (*it) { N_VDestroy_Serial(*it); *it = 0; }
    for (vector<N_Vector>::iterator it=qF_.begin(); it != qF_.end(); ++it)
        if (*it) { N_VDestroy_Serial(*it); *it = 0; }
  }

  void IdasInterface::init() {
    log("IdasInterface::init", "begin");

    // Free memory if already initialized
    if (isInit()) freeIDAS();

    // Call the base class init
    SundialsInterface::init();

    // Get initial conditions for the state derivatives
    if (hasSetOption("init_xdot") && !getOption("init_xdot").isNull()) {
      init_xdot_ = getOption("init_xdot").toDoubleVector();
      casadi_assert_message(
        init_xdot_.size()==nx_,
        "Option \"init_xdot\" has incorrect length. Expecting " << nx_
        << ", but got " << init_xdot_.size()
        << ". Note that this message may actually be generated by the augmented"
        " integrator. In that case, make use of the 'augmented_options' options"
        " to correct 'init_xdot' for the augmented integrator.");
    } else {
      init_xdot_.resize(nx_);
      fill(init_xdot_.begin(), init_xdot_.end(), 0);
    }

    // Read options
    cj_scaling_ = getOption("cj_scaling");
    calc_ic_ = getOption("calc_ic");
    calc_icB_ = hasSetOption("calc_icB") ?  getOption("calc_icB") : getOption("calc_ic");

    // Sundials return flag
    int flag;

    // Create IDAS memory block
    mem_ = IDACreate();
    if (mem_==0) throw CasadiException("IDACreate(): Creation failed");

    // Allocate n-vectors for ivp
    xz_ = N_VNew_Serial(nx_+nz_);
    xzdot_ = N_VNew_Serial(nx_+nz_);

    // Initialize Idas
    double t0 = 0;
    N_VConst(0.0, xz_);
    N_VConst(0.0, xzdot_);
    IDAInit(mem_, res_wrapper, t0, xz_, xzdot_);
    log("IdasInterface::init", "IDA initialized");

    // Disable internal warning messages?
    disable_internal_warnings_ = getOption("disable_internal_warnings");

    // Set error handler function
    flag = IDASetErrHandlerFn(mem_, ehfun_wrapper, this);
    casadi_assert_message(flag == IDA_SUCCESS, "IDASetErrHandlerFn");

    // Include algebraic variables in error testing
    flag = IDASetSuppressAlg(mem_, getOption("suppress_algebraic").toInt());
    casadi_assert_message(flag == IDA_SUCCESS, "IDASetSuppressAlg");

    // Maxinum order for the multistep method
    flag = IDASetMaxOrd(mem_, getOption("max_multistep_order").toInt());
    casadi_assert_message(flag == IDA_SUCCESS, "IDASetMaxOrd");

    // Set user data
    flag = IDASetUserData(mem_, this);
    casadi_assert_message(flag == IDA_SUCCESS, "IDASetUserData");

    // Set maximum step size
    flag = IDASetMaxStep(mem_, getOption("max_step_size").toDouble());
    casadi_assert_message(flag == IDA_SUCCESS, "IDASetMaxStep");

    if (hasSetOption("abstolv")) {
      // Vector absolute tolerances
      vector<double> abstolv = getOption("abstolv").toDoubleVector();
      N_Vector nv_abstol = N_VMake_Serial(abstolv.size(), getPtr(abstolv));
      flag = IDASVtolerances(mem_, reltol_, nv_abstol);
      casadi_assert_message(flag == IDA_SUCCESS, "IDASVtolerances");
      N_VDestroy_Serial(nv_abstol);
    } else {
      // Scalar absolute tolerances
      flag = IDASStolerances(mem_, reltol_, abstol_);
      casadi_assert_message(flag == IDA_SUCCESS, "IDASStolerances");
    }

    // Maximum number of steps
    IDASetMaxNumSteps(mem_, getOption("max_num_steps").toInt());
    if (flag != IDA_SUCCESS) idas_error("IDASetMaxNumSteps", flag);

    // Set algebraic components
    N_Vector id = N_VNew_Serial(nx_+nz_);
    fill_n(NV_DATA_S(id), nx_, 1);
    fill_n(NV_DATA_S(id)+nx_, nz_, 0);

    // Pass this information to IDAS
    flag = IDASetId(mem_, id);
    if (flag != IDA_SUCCESS) idas_error("IDASetId", flag);

    // Delete the allocated memory
    N_VDestroy_Serial(id);

    // attach a linear solver
    switch (linsol_f_) {
    case SD_DENSE:
      initDenseLinearSolver();
      break;
    case SD_BANDED:
      initBandedLinearSolver();
      break;
    case SD_ITERATIVE:
      initIterativeLinearSolver();
      break;
    case SD_USER_DEFINED:
      initUserDefinedLinearSolver();
      break;
    default: casadi_error("Uncaught switch");
    }

    // Quadrature equations
    if (nq_>0) {

      // Allocate n-vectors for quadratures
      q_ = N_VMake_Serial(nq_, &output(INTEGRATOR_QF).front());

      // Initialize quadratures in IDAS
      N_VConst(0.0, q_);
      flag = IDAQuadInit(mem_, rhsQ_wrapper, q_);
      if (flag != IDA_SUCCESS) idas_error("IDAQuadInit", flag);

      // Should the quadrature errors be used for step size control?
      if (getOption("quad_err_con").toInt()) {
        flag = IDASetQuadErrCon(mem_, true);
        casadi_assert_message(flag == IDA_SUCCESS, "IDASetQuadErrCon");

        // Quadrature error tolerances
        // TODO(Joel): vector absolute tolerances
        flag = IDAQuadSStolerances(mem_, reltol_, abstol_);
        if (flag != IDA_SUCCESS) idas_error("IDAQuadSStolerances", flag);
      }
    }

    log("IdasInterface::init", "attached linear solver");

    // // Forward sensitivity problem
    // if (nfdir_>0) {
    //   // Allocate n-vectors
    //   xzF_.resize(nfdir_, 0);
    //   xzdotF_.resize(nfdir_, 0);
    //   qF_.resize(nfdir_, 0);
    //   for (int i=0; i<nfdir_; ++i) {
    //     xzF_[i] = N_VNew_Serial(nx_+nz_);
    //     xzdotF_[i] = N_VNew_Serial(nx_+nz_);
    //     if (nq_>0) {
    //       qF_[i] = N_VMake_Serial(nq_, &fwdSens(INTEGRATOR_QF, i).front());
    //     }
    //   }

    //   // Get the sensitivity method
    //   if (getOption("sensitivity_method")== "simultaneous") ism_ = IDA_SIMULTANEOUS;
    //   else if (getOption("sensitivity_method")=="staggered") ism_ = IDA_STAGGERED;
    //   else throw CasadiException("IDAS: Unknown sensitivity method");

    //   // Copy the forward seeds
    //   for (int i=0; i<nfdir_; ++i) {
    //     const Matrix<double>& x = fwdSeed(INTEGRATOR_X0, i);
    //     copy(x.begin(), x.begin()+nx_, NV_DATA_S(xzF_[i]));
    //     N_VConst(0.0, xzdotF_[i]);
    //   }

    //   // Initialize forward sensitivities
    //   if (finite_difference_fsens_) {

    //     // Use finite differences to calculate the residual in the forward sensitivity equations
    //     flag = IDASensInit(mem_, nfdir_, ism_, 0, getPtr(xzF_), getPtr(xzdotF_));
    //     if (flag != IDA_SUCCESS) idas_error("IDASensInit", flag);

    //     // Scaling factors
    //     double* pbar = 0;
    //     vector<double> pbar_vec;
    //     if (hasSetOption("fsens_scaling_factors")) {
    //       pbar_vec = getOption("fsens_scaling_factors").toDoubleVector();
    //       pbar = getPtr(pbar_vec);
    //     }

    //     // Which parameters should be used to estimate the sensitivity equations
    //     int * plist = 0;
    //     vector<int> plist_vec;
    //     if (hasSetOption("fsens_sensitivity_parameters")) {
    //       plist_vec = getOption("fsens_sensitiviy_parameters").toIntVector();
    //       plist = getPtr(plist_vec);
    //     }

    //     // Specify parameters
    //     flag = IDASetSensParams(mem_, &input(INTEGRATOR_P).front(), pbar, plist);
    //     if (flag != IDA_SUCCESS) idas_error("IDASetSensParams", flag);

    //     //  IDASetSensDQMethod
    //     // ?

    //   } else {
    //     // Use AD to calculate the residual in the forward sensitivity equations
    //     flag = IDASensInit(mem_, nfdir_, ism_, resS_wrapper,&xzF_.front(), &xzdotF_.front());
    //     if (flag != IDA_SUCCESS) idas_error("IDASensInit", flag);
    //   }

    //   // IDAS bugfix
    //   IDAMem IDA_mem = IDAMem(mem_);
    //   int max_multistep_order = getOption("max_multistep_order");
    //   for (int i=0; i<=max_multistep_order; ++i) {
    //     for (int iS=0; iS<nfdir_; ++iS) {
    //       N_VConst(0.0, IDA_mem->ida_phiS[i][iS]);
    //     }
    //   }

    //   // Vector absolute tolerances
    //   vector<double> fsens_abstol(nfdir_, fsens_abstol_);

    //   // Set tolerances
    //   if (hasSetOption("fsens_abstolv")) {
    //     // quick hack
    //     vector<double> fsens_abstolv = getOption("fsens_abstolv");
    //     N_Vector nv_abstol = N_VMake_Serial(fsens_abstolv.size(), getPtr(fsens_abstolv));
    //     vector<N_Vector> nv_abstol_all(nfdir_, nv_abstol);
    //     flag = IDASensSVtolerances(mem_, fsens_reltol_, getPtr(nv_abstol_all));
    //     if (flag != IDA_SUCCESS) idas_error("IDASensSVtolerances", flag);
    //     N_VDestroy_Serial(nv_abstol);
    //   } else {
    //     flag = IDASensSStolerances(mem_, fsens_reltol_, getPtr(fsens_abstol));
    //     if (flag != IDA_SUCCESS) idas_error("IDASensSStolerances", flag);
    //   }

    //   // Set optional inputs
    //   bool errconS = getOption("fsens_err_con");
    //   flag = IDASetSensErrCon(mem_, errconS);
    //   if (flag != IDA_SUCCESS) idas_error("IDASetSensErrCon", flag);

    //   // Quadrature equations
    //   if (nq_>0) {
    //     for (vector<N_Vector>::iterator it=qF_.begin(); it!=qF_.end(); ++it) N_VConst(0.0, *it);
    //     flag = IDAQuadSensInit(mem_, rhsQS_wrapper, getPtr(qF_));
    //     if (flag != IDA_SUCCESS) idas_error("IDAQuadSensInit", flag);

    //     // Set tolerances
    //     flag = IDAQuadSensSStolerances(mem_, fsens_reltol_, getPtr(fsens_abstol));
    //     if (flag != IDA_SUCCESS) idas_error("IDAQuadSensSStolerances", flag);
    //   }

    //   log("IdasInterface::init", "initialized forward sensitivities");
    // } // enable fsens

    // Adjoint sensitivity problem
    if (!g_.isNull()) {

      // Allocate n-vectors
      rxz_ = N_VNew_Serial(nrx_+nrz_);
      rxzdot_ = N_VNew_Serial(nrx_+nrz_);
      N_VConst(0.0, rxz_);
      N_VConst(0.0, rxzdot_);

      // Allocate n-vectors for quadratures
      rq_ = N_VMake_Serial(nrq_, output(INTEGRATOR_RQF).ptr());
    }
    log("IdasInterface::init", "initialized adjoint sensitivities");

    isInitTaping_ = false;
    isInitAdj_ = false;
    log("IdasInterface::init", "end");
  }

  void IdasInterface::initTaping() {
    casadi_assert(!isInitTaping_);
    int flag;

    // Get the number of steos per checkpoint
    int Nd = getOption("steps_per_checkpoint");

    // Get the interpolation type
    int interpType;
    if (getOption("interpolation_type")=="hermite")
      interpType = IDA_HERMITE;
    else if (getOption("interpolation_type")=="polynomial")
      interpType = IDA_POLYNOMIAL;
    else
      throw CasadiException("\"interpolation_type\" must be \"hermite\" or \"polynomial\"");

    // Initialize adjoint sensitivities
    flag = IDAAdjInit(mem_, Nd, interpType);
    if (flag != IDA_SUCCESS) idas_error("IDAAdjInit", flag);

    isInitTaping_ = true;
  }

  void IdasInterface::initAdj() {
    log("IdasInterface::initAdj", "start");

    casadi_assert(!isInitAdj_);
    int flag;

    // Create backward problem
    flag = IDACreateB(mem_, &whichB_);
    if (flag != IDA_SUCCESS) idas_error("IDACreateB", flag);

    // Initialize the backward problem
    double tB0 = tf_;
    flag = IDAInitB(mem_, whichB_, resB_wrapper, tB0, rxz_, rxzdot_);
    if (flag != IDA_SUCCESS) idas_error("IDAInitB", flag);

    // Set tolerances
    flag = IDASStolerancesB(mem_, whichB_, reltolB_, abstolB_);
    if (flag!=IDA_SUCCESS) idas_error("IDASStolerancesB", flag);

    // User data
    flag = IDASetUserDataB(mem_, whichB_, this);
    if (flag != IDA_SUCCESS) idas_error("IDASetUserDataB", flag);

    // Maximum number of steps
    IDASetMaxNumStepsB(mem_, whichB_, getOption("max_num_steps").toInt());
    if (flag != IDA_SUCCESS) idas_error("IDASetMaxNumStepsB", flag);

    // Set algebraic components
    N_Vector id = N_VNew_Serial(nrx_+nrz_);
    fill_n(NV_DATA_S(id), nrx_, 1);
    fill_n(NV_DATA_S(id)+nrx_, nrz_, 0);

    // Pass this information to IDAS
    flag = IDASetIdB(mem_, whichB_, id);
    if (flag != IDA_SUCCESS) idas_error("IDASetIdB", flag);

    // Delete the allocated memory
    N_VDestroy_Serial(id);

    // attach linear solver
    switch (linsol_g_) {
    case SD_DENSE:
      initDenseLinearSolverB();
      break;
    case SD_BANDED:
      initBandedLinearSolverB();
      break;
    case SD_ITERATIVE:
      initIterativeLinearSolverB();
      break;
    case SD_USER_DEFINED:
      initUserDefinedLinearSolverB();
      break;
    default: casadi_error("Uncaught switch");
    }

    // Quadratures for the adjoint problem
    N_VConst(0.0, rq_);
    flag = IDAQuadInitB(mem_, whichB_, rhsQB_wrapper, rq_);
    if (flag!=IDA_SUCCESS) idas_error("IDAQuadInitB", flag);

    // Quadrature error control
    if (getOption("quad_err_con").toInt()) {
      flag = IDASetQuadErrConB(mem_, whichB_, true);
      if (flag != IDA_SUCCESS) idas_error("IDASetQuadErrConB", flag);

      flag = IDAQuadSStolerancesB(mem_, whichB_, reltolB_, abstolB_);
      if (flag != IDA_SUCCESS) idas_error("IDAQuadSStolerancesB", flag);
    }

    // Mark initialized
    isInitAdj_ = true;

    log("IdasInterface::initAdj", "end");
  }


  void IdasInterface::res(double t, const double* xz, const double* xzdot, double* r) {
    log("IdasInterface::res", "begin");

    // Get time
    time1 = clock();

    // Pass input
    f_.setInput(&t, DAE_T);
    f_.setInput(xz, DAE_X);
    f_.setInput(xz+nx_, DAE_Z);
    f_.setInput(input(INTEGRATOR_P), DAE_P);

    if (monitored("res")) {
      cout << "DAE_T    = " << t << endl;
      cout << "DAE_X    = " << f_.input(DAE_X) << endl;
      cout << "DAE_Z    = " << f_.input(DAE_Z) << endl;
      cout << "DAE_P    = " << f_.input(DAE_P) << endl;
    }

    // Evaluate
    f_.evaluate();

    // Get results
    f_.getOutput(r, DAE_ODE);
    f_.getOutput(r+nx_, DAE_ALG);

    if (monitored("res")) {
      cout << "ODE rhs  = " << f_.output(DAE_ODE) << endl;
      cout << "ALG rhs  = " << f_.output(DAE_ALG) << endl;
    }

    if (regularity_check_) {
      casadi_assert_message(isRegular(f_.output(DAE_ODE).data()),
                            "IdasInterface::res: f.output(DAE_ODE) is not regular.");
      casadi_assert_message(isRegular(f_.output(DAE_ALG).data()),
                            "IdasInterface::res: f.output(DAE_ALG) is not regular.");
    }

    // Subtract state derivative to get residual
    for (int i=0; i<nx_; ++i) {
      r[i] -= xzdot[i];
    }

    time2 = clock();
    t_res += static_cast<double>(time2-time1)/CLOCKS_PER_SEC;
    log("IdasInterface::res", "end");
  }

  int IdasInterface::res_wrapper(double t, N_Vector xz, N_Vector xzdot,
                                N_Vector rr, void *user_data) {
    try {
      IdasInterface *this_ = static_cast<IdasInterface*>(user_data);
      this_->res(t, NV_DATA_S(xz), NV_DATA_S(xzdot), NV_DATA_S(rr));
      return 0;
    } catch(int flag) { // recoverable error
      return flag;
    } catch(exception& e) { // non-recoverable error
      cerr << "res failed: " << e.what() << endl;
      return -1;
    }
  }

  void IdasInterface::ehfun_wrapper(int error_code, const char *module, const char *function,
                                   char *msg, void *eh_data) {
    try {
      IdasInterface *this_ = static_cast<IdasInterface*>(eh_data);
      this_->ehfun(error_code, module, function, msg);
    } catch(exception& e) {
      cerr << "ehfun failed: " << e.what() << endl;
    }
  }

  void IdasInterface::ehfun(int error_code, const char *module, const char *function, char *msg) {
    cerr << msg << endl;
  }

  void IdasInterface::jtimes(double t, const double *xz, const double *xzdot, const double *rr,
                            const double *v, double *Jv, double cj,
                            double *tmp1, double *tmp2) {
    log("IdasInterface::jtimes", "begin");
    // Get time
    time1 = clock();

    // Pass input
    f_fwd_.setInput(&t,                  DAE_T);
    f_fwd_.setInput(xz,                  DAE_X);
    f_fwd_.setInput(xz+nx_,              DAE_Z);
    f_fwd_.setInput(input(INTEGRATOR_P), DAE_P);

    // Pass seeds of the state vectors
    f_fwd_.setInput(v,      DAE_NUM_IN + DAE_X);
    f_fwd_.setInput(v+nx_,  DAE_NUM_IN + DAE_Z);

    // Evaluate the AD forward algorithm
    f_fwd_.evaluate();

    // Get the output seeds
    f_fwd_.getOutput(Jv,    DAE_NUM_OUT + DAE_ODE);
    f_fwd_.getOutput(Jv+nx_, DAE_NUM_OUT + DAE_ALG);

    // Subtract state derivative to get residual
    for (int i=0; i<nx_; ++i) {
      Jv[i] -= cj*v[i];
    }

    // Log time duration
    time2 = clock();
    t_jac += static_cast<double>(time2-time1)/CLOCKS_PER_SEC;
    log("IdasInterface::jtimes", "end");
  }

  int IdasInterface::jtimes_wrapper(double t, N_Vector xz, N_Vector xzdot, N_Vector rr, N_Vector v,
                                   N_Vector Jv, double cj, void *user_data,
                                   N_Vector tmp1, N_Vector tmp2) {
    try {
      IdasInterface *this_ = static_cast<IdasInterface*>(user_data);
      this_->jtimes(t, NV_DATA_S(xz), NV_DATA_S(xzdot), NV_DATA_S(rr), NV_DATA_S(v),
                    NV_DATA_S(Jv), cj, NV_DATA_S(tmp1), NV_DATA_S(tmp2));
      return 0;
    } catch(exception& e) {
      cerr << "jtimes failed: " << e.what() << endl;
      return 1;
    }
  }

  void IdasInterface::jtimesB(double t, const double *xz, const double *xzdot, const double *xzB,
                             const double *xzdotB, const double *resvalB, const double *vB,
                             double *JvB, double cjB, double * tmp1B, double * tmp2B) {
    log("IdasInterface::jtimesB", "begin");
    // Get time
    time1 = clock();

    // Pass input
    g_fwd_.setInput(&t,                  RDAE_T);
    g_fwd_.setInput(xz,                  RDAE_X);
    g_fwd_.setInput(xz+nx_,              RDAE_Z);
    g_fwd_.setInput(input(INTEGRATOR_P), RDAE_P);

    g_fwd_.setInput(xzB,                 RDAE_RX);
    g_fwd_.setInput(xzB+nrx_,            RDAE_RZ);
    g_fwd_.setInput(input(INTEGRATOR_RP), RDAE_RP);

    // Pass seeds of the state vectors
    g_fwd_.setInput(0.0,    RDAE_NUM_IN + RDAE_T);
    g_fwd_.setInput(0.0,    RDAE_NUM_IN + RDAE_X);
    g_fwd_.setInput(0.0,    RDAE_NUM_IN + RDAE_Z);
    g_fwd_.setInput(0.0,    RDAE_NUM_IN + RDAE_P);
    g_fwd_.setInput(vB,     RDAE_NUM_IN + RDAE_RX);
    g_fwd_.setInput(vB+nx_, RDAE_NUM_IN + RDAE_RZ);
    g_fwd_.setInput(0.0,    RDAE_NUM_IN + RDAE_RP);

    if (monitored("jtimesB")) {
      cout << "RDAE_T    = " << t << endl;
      cout << "RDAE_X    = " << g_fwd_.input(RDAE_X) << endl;
      cout << "RDAE_Z    = " << g_fwd_.input(RDAE_Z) << endl;
      cout << "RDAE_P    = " << g_fwd_.input(RDAE_P) << endl;
      cout << "RDAE_XDOT  = ";
      for (int k=0;k<nx_;++k) {
        cout << xzdot[k] << " " ;
      }
      cout << endl;
      cout << "RDAE_RX    = " << g_fwd_.input(RDAE_RX) << endl;
      cout << "RDAE_RZ    = " << g_fwd_.input(RDAE_RZ) << endl;
      cout << "RDAE_RP    = " << g_fwd_.input(RDAE_RP) << endl;
      cout << "RDAE_RXDOT  = ";
      for (int k=0;k<nrx_;++k) {
        cout << xzdotB[k] << " " ;
      }
      cout << endl;
      cout << "fwdSeed(RDAE_RX) = " << g_fwd_.input(RDAE_NUM_IN + RDAE_RX) << endl;
      cout << "fwdSeed(RDAE_RZ) = " << g_fwd_.input(RDAE_NUM_IN + RDAE_RZ) << endl;
    }

    // Evaluate the AD forward algorithm
    g_fwd_.evaluate();

    if (monitored("jtimesB")) {
      cout << "fwdSens(RDAE_ODE) = " << g_fwd_.output(RDAE_NUM_OUT + RDAE_ODE) << endl;
      cout << "fwdSens(RDAE_ALG) = " << g_fwd_.output(RDAE_NUM_OUT + RDAE_ALG) << endl;
    }

    // Get the output seeds
    g_fwd_.getOutput(JvB,     RDAE_NUM_OUT + RDAE_ODE);
    g_fwd_.getOutput(JvB+nx_, RDAE_NUM_OUT + RDAE_ALG);

    // Subtract state derivative to get residual
    for (int i=0; i<nrx_; ++i) {
      JvB[i] += cjB*vB[i];
    }

    if (monitored("jtimesB")) {
      g_fwd_.setOutput(JvB,     RDAE_NUM_OUT + RDAE_ODE);
      g_fwd_.setOutput(JvB+nx_, RDAE_NUM_OUT + RDAE_ALG);
      cout << "res fwdSens(RDAE_ODE)    = " << g_fwd_.output(RDAE_NUM_OUT + RDAE_ODE) << endl;
      cout << "res fwdSens(RDAE_ALG)    = " << g_fwd_.output(RDAE_NUM_OUT + RDAE_ALG) << endl;
    }

    // Log time duration
    time2 = clock();
    t_jac += static_cast<double>(time2-time1)/CLOCKS_PER_SEC;
    log("IdasInterface::jtimesB", "end");
  }

  int IdasInterface::jtimesB_wrapper(double t, N_Vector xz, N_Vector xzdot, N_Vector xzB,
                                    N_Vector xzdotB, N_Vector resvalB, N_Vector vB, N_Vector JvB,
                                    double cjB, void *user_data,
                                    N_Vector tmp1B, N_Vector tmp2B) {
    try {
      IdasInterface *this_ = static_cast<IdasInterface*>(user_data);
      this_->jtimesB(t, NV_DATA_S(xz), NV_DATA_S(xzdot), NV_DATA_S(xzB),
                     NV_DATA_S(xzdotB), NV_DATA_S(resvalB), NV_DATA_S(vB),
                     NV_DATA_S(JvB), cjB, NV_DATA_S(tmp1B), NV_DATA_S(tmp2B));
      return 0;
    } catch(exception& e) {
      cerr << "jtimesB failed: " << e.what() << endl;
      return 1;
    }
  }

  void IdasInterface::resS(int Ns, double t, const double* xz, const double* xzdot,
                          const double *resval, N_Vector *xzF, N_Vector* xzdotF, N_Vector *rrF,
                          double *tmp1, double *tmp2, double *tmp3) {
    log("IdasInterface::resS", "begin");
    //    casadi_assert(Ns==nfdir_);

    // Record the current cpu time
    time1 = clock();

    // Commented out since a new implementation currently cannot be tested
    casadi_error("Commented out, #884, #794.");

    // // Pass input
    // f_.setInput(&t, DAE_T);
    // f_.setInput(xz, DAE_X);
    // f_.setInput(xz+nx_, DAE_Z);
    // f_.setInput(input(INTEGRATOR_P), DAE_P);

    // // Calculate the forward sensitivities, nfdir_f_ directions at a time
    // for (int offset=0; offset<nfdir_; offset += nfdir_f_) {
    //   // Number of directions in this batch
    //   int nfdir_batch = std::min(nfdir_-offset, nfdir_f_);
    //   for (int dir=0; dir<nfdir_batch; ++dir) {
    //     // Pass forward seeds
    //     f_.fwdSeed(DAE_T, dir).setZero();
    //     f_.setFwdSeed(NV_DATA_S(xzF[offset+dir]), DAE_X, dir);
    //     f_.setFwdSeed(NV_DATA_S(xzF[offset+dir])+nx_, DAE_Z, dir);
    //     f_.setFwdSeed(fwdSeed(INTEGRATOR_P, offset+dir), DAE_P, dir);
    //   }

    //   // Evaluate the AD forward algorithm
    //   f_.evaluateOld(nfdir_batch, 0);

    //   // Get the output seeds
    //   for (int dir=0; dir<nfdir_batch; ++dir) {
    //     f_.getFwdSens(NV_DATA_S(rrF[offset+dir]), DAE_ODE, dir);
    //     f_.getFwdSens(NV_DATA_S(rrF[offset+dir])+nx_, DAE_ALG, dir);

    //     // Subtract state derivative to get residual
    //     for (int i=0; i<nx_; ++i) {
    //       NV_DATA_S(rrF[offset+dir])[i] -= NV_DATA_S(xzdotF[offset+dir])[i];
    //     }
    //   }
    // }

    // Record timings
    time2 = clock();
    t_fres += static_cast<double>(time2-time1)/CLOCKS_PER_SEC;
    log("IdasInterface::resS", "end");
  }

  int IdasInterface::resS_wrapper(int Ns, double t, N_Vector xz, N_Vector xzdot, N_Vector resval,
                                 N_Vector *xzF, N_Vector *xzdotF, N_Vector *rrF, void *user_data,
                                 N_Vector tmp1, N_Vector tmp2, N_Vector tmp3) {
    try {
      IdasInterface *this_ = static_cast<IdasInterface*>(user_data);
      this_->resS(Ns, t, NV_DATA_S(xz), NV_DATA_S(xzdot), NV_DATA_S(resval),
                  xzF, xzdotF, rrF, NV_DATA_S(tmp1), NV_DATA_S(tmp2), NV_DATA_S(tmp3));
      return 0;
    } catch(exception& e) {
      cerr << "resS failed: " << e.what() << endl;
      return 1;
    }
  }

  void IdasInterface::reset() {
    log("IdasInterface::reset", "begin");

    // Reset the base classes
    SundialsInterface::reset();

    if (nrx_>0 && !isInitTaping_)
      initTaping();

    // If we have forward sensitivities, rest one extra time without forward sensitivities to get
    // a consistent initial guess
    //   if (nfdir>0 && getOption("extra_fsens_calc_ic").toInt())
    //     reset(0);

    // Reset timers
    t_res = t_fres = t_jac = t_jacB = t_lsolve = t_lsetup_jac = t_lsetup_fac = 0;

    // Return flag
    int flag;

    // Copy to N_Vectors
    const Matrix<double>& x = output(INTEGRATOR_XF);
    copy(x.begin(), x.begin()+nx_, NV_DATA_S(xz_));
    const Matrix<double>& z = output(INTEGRATOR_ZF);
    copy(z.begin(), z.end(), NV_DATA_S(xz_)+nx_);
    copy(init_xdot_.begin(), init_xdot_.end(), NV_DATA_S(xzdot_));

    // Re-initialize
    flag = IDAReInit(mem_, t0_, xz_, xzdot_);
    if (flag != IDA_SUCCESS) idas_error("IDAReInit", flag);
    log("IdasInterface::reset", "re-initialized IVP solution");


    // Re-initialize quadratures
    if (nq_>0) {
      N_VConst(0.0, q_);
      flag = IDAQuadReInit(mem_, q_);
      if (flag != IDA_SUCCESS) idas_error("IDAQuadReInit", flag);
      log("IdasInterface::reset", "re-initialized quadratures");
    }

    // if (nsens>0) {
    //   // Get the forward seeds
    //   for (int i=0; i<nfdir_; ++i) {
    //     const Matrix<double>& x0_seed = fwdSeed(INTEGRATOR_X0, i);
    //     copy(x0_seed.begin(), x0_seed.begin()+nx_, NV_DATA_S(xzF_[i]));
    //     N_VConst(0.0, xzdotF_[i]);
    //   }

    //   // Re-initialize sensitivities
    //   flag = IDASensReInit(mem_, ism_, getPtr(xzF_), getPtr(xzdotF_));
    //   if (flag != IDA_SUCCESS) idas_error("IDASensReInit", flag);
    //   log("IdasInterface::reset", "re-initialized forward sensitivity solution");

    //   if (nq_>0) {
    //     for (vector<N_Vector>::iterator it=qF_.begin(); it!=qF_.end(); ++it) N_VConst(0.0, *it);
    //     flag = IDAQuadSensReInit(mem_, getPtr(qF_));
    //     if (flag != IDA_SUCCESS) idas_error("IDAQuadSensReInit", flag);
    //     log("IdasInterface::reset", "re-initialized forward sensitivity dependent quadratures");
    //   }
    // } else {
    // Turn of sensitivities
    flag = IDASensToggleOff(mem_);
    if (flag != IDA_SUCCESS) idas_error("IDASensToggleOff", flag);
    // }

    // Correct initial conditions, if necessary
    if (calc_ic_) {
      correctInitialConditions();
    }

    // Re-initialize backward integration
    if (nrx_>0) {
      flag = IDAAdjReInit(mem_);
      if (flag != IDA_SUCCESS) idas_error("IDAAdjReInit", flag);
    }

    // Set the stop time of the integration -- don't integrate past this point
    if (stop_at_end_) setStopTime(tf_);

    log("IdasInterface::reset", "end");
  }


  void IdasInterface::correctInitialConditions() {
    log("IdasInterface::correctInitialConditions", "begin");
    if (monitored("correctInitialConditions")) {
      cout << "initial guess: " << endl;
      cout << "p = " << input(INTEGRATOR_P) << endl;
      cout << "x0 = " << input(INTEGRATOR_X0) << endl;
      // if (nsens_>0) {
      //   for (int dir=0; dir<nfdir_; ++dir) {
      //     cout << "forward seed guess, direction " << dir << ": " << endl;
      //     cout << "p_seed = " << fwdSeed(INTEGRATOR_P, dir) << endl;
      //     cout << "x0_seed = " << fwdSeed(INTEGRATOR_X0, dir) << endl;
      //   }
      // }
    }

    int icopt = IDA_YA_YDP_INIT; // calculate z and xdot given x
    // int icopt = IDA_Y_INIT; // calculate z and x given zdot and xdot (e.g. start in stationary)

    double t_first =
        hasSetOption("first_time") ? static_cast<double>(getOption("first_time")) : tf_;
    int flag = IDACalcIC(mem_, icopt , t_first);
    if (flag != IDA_SUCCESS) idas_error("IDACalcIC", flag);

    // Retrieve the initial values
    flag = IDAGetConsistentIC(mem_, xz_, xzdot_);
    if (flag != IDA_SUCCESS) idas_error("IDAGetConsistentIC", flag);

    // Print progress
    log("IdasInterface::correctInitialConditions", "found consistent initial values");
    if (monitored("correctInitialConditions")) {
      cout << "p = " << input(INTEGRATOR_P) << endl;
      cout << "x0 = " << input(INTEGRATOR_X0) << endl;
      // if (nsens_>0) {
      //   for (int dir=0; dir<nfdir_; ++dir) {
      //     cout << "forward seed, direction " << dir << ": " << endl;
      //     cout << "p_seed = " << fwdSeed(INTEGRATOR_P, dir) << endl;
      //     cout << "x0_seed = " << fwdSeed(INTEGRATOR_X0, dir) << endl;
      //   }
      // }
    }
    log("IdasInterface::correctInitialConditions", "end");
  }

  void IdasInterface::integrate(double t_out) {
    casadi_log("IdasInterface::integrate(" << t_out << ") begin");

    casadi_assert_message(t_out>=t0_, "IdasInterface::integrate(" << t_out << "): "
                          "Cannot integrate to a time earlier than t0 (" << t0_ << ")");
    casadi_assert_message(t_out<=tf_ || !stop_at_end_, "IdasInterface::integrate(" << t_out << "): "
                          "Cannot integrate past a time later than tf (" << tf_ << ") "
                          "unless stop_at_end is set to False.");

    int flag;

    // Check if we are already at the output time
    double ttol = 1e-9;   // tolerance
    if (fabs(t_-t_out)<ttol) {
      // No integration necessary
      log("IdasInterface::integrate", "already at the end of the horizon end");

    } else {
      // Integrate ...
      if (nrx_>0) {
        // ... with taping
        log("IdasInterface::integrate", "integration with taping");
        flag = IDASolveF(mem_, t_out, &t_, xz_, xzdot_, IDA_NORMAL, &ncheck_);
        if (flag != IDA_SUCCESS && flag != IDA_TSTOP_RETURN) idas_error("IDASolveF", flag);
      } else {
        // ... without taping
        log("IdasInterface::integrate", "integration without taping");
        flag = IDASolve(mem_, t_out, &t_, xz_, xzdot_, IDA_NORMAL);
        if (flag != IDA_SUCCESS && flag != IDA_TSTOP_RETURN) idas_error("IDASolve", flag);
      }
      log("IdasInterface::integrate", "integration complete");

      // Get quadrature states
      if (nq_>0) {
        double tret;
        flag = IDAGetQuad(mem_, &tret, q_);
        if (flag != IDA_SUCCESS) idas_error("IDAGetQuad", flag);
      }

      // if (nsens_>0) {
      //   // Get the sensitivities
      //   flag = IDAGetSens(mem_, &t_, getPtr(xzF_));
      //   if (flag != IDA_SUCCESS) idas_error("IDAGetSens", flag);

      //   if (nq_>0) {
      //     double tret;
      //     flag = IDAGetQuadSens(mem_, &tret, getPtr(qF_));
      //     if (flag != IDA_SUCCESS) idas_error("IDAGetQuadSens", flag);
      //   }
      // }
    }

    // Save the final state
    copy(NV_DATA_S(xz_), NV_DATA_S(xz_)+nx_, output(INTEGRATOR_XF).begin());
    // if (nsens_>0) {
    //   for (int i=0; i<nfdir_; ++i) {
    //     copy(NV_DATA_S(xzF_[i]), NV_DATA_S(xzF_[i])+nx_, fwdSens(INTEGRATOR_XF, i).begin());
    //   }
    // }

    // Save the final algebraic variable
    copy(NV_DATA_S(xz_)+nx_, NV_DATA_S(xz_)+nx_+nz_, output(INTEGRATOR_ZF).begin());

    // Print statistics
    if (getOption("print_stats")) printStats(std::cout);

    if (gather_stats_) {
      long nsteps, nfevals, nlinsetups, netfails;
      int qlast, qcur;
      double hinused, hlast, hcur, tcur;
      int flag = IDAGetIntegratorStats(mem_, &nsteps, &nfevals, &nlinsetups, &netfails,
                                       &qlast, &qcur, &hinused, &hlast, &hcur, &tcur);
      if (flag!=IDA_SUCCESS) idas_error("IDAGetIntegratorStats", flag);

      stats_["nsteps"] = 1.0*nsteps;
      stats_["nlinsetups"] = 1.0*nlinsetups;

    }


    casadi_log("IdasInterface::integrate(" << t_out << ") end");
  }

  void IdasInterface::resetB() {
    log("IdasInterface::resetB", "begin");

    int flag;

    // Reset adjoint sensitivities for the parameters
    N_VConst(0.0, rq_);

    // Get the backward state
    const Matrix<double> &rx0 = input(INTEGRATOR_RX0);
    copy(rx0.begin(), rx0.end(), NV_DATA_S(rxz_));

    if (isInitAdj_) {
      flag = IDAReInitB(mem_, whichB_, tf_, rxz_, rxzdot_);
      if (flag != IDA_SUCCESS) idas_error("IDAReInitB", flag);

      if (nrq_>0) {
        N_VConst(0.0, rq_);
        flag = IDAQuadReInit(IDAGetAdjIDABmem(mem_, whichB_), rq_);
        // flag = IDAQuadReInitB(mem_, whichB_[dir], rq_[dir]); // BUG in Sundials
        //                                                      // do not use this!
      }
      if (flag!=IDA_SUCCESS) idas_error("IDAQuadReInitB", flag);
    } else {
      // Initialize the adjoint integration
      initAdj();
    }

    // Correct initial values for the integration if necessary
    if (calc_icB_) {
      log("IdasInterface::resetB", "IDACalcICB begin");
      flag = IDACalcICB(mem_, whichB_, t0_, xz_, xzdot_);
      if (flag != IDA_SUCCESS) idas_error("IDACalcICB", flag);
      log("IdasInterface::resetB", "IDACalcICB end");

      // Retrieve the initial values
      flag = IDAGetConsistentICB(mem_, whichB_, rxz_, rxzdot_);
      if (flag != IDA_SUCCESS) idas_error("IDAGetConsistentICB", flag);

    }

    log("IdasInterface::resetB", "end");

  }

  void IdasInterface::integrateB(double t_out) {
    casadi_log("IdasInterface::integrateB(" << t_out << ") begin");
    int flag;
    // Integrate backwards to t_out
    flag = IDASolveB(mem_, t_out, IDA_NORMAL);
    if (flag<IDA_SUCCESS) idas_error("IDASolveB", flag);

    // Get the sensitivities
    double tret;
    flag = IDAGetB(mem_, whichB_, &tret, rxz_, rxzdot_);
    if (flag!=IDA_SUCCESS) idas_error("IDAGetB", flag);

    if (nrq_>0) {
      flag = IDAGetQuadB(mem_, whichB_, &tret, rq_);
      if (flag!=IDA_SUCCESS) idas_error("IDAGetQuadB", flag);
    }

    // Save the backward state and algebraic variable
    const double *rxz = NV_DATA_S(rxz_);
    copy(rxz, rxz+nrx_, output(INTEGRATOR_RXF).begin());
    copy(rxz+nrx_, rxz+nrx_+nrz_, output(INTEGRATOR_RZF).begin());

    if (gather_stats_) {
      long nsteps, nfevals, nlinsetups, netfails;
      int qlast, qcur;
      double hinused, hlast, hcur, tcur;

      IDAMem IDA_mem = IDAMem(mem_);
      IDAadjMem IDAADJ_mem = IDA_mem->ida_adj_mem;
      IDABMem IDAB_mem = IDAADJ_mem->IDAB_mem;

      int flag = IDAGetIntegratorStats(IDAB_mem->IDA_mem, &nsteps, &nfevals, &nlinsetups,
                                       &netfails, &qlast, &qcur, &hinused, &hlast, &hcur, &tcur);
      if (flag!=IDA_SUCCESS) idas_error("IDAGetIntegratorStatsB", flag);

      stats_["nstepsB"] = 1.0*nsteps;
      stats_["nlinsetupsB"] = 1.0*nlinsetups;
    }
    casadi_log("IdasInterface::integrateB(" << t_out << ") end");
  }

  void IdasInterface::printStats(std::ostream &stream) const {
    long nsteps, nfevals, nlinsetups, netfails;
    int qlast, qcur;
    double hinused, hlast, hcur, tcur;
    int flag = IDAGetIntegratorStats(mem_, &nsteps, &nfevals, &nlinsetups, &netfails, &qlast,
                                     &qcur, &hinused, &hlast, &hcur, &tcur);
    if (flag!=IDA_SUCCESS) idas_error("IDAGetIntegratorStats", flag);

    // Get the number of right hand side evaluations in the linear solver
    long nfevals_linsol=0;
    switch (linsol_f_) {
    case SD_DENSE:
    case SD_BANDED:
      flag = IDADlsGetNumResEvals(mem_, &nfevals_linsol);
      if (flag!=IDA_SUCCESS) idas_error("IDADlsGetNumResEvals", flag);
      break;
    case SD_ITERATIVE:
      flag = IDASpilsGetNumResEvals(mem_, &nfevals_linsol);
      if (flag!=IDA_SUCCESS) idas_error("IDASpilsGetNumResEvals", flag);
      break;
    default:
      nfevals_linsol = 0;
    }

    stream << "number of steps taken by IDAS:            " << nsteps << std::endl;
    stream << "number of calls to the user's f function: " << (nfevals + nfevals_linsol)
           << std::endl;
    stream << "   step calculation:                      " << nfevals << std::endl;
    stream << "   linear solver:                         " << nfevals_linsol << std::endl;
    stream << "number of calls made to the linear solver setup function: " << nlinsetups
           << std::endl;
    stream << "number of error test failures: " << netfails << std::endl;
    stream << "method order used on the last internal step: " << qlast << std::endl;
    stream << "method order to be used on the next internal step: " << qcur << std::endl;
    stream << "actual value of initial step size: " << hinused << std::endl;
    stream << "step size taken on the last internal step: " << hlast << std::endl;
    stream << "step size to be attempted on the next internal step: " << hcur << std::endl;
    stream << "current internal time reached: " << tcur << std::endl;
    stream << std::endl;

    stream << "number of checkpoints stored: " << ncheck_ << endl;
    stream << std::endl;

    stream << "Time spent in the DAE residual: " << t_res << " s." << endl;
    stream << "Time spent in the forward sensitivity residual: " << t_fres << " s." << endl;
    stream << "Time spent in the jacobian function or jacobian times vector function: "
           << t_jac << " s." << endl;
    stream << "Time spent in the linear solver solve function: " << t_lsolve << " s." << endl;
    stream << "Time spent to generate the jacobian in the linear solver setup function: "
           << t_lsetup_jac << " s." << endl;
    stream << "Time spent to factorize the jacobian in the linear solver setup function: "
           << t_lsetup_fac << " s." << endl;
    stream << std::endl;
  }

  void IdasInterface::idas_error(const string& module, int flag) {
    // Find the error
    char* flagname = IDAGetReturnFlagName(flag);
    stringstream ss;
    ss << "Module \"" << module << "\" returned flag " << flag << " (\"" << flagname << "\").";
    ss << " Consult Idas documentation." << std::endl;
    free(flagname);

    // Heuristics
    if (
        (module=="IDACalcIC" && (flag==IDA_CONV_FAIL || flag==IDA_NO_RECOVERY ||
                                 flag==IDA_LINESEARCH_FAIL)) ||
        (module=="IDASolve" && flag ==IDA_ERR_FAIL)
        ) {
      ss << "Some common causes for this error: " << std::endl;
      ss << "  - providing an initial guess for which 0=g(y, z, t) is not invertible wrt y. "
         << std::endl;
      ss << "  - having a DAE-index higher than 1 such that 0=g(y, z, t) is not invertible wrt y "
          "over the whole domain." << std::endl;
      ss << "  - having set abstol or reltol too small." << std::endl;
      ss << "  - using 'calcic'=True for systems that are not semi-explicit index-one. "
          "You must provide consistent initial conditions yourself in this case. " << std::endl;
      ss << "  - your problem is too hard for IDAcalcIC to solve. Provide consistent "
          "initial conditions yourself." << std::endl;
    }

    casadi_error(ss.str());
  }

  int IdasInterface::rhsQ_wrapper(double t, N_Vector xz, N_Vector xzdot, N_Vector rhsQ,
                                 void *user_data) {
    try {
      IdasInterface *this_ = static_cast<IdasInterface*>(user_data);
      this_->rhsQ(t, NV_DATA_S(xz), NV_DATA_S(xzdot), NV_DATA_S(rhsQ));
      return 0;
    } catch(exception& e) {
      cerr << "rhsQ failed: " << e.what() << endl;
      return 1;
    }
  }

  void IdasInterface::rhsQ(double t, const double* xz, const double* xzdot, double* rhsQ) {
    log("IdasInterface::rhsQ", "begin");
    // Pass input
    f_.setInput(&t, DAE_T);
    f_.setInput(xz, DAE_X);
    f_.setInput(xz+nx_, DAE_Z);
    f_.setInput(input(INTEGRATOR_P), DAE_P);

    // Evaluate
    f_.evaluate();

    // Get results
    f_.getOutput(rhsQ, DAE_QUAD);
    log("IdasInterface::rhsQ", "end");
  }

  void IdasInterface::rhsQS(int Ns, double t, N_Vector xz, N_Vector xzdot, N_Vector *xzF,
                           N_Vector *xzdotF, N_Vector rrQ, N_Vector *qdotF,
                           N_Vector tmp1, N_Vector tmp2, N_Vector tmp3) {

    log("IdasInterface::rhsQS", "enter");
    //    casadi_assert(Ns==nfdir_);

    // Commented out since a new implementation currently cannot be tested
    casadi_error("Commented out, #884, #794.");

    // // Pass input
    //  f_.setInput(&t, DAE_T);
    //  f_.setInput(NV_DATA_S(xz), DAE_X);
    //  f_.setInput(NV_DATA_S(xz)+nx_, DAE_Z);
    //  f_.setInput(input(INTEGRATOR_P), DAE_P);

    //  // Pass forward seeds
    // for (int i=0; i<nfdir_; ++i) {
    //   f_.fwdSeed(DAE_T).setZero();
    //   f_.setFwdSeed(NV_DATA_S(xzF[i]), DAE_X);
    //   f_.setFwdSeed(NV_DATA_S(xzF[i])+nx_, DAE_Z);
    //   f_.setFwdSeed(fwdSeed(INTEGRATOR_P, i), DAE_P);

    //   // Evaluate the AD forward algorithm
    //   f_.evaluateOld(1, 0);

    //   // Get the output seeds
    //   f_.getFwdSens(NV_DATA_S(qdotF[i]), DAE_QUAD);
    //  }
    log("IdasInterface::rhsQS", "end");
  }

  int IdasInterface::rhsQS_wrapper(int Ns, double t, N_Vector xz, N_Vector xzdot, N_Vector *xzF,
                                  N_Vector *xzdotF, N_Vector rrQ, N_Vector *qdotF, void *user_data,
                                  N_Vector tmp1, N_Vector tmp2, N_Vector tmp3) {

    try {
      IdasInterface *this_ = static_cast<IdasInterface*>(user_data);
      this_->rhsQS(Ns, t, xz, xzdot, xzF, xzdotF, rrQ, qdotF, tmp1, tmp2, tmp3);
      return 0;
    } catch(exception& e) {
      cerr << "rhsQS failed: " << e.what() << endl;
      return 1;
    }
  }

  void IdasInterface::resB(double t, const double* xz, const double* xzdot, const double* xzA,
                          const double* xzdotA, double* rrA) {
    log("IdasInterface::resB", "begin");

    // Pass inputs
    g_.setInput(&t, RDAE_T);
    g_.setInput(xz, RDAE_X);
    g_.setInput(xz+nx_, RDAE_Z);
    g_.setInput(input(INTEGRATOR_P), RDAE_P);
    g_.setInput(input(INTEGRATOR_RP), RDAE_RP);
    g_.setInput(xzA, RDAE_RX);
    g_.setInput(xzA+nrx_, RDAE_RZ);

    if (monitored("resB")) {
      cout << "RDAE_T    = " << t << endl;
      cout << "RDAE_X    = " << g_.input(RDAE_X) << endl;
      cout << "RDAE_Z    = " << g_.input(RDAE_Z) << endl;
      cout << "RDAE_P    = " << g_.input(RDAE_P) << endl;
      cout << "RDAE_XDOT  = ";
      for (int k=0;k<nx_;++k) {
        cout << xzdot[k] << " " ;
      }
      cout << endl;
      cout << "RDAE_RX    = " << g_.input(RDAE_RX) << endl;
      cout << "RDAE_RZ    = " << g_.input(RDAE_RZ) << endl;
      cout << "RDAE_RP    = " << g_.input(RDAE_RP) << endl;
      cout << "RDAE_RXDOT  = ";
      for (int k=0;k<nrx_;++k) {
        cout << xzdotA[k] << " " ;
      }
      cout << endl;
    }

    // Evaluate
    g_.evaluate();

    // Save to output
    g_.getOutput(rrA, RDAE_ODE);
    g_.getOutput(rrA+nrx_, RDAE_ALG);

    if (monitored("resB")) {
      cout << "RDAE_ODE    = " << g_.output(RDAE_ODE) << endl;
      cout << "RDAE_ALG    = " << g_.output(RDAE_ALG) << endl;
    }

    // Add state derivative to get residual (note definition of g)
    for (int i=0; i<nrx_; ++i) {
      rrA[i] += xzdotA[i];
    }

    if (monitored("resB")) {
      g_.setOutput(rrA, RDAE_ODE);
      g_.setOutput(rrA+nrx_, RDAE_ALG);
      cout << "res ODE    = " << g_.output(RDAE_ODE) << endl;
      cout << "res ALG    = " << g_.output(RDAE_ALG) << endl;
    }


    log("IdasInterface::resB", "end");
  }

  int IdasInterface::resB_wrapper(double t, N_Vector xz, N_Vector xzdot, N_Vector xzA,
                                 N_Vector xzdotA, N_Vector rrA, void *user_data) {
    try {
      IdasInterface *this_ = static_cast<IdasInterface*>(user_data);
      this_->resB(t, NV_DATA_S(xz), NV_DATA_S(xzdot), NV_DATA_S(xzA),
                  NV_DATA_S(xzdotA), NV_DATA_S(rrA));
      return 0;
    } catch(exception& e) {
      cerr << "resB failed: " << e.what() << endl;
      return 1;
    }
  }

  void IdasInterface::rhsQB(double t, const double* xz, const double* xzdot, const double* xzA,
                           const double* xzdotA, double *qdotA) {
    log("IdasInterface::rhsQB", "begin");

    // Pass inputs
    g_.setInput(&t, RDAE_T);
    g_.setInput(xz, RDAE_X);
    g_.setInput(xz+nx_, RDAE_Z);
    g_.setInput(input(INTEGRATOR_P), RDAE_P);
    g_.setInput(input(INTEGRATOR_RP), RDAE_RP);
    g_.setInput(xzA, RDAE_RX);
    g_.setInput(xzA+nrx_, RDAE_RZ);

    // Evaluate
    g_.evaluate();

    // Save to output
    g_.getOutput(qdotA, RDAE_QUAD);

    if (monitored("rhsQB")) {
      cout << "RDAE_T    = " << t << endl;
      cout << "RDAE_X    = " << g_.input(RDAE_X) << endl;
      cout << "RDAE_Z    = " << g_.input(RDAE_Z) << endl;
      cout << "RDAE_P    = " << g_.input(RDAE_P) << endl;
      cout << "RDAE_RX    = " << g_.input(RDAE_RX) << endl;
      cout << "RDAE_RZ    = " << g_.input(RDAE_RZ) << endl;
      cout << "RDAE_RP    = " << g_.input(RDAE_RP) << endl;
      cout << "rhs = " << g_.output(RDAE_QUAD) << endl;
    }

    // Negate (note definition of g)
    for (int i=0; i<nrq_; ++i)
      qdotA[i] *= -1;

    log("IdasInterface::rhsQB", "end");
  }

  int IdasInterface::rhsQB_wrapper(double t, N_Vector y, N_Vector xzdot, N_Vector xzA,
                                  N_Vector xzdotA, N_Vector qdotA, void *user_data) {
    try {
      IdasInterface *this_ = static_cast<IdasInterface*>(user_data);
      this_->rhsQB(t, NV_DATA_S(y), NV_DATA_S(xzdot), NV_DATA_S(xzA),
                   NV_DATA_S(xzdotA), NV_DATA_S(qdotA));
      return 0;
    } catch(exception& e) {
      cerr << "rhsQB failed: " << e.what() << endl;
      return 1;
    }
  }

  void IdasInterface::djac(long Neq, double t, double cj, N_Vector xz, N_Vector xzdot,
                          N_Vector rr, DlsMat Jac,
                          N_Vector tmp1, N_Vector tmp2, N_Vector tmp3) {
    log("IdasInterface::djac", "begin");

    // Get time
    time1 = clock();

    // Pass input to the Jacobian function
    jac_.setInput(&t, DAE_T);
    jac_.setInput(NV_DATA_S(xz), DAE_X);
    jac_.setInput(NV_DATA_S(xz)+nx_, DAE_Z);
    jac_.setInput(input(INTEGRATOR_P), DAE_P);
    jac_.setInput(cj, DAE_NUM_IN);

    // Evaluate Jacobian
    jac_.evaluate();

    // Get sparsity and non-zero elements
    const vector<int>& colind = jac_.output().colind();
    const vector<int>& row = jac_.output().row();
    const vector<double>& val = jac_.output().data();

    // Loop over columns
    for (int cc=0; cc<colind.size()-1; ++cc) {
      // Loop over non-zero entries
      for (int el=colind[cc]; el<colind[cc+1]; ++el) {

        // Get row
        int rr = row[el];

        // Add to the element
        DENSE_ELEM(Jac, rr, cc) = val[el];
      }
    }

    // Log time duration
    time2 = clock();
    t_jac += static_cast<double>(time2-time1)/CLOCKS_PER_SEC;
    log("IdasInterface::djac", "end");
  }

  int IdasInterface::djac_wrapper(long Neq, double t, double cj, N_Vector xz, N_Vector xzdot,
                                 N_Vector rr, DlsMat Jac, void *user_data,
                                 N_Vector tmp1, N_Vector tmp2, N_Vector tmp3) {
    try {
      IdasInterface *this_ = static_cast<IdasInterface*>(user_data);
      this_->djac(Neq, t, cj, xz, xzdot, rr, Jac, tmp1, tmp2, tmp3);
      return 0;
    } catch(exception& e) {
      cerr << "djac failed: " << e.what() << endl;
      return 1;
    }
  }

  void IdasInterface::djacB(long int NeqB, double t, double cjB, N_Vector xz, N_Vector xzdot,
                           N_Vector xzB, N_Vector xzdotB, N_Vector rrB, DlsMat JacB,
                           N_Vector tmp1B, N_Vector tmp2B, N_Vector tmp3B) {
    log("IdasInterface::djacB", "begin");

    // Get time
    time1 = clock();

    // Pass input to the Jacobian function
    jacB_.setInput(&t, RDAE_T);
    jacB_.setInput(NV_DATA_S(xz), RDAE_X);
    jacB_.setInput(NV_DATA_S(xz)+nx_, RDAE_Z);
    jacB_.setInput(input(INTEGRATOR_P), RDAE_P);
    jacB_.setInput(input(INTEGRATOR_RP), RDAE_RP);
    jacB_.setInput(NV_DATA_S(xzB), RDAE_RX);
    jacB_.setInput(NV_DATA_S(xzB)+nrx_, RDAE_RZ);
    jacB_.setInput(cjB, RDAE_NUM_IN);

    if (monitored("djacB")) {
      cout << "RDAE_T    = " << t << endl;
      cout << "RDAE_X    = " << jacB_.input(RDAE_X) << endl;
      cout << "RDAE_Z    = " << jacB_.input(RDAE_Z) << endl;
      cout << "RDAE_P    = " << jacB_.input(RDAE_P) << endl;
      cout << "RDAE_XDOT  = ";
      for (int k=0;k<nx_;++k) {
        cout << NV_DATA_S(xzdot)[k] << " " ;
      }
      cout << endl;
      cout << "RDAE_RX    = " << jacB_.input(RDAE_RX) << endl;
      cout << "RDAE_RZ    = " << jacB_.input(RDAE_RZ) << endl;
      cout << "RDAE_RP    = " << jacB_.input(RDAE_RP) << endl;
      cout << "RDAE_RXDOT  = ";
      for (int k=0;k<nrx_;++k) {
        cout << NV_DATA_S(xzdotB)[k] << " " ;
      }
      cout << endl;
      cout << "cjB = " << cjB << endl;
    }

    // Evaluate Jacobian
    jacB_.evaluate();

    if (monitored("djacB")) {
      cout << "jacB = " << jacB_.output() << endl;
    }

    // Get sparsity and non-zero elements
    const vector<int>& colind = jacB_.output().colind();
    const vector<int>& row = jacB_.output().row();
    const vector<double>& val = jacB_.output().data();

    // Loop over columns
    for (int cc=0; cc<colind.size()-1; ++cc) {
      // Loop over non-zero entries
      for (int el=colind[cc]; el<colind[cc+1]; ++el) {

        // Get row
        int rr = row[el];

        // Add to the element
        DENSE_ELEM(JacB, rr, cc) = val[el];
      }
    }

    // Log time duration
    time2 = clock();
    t_jacB += static_cast<double>(time2-time1)/CLOCKS_PER_SEC;
    log("IdasInterface::djacB", "end");
  }

  int IdasInterface::djacB_wrapper(long int NeqB, double t, double cjB, N_Vector xz, N_Vector xzdot,
                                  N_Vector xzB, N_Vector xzdotB, N_Vector rrB, DlsMat JacB,
                                  void *user_data,
                                  N_Vector tmp1B, N_Vector tmp2B, N_Vector tmp3B) {
    try {
      IdasInterface *this_ = static_cast<IdasInterface*>(user_data);
      this_->djacB(NeqB, t, cjB, xz, xzdot, xzB, xzdotB, rrB, JacB, tmp1B, tmp2B, tmp3B);
      return 0;
    } catch(exception& e) {
      cerr << "djacB failed: " << e.what() << endl;
      return 1;
    }
  }

  void IdasInterface::bjac(long Neq, long mupper, long mlower, double t, double cj, N_Vector xz,
                          N_Vector xzdot, N_Vector rr, DlsMat Jac,
                          N_Vector tmp1, N_Vector tmp2, N_Vector tmp3) {
    log("IdasInterface::bjac", "begin");
    // Get time
    time1 = clock();

    // Pass input to the jacobian function
    jac_.setInput(&t, DAE_T);
    jac_.setInput(NV_DATA_S(xz), DAE_X);
    jac_.setInput(NV_DATA_S(xz)+nx_, DAE_Z);
    jac_.setInput(input(INTEGRATOR_P), DAE_P);
    jac_.setInput(cj, DAE_NUM_IN);

    // Evaluate jacobian
    jac_.evaluate();

    // Get sparsity and non-zero elements
    const vector<int>& colind = jac_.output().colind();
    const vector<int>& row = jac_.output().row();
    const vector<double>& val = jac_.output().data();

    // Loop over columns
    for (int cc=0; cc<colind.size()-1; ++cc) {
      // Loop over non-zero entries
      for (int el=colind[cc]; el<colind[cc+1]; ++el) {
        // Get row
        int rr = row[el];

        // Set the element
        if (cc-rr<=mupper && rr-cc<=mlower)
          BAND_ELEM(Jac, rr, cc) = val[el];
      }
    }

    // Log time duration
    time2 = clock();
    t_jac += static_cast<double>(time2-time1)/CLOCKS_PER_SEC;
    log("IdasInterface::bjac", "end");
  }

  int IdasInterface::bjac_wrapper(long Neq, long mupper, long mlower, double t, double cj,
                                 N_Vector xz, N_Vector xzdot, N_Vector rr,
                                 DlsMat Jac, void *user_data,
                                 N_Vector tmp1, N_Vector tmp2, N_Vector tmp3) {
    try {
      IdasInterface *this_ = static_cast<IdasInterface*>(user_data);
      this_->bjac(Neq, mupper, mlower, t, cj, xz, xzdot, rr, Jac, tmp1, tmp2, tmp3);
      return 0;
    } catch(exception& e) {
      cerr << "bjac failed: " << e.what() << endl;
      return 1;
    }
  }

  void IdasInterface::bjacB(long NeqB, long mupperB, long mlowerB, double t, double cjB,
                           N_Vector xz, N_Vector xzdot, N_Vector xzB, N_Vector xzdotB,
                           N_Vector resvalB, DlsMat JacB, N_Vector tmp1B, N_Vector tmp2B,
                           N_Vector tmp3B) {
    log("IdasInterface::bjacB", "begin");

    // Get time
    time1 = clock();

    // Pass input to the Jacobian function
    jacB_.setInput(&t, RDAE_T);
    jacB_.setInput(NV_DATA_S(xz), RDAE_X);
    jacB_.setInput(NV_DATA_S(xz)+nx_, RDAE_Z);
    jacB_.setInput(input(INTEGRATOR_P), RDAE_P);
    jacB_.setInput(input(INTEGRATOR_RP), RDAE_RP);
    jacB_.setInput(NV_DATA_S(xzB), RDAE_RX);
    jacB_.setInput(NV_DATA_S(xzB)+nrx_, RDAE_RZ);
    jacB_.setInput(cjB, RDAE_NUM_IN);

    if (monitored("bjacB")) {
      cout << "RDAE_T    = " << t << endl;
      cout << "RDAE_X    = " << jacB_.input(RDAE_X) << endl;
      cout << "RDAE_Z    = " << jacB_.input(RDAE_Z) << endl;
      cout << "RDAE_P    = " << jacB_.input(RDAE_P) << endl;
      cout << "RDAE_XDOT  = ";
      for (int k=0;k<nx_;++k) {
        cout << NV_DATA_S(xzdot)[k] << " " ;
      }
      cout << endl;
      cout << "RDAE_RX    = " << jacB_.input(RDAE_RX) << endl;
      cout << "RDAE_RZ    = " << jacB_.input(RDAE_RZ) << endl;
      cout << "RDAE_RP    = " << jacB_.input(RDAE_RP) << endl;
      cout << "RDAE_RXDOT  = ";
      for (int k=0;k<nrx_;++k) {
        cout << NV_DATA_S(xzdotB)[k] << " " ;
      }
      cout << endl;
      cout << "cjB = " << cjB << endl;
    }

    // Evaluate Jacobian
    jacB_.evaluate();

    if (monitored("bjacB")) {
      cout << "jacB = " << jacB_.output() << endl;
    }

    // Get sparsity and non-zero elements
    const vector<int>& colind = jacB_.output().colind();
    const vector<int>& row = jacB_.output().row();
    const vector<double>& val = jacB_.output().data();

    // Loop over columns
    for (int cc=0; cc<colind.size()-1; ++cc) {
      // Loop over non-zero entries
      for (int el=colind[cc]; el<colind[cc+1]; ++el) {
        // Get row
        int rr = row[el];

        // Set the element
        if (cc-rr<=mupperB && rr-cc<=mlowerB)
          BAND_ELEM(JacB, rr, cc) = val[el];
      }
    }

    // Log time duration
    time2 = clock();
    t_jacB += static_cast<double>(time2-time1)/CLOCKS_PER_SEC;
    log("IdasInterface::bjacB", "end");
  }

  int IdasInterface::bjacB_wrapper(
    long NeqB, long mupperB, long mlowerB, double t, double cjB,
    N_Vector xz, N_Vector xzdot, N_Vector xzB, N_Vector xzdotB,
    N_Vector resvalB, DlsMat JacB, void *user_data, N_Vector tmp1B,
    N_Vector tmp2B, N_Vector tmp3B) {
    try {
      IdasInterface *this_ = static_cast<IdasInterface*>(user_data);
      this_->bjacB(NeqB, mupperB, mlowerB, t, cjB, xz, xzdot, xzB, xzdotB,
                   resvalB, JacB, tmp1B, tmp2B, tmp3B);
      return 0;
    } catch(exception& e) {
      cerr << "bjacB failed: " << e.what() << endl;
      return 1;
    }
  }

  void IdasInterface::setStopTime(double tf) {
    // Set the stop time of the integration -- don't integrate past this point
    int flag = IDASetStopTime(mem_, tf);
    if (flag != IDA_SUCCESS) idas_error("IDASetStopTime", flag);
  }

  int IdasInterface::psolve_wrapper(double t, N_Vector xz, N_Vector xzdot, N_Vector rr,
                                   N_Vector rvec, N_Vector zvec, double cj, double delta,
                                   void *user_data, N_Vector tmp) {
    try {
      IdasInterface *this_ = static_cast<IdasInterface*>(user_data);
      casadi_assert(this_);
      this_->psolve(t, xz, xzdot, rr, rvec, zvec, cj, delta, tmp);
      return 0;
    } catch(exception& e) {
      cerr << "psolve failed: " << e.what() << endl;
      return 1;
    }
  }

  int IdasInterface::psolveB_wrapper(double t, N_Vector xz, N_Vector xzdot, N_Vector xzB,
                                    N_Vector xzdotB, N_Vector resvalB, N_Vector rvecB,
                                    N_Vector zvecB, double cjB, double deltaB,
                                    void *user_data, N_Vector tmpB) {
    try {
      IdasInterface *this_ = static_cast<IdasInterface*>(user_data);
      casadi_assert(this_);
      this_->psolveB(t, xz, xzdot, xzB, xzdotB, resvalB, rvecB, zvecB, cjB, deltaB, tmpB);
      return 0;
    } catch(exception& e) {
      cerr << "psolveB failed: " << e.what() << endl;
      return 1;
    }
  }

  int IdasInterface::psetup_wrapper(double t, N_Vector xz, N_Vector xzdot, N_Vector rr,
                                   double cj, void* user_data,
                                   N_Vector tmp1, N_Vector tmp2, N_Vector tmp3) {
    try {
      IdasInterface *this_ = static_cast<IdasInterface*>(user_data);
      casadi_assert(this_);
      this_->psetup(t, xz, xzdot, rr, cj, tmp1, tmp2, tmp3);
      return 0;
    } catch(exception& e) {
      cerr << "psetup failed: " << e.what() << endl;
      return 1;
    }
  }

  int IdasInterface::psetupB_wrapper(double t, N_Vector xz, N_Vector xzdot,
                                    N_Vector xzB, N_Vector xzdotB,
                                    N_Vector resvalB, double cjB, void *user_data,
                                    N_Vector tmp1B, N_Vector tmp2B, N_Vector tmp3B) {
    try {
      IdasInterface *this_ = static_cast<IdasInterface*>(user_data);
      casadi_assert(this_);
      this_->psetupB(t, xz, xzdot, xzB, xzdotB, resvalB, cjB, tmp1B, tmp2B, tmp3B);
      return 0;
    } catch(exception& e) {
      cerr << "psetupB failed: " << e.what() << endl;
      return 1;
    }
  }

  void IdasInterface::psolve(double t, N_Vector xz, N_Vector xzdot, N_Vector rr, N_Vector rvec,
                            N_Vector zvec, double cj, double delta, N_Vector tmp) {
    log("IdasInterface::psolve", "begin");

    // Get time
    time1 = clock();

    // Copy input to output, if necessary
    if (rvec!=zvec) {
      N_VScale(1.0, rvec, zvec);
    }

    // Solve the (possibly factorized) system
    casadi_assert_message(linsol_.output().size() == NV_LENGTH_S(zvec), "Assertion error: "
                          << linsol_.output().size() << " == " << NV_LENGTH_S(zvec));
    linsol_.solve(NV_DATA_S(zvec), 1, false);

    // Log time duration
    time2 = clock();
    t_lsolve += static_cast<double>(time2-time1)/CLOCKS_PER_SEC;
    log("IdasInterface::psolve", "end");
  }


  void IdasInterface::psolveB(double t, N_Vector xz, N_Vector xzdot, N_Vector xzB, N_Vector xzdotB,
                             N_Vector resvalB, N_Vector rvecB, N_Vector zvecB,
                             double cjB, double deltaB, N_Vector tmpB) {
    log("IdasInterface::psolveB", "begin");

    // Get time
    time1 = clock();

    // Copy input to output, if necessary
    if (rvecB!=zvecB) {
      N_VScale(1.0, rvecB, zvecB);
    }

    casadi_assert(!linsolB_.isNull());

    // Solve the (possibly factorized) system
    casadi_assert_message(linsolB_.output().size() == NV_LENGTH_S(zvecB),
                          "Assertion error: " << linsolB_.output().size()
                          << " == " << NV_LENGTH_S(zvecB));
    if (monitored("psolveB")) {
      cout << "zvecB = " << std::endl;
      for (int k=0;k<NV_LENGTH_S(zvecB);++k) {
        cout << NV_DATA_S(zvecB)[k] << " " ;
      }
      cout << endl;
    }

    linsolB_.solve(NV_DATA_S(zvecB), 1, false);

    if (monitored("psolveB")) {
      cout << "zvecB sol = " << std::endl;
      for (int k=0;k<NV_LENGTH_S(zvecB);++k) {
        cout << NV_DATA_S(zvecB)[k] << " " ;
      }
      cout << endl;
    }

    // Log time duration
    time2 = clock();
    t_lsolve += static_cast<double>(time2-time1)/CLOCKS_PER_SEC;
    log("IdasInterface::psolveB", "end");
  }

  void IdasInterface::psetup(double t, N_Vector xz, N_Vector xzdot, N_Vector rr, double cj,
                            N_Vector tmp1, N_Vector tmp2, N_Vector tmp3) {
    log("IdasInterface::psetup", "begin");

    // Get time
    time1 = clock();

    // Pass input to the jacobian function
    jac_.setInput(&t, DAE_T);
    jac_.setInput(NV_DATA_S(xz), DAE_X);
    jac_.setInput(NV_DATA_S(xz)+nx_, DAE_Z);
    jac_.setInput(input(INTEGRATOR_P), DAE_P);
    jac_.setInput(cj, DAE_NUM_IN);

    if (monitored("psetup")) {
      cout << "DAE_T    = " << t << endl;
      cout << "DAE_X    = " << jac_.input(DAE_X) << endl;
      cout << "DAE_Z    = " << jac_.input(DAE_Z) << endl;
      cout << "DAE_P    = " << jac_.input(DAE_P) << endl;
      cout << "cj = " << cj << endl;
    }

    // Evaluate jacobian
    jac_.evaluate();

    if (monitored("psetup")) {
      cout << "psetup = " << jac_.output() << endl;
    }

    // Log time duration
    time2 = clock();
    t_lsetup_jac += static_cast<double>(time2-time1)/CLOCKS_PER_SEC;

    // Pass non-zero elements to the linear solver
    linsol_.setInput(jac_.output(), 0);

    // Prepare the solution of the linear system (e.g. factorize)
    // -- only if the linear solver inherits from LinearSolver
    linsol_.prepare();

    // Log time duration
    time1 = clock();
    t_lsetup_fac += static_cast<double>(time1-time2)/CLOCKS_PER_SEC;

    log("IdasInterface::psetup", "end");

  }


  void IdasInterface::psetupB(double t, N_Vector xz, N_Vector xzdot, N_Vector xzB, N_Vector xzdotB,
                             N_Vector resvalB, double cjB,
                             N_Vector tmp1B, N_Vector tmp2B, N_Vector tmp3B) {
    log("IdasInterface::psetupB", "begin");

    // Get time
    time1 = clock();

    // Pass input to the Jacobian function
    jacB_.setInput(&t, RDAE_T);
    jacB_.setInput(NV_DATA_S(xz), RDAE_X);
    jacB_.setInput(NV_DATA_S(xz)+nx_, RDAE_Z);
    jacB_.setInput(input(INTEGRATOR_P), RDAE_P);
    jacB_.setInput(input(INTEGRATOR_RP), RDAE_RP);
    jacB_.setInput(NV_DATA_S(xzB), RDAE_RX);
    jacB_.setInput(NV_DATA_S(xzB)+nrx_, RDAE_RZ);
    jacB_.setInput(cjB, RDAE_NUM_IN);

    if (monitored("psetupB")) {
      cout << "RDAE_T    = " << t << endl;
      cout << "RDAE_X    = " << jacB_.input(RDAE_X) << endl;
      cout << "RDAE_Z    = " << jacB_.input(RDAE_Z) << endl;
      cout << "RDAE_P    = " << jacB_.input(RDAE_P) << endl;
      cout << "RDAE_XDOT  = ";
      for (int k=0;k<nx_;++k) {
        cout << NV_DATA_S(xzdot)[k] << " " ;
      }
      cout << endl;
      cout << "RDAE_RX    = " << jacB_.input(RDAE_RX) << endl;
      cout << "RDAE_RZ    = " << jacB_.input(RDAE_RZ) << endl;
      cout << "RDAE_RP    = " << jacB_.input(RDAE_RP) << endl;
      cout << "RDAE_RXDOT  = ";
      for (int k=0;k<nrx_;++k) {
        cout << NV_DATA_S(xzdotB)[k] << " " ;
      }
      cout << endl;
      cout << "cjB = " << cjB << endl;
    }

    // Evaluate jacobian
    jacB_.evaluate();

    if (monitored("psetupB")) {
      cout << "psetupB = " << jacB_.output() << endl;
    }

    // Log time duration
    time2 = clock();
    t_lsetup_jac += static_cast<double>(time2-time1)/CLOCKS_PER_SEC;

    // Pass non-zero elements to the linear solver
    linsolB_.setInput(jacB_.output(), 0);

    // Prepare the solution of the linear system (e.g. factorize)
    // -- only if the linear solver inherits from LinearSolver
    linsolB_.prepare();

    // Log time duration
    time1 = clock();
    t_lsetup_fac += static_cast<double>(time1-time2)/CLOCKS_PER_SEC;

    log("IdasInterface::psetupB", "end");

  }


  int IdasInterface::lsetup_wrapper(IDAMem IDA_mem, N_Vector xz, N_Vector xzdot, N_Vector resp,
                                   N_Vector vtemp1, N_Vector vtemp2, N_Vector vtemp3) {
    try {
      IdasInterface *this_ = static_cast<IdasInterface*>(IDA_mem->ida_lmem);
      casadi_assert(this_);
      this_->lsetup(IDA_mem, xz, xzdot, resp, vtemp1, vtemp2, vtemp3);
      return 0;
    } catch(exception& e) {
      cerr << "lsetup failed: " << e.what() << endl;
      return -1;
    }
  }

  int IdasInterface::lsetupB_wrapper(IDAMem IDA_mem, N_Vector xzB, N_Vector xzdotB, N_Vector respB,
                                    N_Vector vtemp1B, N_Vector vtemp2B, N_Vector vtemp3B) {
    try {
      IdasInterface *this_ = static_cast<IdasInterface*>(IDA_mem->ida_lmem);
      casadi_assert(this_);
      IDAadjMem IDAADJ_mem;
      //IDABMem IDAB_mem;
      int flag;

      // Current time
      double t = IDA_mem->ida_tn; // TODO(Joel): is this correct?
      // Multiple of df_dydot to be added to the matrix
      double cj = IDA_mem->ida_cj;

      IDA_mem = static_cast<IDAMem>(IDA_mem->ida_user_data);

      IDAADJ_mem = IDA_mem->ida_adj_mem;
      //IDAB_mem = IDAADJ_mem->ia_bckpbCrt;

      // Get FORWARD solution from interpolation.
      if (IDAADJ_mem->ia_noInterp==FALSE) {
        flag = IDAADJ_mem->ia_getY(IDA_mem, t, IDAADJ_mem->ia_yyTmp, IDAADJ_mem->ia_ypTmp,
                                   NULL, NULL);
        if (flag != IDA_SUCCESS) casadi_error("Could not interpolate forward states");
      }
      this_->lsetupB(t, cj, IDAADJ_mem->ia_yyTmp,
                     IDAADJ_mem->ia_ypTmp, xzB, xzdotB, respB, vtemp1B, vtemp2B, vtemp3B);
      return 0;
    } catch(exception& e) {
      cerr << "lsetupB failed: " << e.what() << endl;
      return -1;
    }
  }


  int IdasInterface::lsolve_wrapper(IDAMem IDA_mem, N_Vector b, N_Vector weight, N_Vector xz,
                                   N_Vector xzdot, N_Vector rr) {
    try {
      IdasInterface *this_ = static_cast<IdasInterface*>(IDA_mem->ida_lmem);
      casadi_assert(this_);
      this_->lsolve(IDA_mem, b, weight, xz, xzdot, rr);
      return 0;
    } catch(int wrn) {
      /*    cerr << "warning: " << wrn << endl;*/
      return wrn;
    } catch(exception& e) {
      cerr << "lsolve failed: " << e.what() << endl;
      return -1;
    }
  }

  int IdasInterface::lsolveB_wrapper(IDAMem IDA_mem, N_Vector b, N_Vector weight, N_Vector xzB,
                                    N_Vector xzdotB, N_Vector rrB) {
    try {
      IdasInterface *this_ = static_cast<IdasInterface*>(IDA_mem->ida_lmem);
      casadi_assert(this_);
      IDAadjMem IDAADJ_mem;
      //IDABMem IDAB_mem;
      int flag;

      // Current time
      double t = IDA_mem->ida_tn; // TODO(Joel): is this correct?
      // Multiple of df_dydot to be added to the matrix
      double cj = IDA_mem->ida_cj;
      double cjratio = IDA_mem->ida_cjratio;

      IDA_mem = (IDAMem) IDA_mem->ida_user_data;

      IDAADJ_mem = IDA_mem->ida_adj_mem;
      //IDAB_mem = IDAADJ_mem->ia_bckpbCrt;

      // Get FORWARD solution from interpolation.
      if (IDAADJ_mem->ia_noInterp==FALSE) {
        flag = IDAADJ_mem->ia_getY(IDA_mem, t, IDAADJ_mem->ia_yyTmp, IDAADJ_mem->ia_ypTmp,
                                   NULL, NULL);
        if (flag != IDA_SUCCESS) casadi_error("Could not interpolate forward states");
      }
      this_->lsolveB(t, cj, cjratio, b, weight, IDAADJ_mem->ia_yyTmp,
                     IDAADJ_mem->ia_ypTmp, xzB, xzdotB, rrB);
      return 0;
    } catch(int wrn) {
      /*    cerr << "warning: " << wrn << endl;*/
      return wrn;
    } catch(exception& e) {
      cerr << "lsolveB failed: " << e.what() << endl;
      return -1;
    }
  }

  void IdasInterface::lsetup(IDAMem IDA_mem, N_Vector xz, N_Vector xzdot, N_Vector resp,
                            N_Vector vtemp1, N_Vector vtemp2, N_Vector vtemp3) {
    log("IdasInterface::lsetup", "begin");

    // Current time
    double t = IDA_mem->ida_tn;

    // Multiple of df_dydot to be added to the matrix
    double cj = IDA_mem->ida_cj;

    // Call the preconditioner setup function (which sets up the linear solver)
    psetup(t, xz, xzdot, 0, cj, vtemp1, vtemp1, vtemp3);
    log("IdasInterface::lsetup", "end");
  }

  void IdasInterface::lsetupB(double t, double cj, N_Vector xz, N_Vector xzdot, N_Vector xzB,
                             N_Vector xzdotB, N_Vector resp,
                             N_Vector vtemp1, N_Vector vtemp2, N_Vector vtemp3) {
    log("IdasInterface::lsetupB", "begin");

    // Call the preconditioner setup function (which sets up the linear solver)
    psetupB(t, xz, xzdot, xzB, xzdotB, 0, cj, vtemp1, vtemp1, vtemp3);
    log("IdasInterface::lsetupB", "end");
  }


  void IdasInterface::lsolve(IDAMem IDA_mem, N_Vector b, N_Vector weight, N_Vector xz,
                            N_Vector xzdot, N_Vector rr) {
    log("IdasInterface::lsolve", "begin");
    // Current time
    double t = IDA_mem->ida_tn;

    // Multiple of df_dydot to be added to the matrix
    double cj = IDA_mem->ida_cj;

    // Accuracy
    double delta = 0.0;

    // Call the preconditioner solve function (which solves the linear system)
    psolve(t, xz, xzdot, rr, b, b, cj, delta, 0);

    // Scale the correction to account for change in cj
    if (cj_scaling_) {
      double cjratio = IDA_mem->ida_cjratio;
      if (cjratio != 1.0) N_VScale(2.0/(1.0 + cjratio), b, b);
    }
    log("IdasInterface::lsolve", "end");
  }

  void IdasInterface::lsolveB(double t, double cj, double cjratio, N_Vector b, N_Vector weight,
                             N_Vector xz, N_Vector xzdot, N_Vector xzB, N_Vector xzdotB,
                             N_Vector rr) {
    log("IdasInterface::lsolveB", "begin");

    // Accuracy
    double delta = 0.0;

    // Call the preconditioner solve function (which solves the linear system)
    psolveB(t, xz, xzdot, xzB, xzdotB, rr, b, b, cj, delta, 0);

    // Scale the correction to account for change in cj
    if (cj_scaling_) {
      if (cjratio != 1.0) N_VScale(2.0/(1.0 + cjratio), b, b);
    }

    log("IdasInterface::lsolveB", "end");
  }


  void IdasInterface::initDenseLinearSolver() {
    // Dense jacobian
    int flag = IDADense(mem_, nx_+nz_);
    if (flag != IDA_SUCCESS) idas_error("IDADense", flag);
    if (exact_jacobian_) {
      flag = IDADlsSetDenseJacFn(mem_, djac_wrapper);
      if (flag!=IDA_SUCCESS) idas_error("IDADlsSetDenseJacFn", flag);
    }
  }

  void IdasInterface::initBandedLinearSolver() {
    // Banded jacobian
    pair<int, int> bw = getBandwidth();
    int flag = IDABand(mem_, nx_+nz_, bw.first, bw.second);
    if (flag != IDA_SUCCESS) idas_error("IDABand", flag);

    // Banded Jacobian information
    if (exact_jacobian_) {
      flag = IDADlsSetBandJacFn(mem_, bjac_wrapper);
      if (flag != IDA_SUCCESS) idas_error("IDADlsSetBandJacFn", flag);
    }
  }

  void IdasInterface::initIterativeLinearSolver() {
    // Attach an iterative solver
    int flag;
    switch (itsol_f_) {
    case SD_GMRES:
      flag = IDASpgmr(mem_, max_krylov_);
      if (flag != IDA_SUCCESS) idas_error("IDASpgmr", flag);
      break;
    case SD_BCGSTAB:
      flag = IDASpbcg(mem_, max_krylov_);
      if (flag != IDA_SUCCESS) idas_error("IDASpbcg", flag);
      break;
    case SD_TFQMR:
      flag = IDASptfqmr(mem_, max_krylov_);
      if (flag != IDA_SUCCESS) idas_error("IDASptfqmr", flag);
      break;
    default: casadi_error("Uncaught switch");
    }

    // Attach functions for jacobian information
    if (exact_jacobian_) {
      // Form the Jacobian-times-vector function
      f_fwd_ = f_.derivative(1, 0);

      flag = IDASpilsSetJacTimesVecFn(mem_, jtimes_wrapper);
      if (flag != IDA_SUCCESS) idas_error("IDASpilsSetJacTimesVecFn", flag);
    }

    // Add a preconditioner
    if (use_preconditioner_) {
      // Make sure that a Jacobian has been provided
      if (jac_.isNull())
          throw CasadiException("IdasInterface::init(): No Jacobian has been provided.");

      // Make sure that a linear solver has been provided
      if (linsol_.isNull())
          throw CasadiException("IdasInterface::init(): "
                                "No user defined linear solver has been provided.");

      // Pass to IDA
      flag = IDASpilsSetPreconditioner(mem_, psetup_wrapper, psolve_wrapper);
      if (flag != IDA_SUCCESS) idas_error("IDASpilsSetPreconditioner", flag);
    }
  }

  void IdasInterface::initUserDefinedLinearSolver() {
    // Make sure that a Jacobian has been provided
    casadi_assert(!jac_.isNull());

    // Make sure that a linear solver has been provided
    casadi_assert(!linsol_.isNull());

    //  Set fields in the IDA memory
    IDAMem IDA_mem = IDAMem(mem_);
    IDA_mem->ida_lmem   = this;
    IDA_mem->ida_lsetup = lsetup_wrapper;
    IDA_mem->ida_lsolve = lsolve_wrapper;
    IDA_mem->ida_setupNonNull = TRUE;
  }

  void IdasInterface::initDenseLinearSolverB() {
    // Dense jacobian
    int flag = IDADenseB(mem_, whichB_, nrx_+nrz_);
    if (flag != IDA_SUCCESS) idas_error("IDADenseB", flag);
    if (exact_jacobianB_) {
      // Pass to IDA
      flag = IDADlsSetDenseJacFnB(mem_, whichB_, djacB_wrapper);
      if (flag!=IDA_SUCCESS) idas_error("IDADlsSetDenseJacFnB", flag);
    }
  }

  void IdasInterface::initBandedLinearSolverB() {
    pair<int, int> bw = getBandwidthB();
    int flag = IDABandB(mem_, whichB_, nrx_+nrz_, bw.first, bw.second);
    if (flag != IDA_SUCCESS) idas_error("IDABand", flag);
    if (exact_jacobianB_) {
      // Pass to IDA
      flag = IDADlsSetBandJacFnB(mem_, whichB_, bjacB_wrapper);
      if (flag!=IDA_SUCCESS) idas_error("IDADlsSetBandJacFnB", flag);
    }
  }

  // Workaround for bug in Sundials, to be removed when Debian/Homebrew
  // Sundials have been patched
#ifdef WITH_SYSTEM_SUNDIALS
  extern "C" {

#define yyTmp        (IDAADJ_mem->ia_yyTmp)
#define ypTmp        (IDAADJ_mem->ia_ypTmp)
#define noInterp     (IDAADJ_mem->ia_noInterp)

    static int IDAAspilsJacTimesVecPatched(realtype tt,
                                    N_Vector yyB, N_Vector ypB, N_Vector rrB,
                                    N_Vector vB, N_Vector JvB,
                                    realtype c_jB, void *ida_mem,
                                    N_Vector tmp1B, N_Vector tmp2B) {
      IDAMem IDA_mem;
      IDAadjMem IDAADJ_mem;
      IDASpilsMemB idaspilsB_mem;
      IDABMem IDAB_mem;
      int flag;

      IDA_mem = (IDAMem) ida_mem;
      IDAADJ_mem = IDA_mem->ida_adj_mem;
      IDAB_mem = IDAADJ_mem->ia_bckpbCrt;
      idaspilsB_mem = (IDASpilsMemB)IDAB_mem->ida_lmem;

      /* Get FORWARD solution from interpolation. */
      if (noInterp==FALSE) {
        flag = IDAADJ_mem->ia_getY(IDA_mem, tt, yyTmp, ypTmp, NULL, NULL);
        if (flag != IDA_SUCCESS) {
          IDAProcessError(IDA_mem, -1, "IDASSPILS", "IDAAspilsJacTimesVec", MSGS_BAD_T);
          return(-1);
        }
      }
      /* Call user's adjoint psolveB routine */
      flag = idaspilsB_mem->s_jtimesB(tt, yyTmp, ypTmp,
                                      yyB, ypB, rrB,
                                      vB, JvB,
                                      c_jB, IDAB_mem->ida_user_data,
                                      tmp1B, tmp2B);
      return(flag);
    }

    int IDASpilsSetJacTimesVecFnBPatched(void *ida_mem, int which, IDASpilsJacTimesVecFnB jtvB) {
      IDAadjMem IDAADJ_mem;
      IDAMem IDA_mem;
      IDABMem IDAB_mem;
      IDASpilsMemB idaspilsB_mem;
      void *ida_memB;
      int flag;

      /* Check if ida_mem is allright. */
      if (ida_mem == NULL) {
        IDAProcessError(NULL, IDASPILS_MEM_NULL, "IDASSPILS",
                        "IDASpilsSetJacTimesVecFnB", MSGS_IDAMEM_NULL);
        return(IDASPILS_MEM_NULL);
      }
      IDA_mem = (IDAMem) ida_mem;

      /* Is ASA initialized? */
      if (IDA_mem->ida_adjMallocDone == FALSE) {
        IDAProcessError(IDA_mem, IDASPILS_NO_ADJ, "IDASSPILS",
                        "IDASpilsSetJacTimesVecFnB",  MSGS_NO_ADJ);
        return(IDASPILS_NO_ADJ);
      }
      IDAADJ_mem = IDA_mem->ida_adj_mem;

      /* Check the value of which */
      if ( which >= IDAADJ_mem->ia_nbckpbs ) {
        IDAProcessError(IDA_mem, IDASPILS_ILL_INPUT, "IDASSPILS",
                        "IDASpilsSetJacTimesVecFnB", MSGS_BAD_WHICH);
        return(IDASPILS_ILL_INPUT);
      }

      /* Find the IDABMem entry in the linked list corresponding to 'which'. */
      IDAB_mem = IDAADJ_mem->IDAB_mem;
      while (IDAB_mem != NULL) {
        if ( which == IDAB_mem->ida_index ) break;
        /* advance */
        IDAB_mem = IDAB_mem->ida_next;
      }
      /* ida_mem corresponding to 'which' problem. */
      ida_memB = reinterpret_cast<void *>(IDAB_mem->IDA_mem);

      if (IDAB_mem->ida_lmem == NULL) {
        IDAProcessError(IDA_mem, IDASPILS_LMEMB_NULL, "IDASSPILS",
                        "IDASpilsSetJacTimesVecFnB", MSGS_LMEMB_NULL);
        return(IDASPILS_ILL_INPUT);
      }

      idaspilsB_mem = (IDASpilsMemB) IDAB_mem->ida_lmem;

      /* Call the corresponding Set* function for the backward problem. */

      idaspilsB_mem->s_jtimesB   = jtvB;

      if (jtvB != NULL) {
        flag = IDASpilsSetJacTimesVecFn(ida_memB, IDAAspilsJacTimesVecPatched);
      } else {
        flag = IDASpilsSetJacTimesVecFn(ida_memB, NULL);
      }
      return(flag);
    }
  }
#endif // WITH_SYSTEM_SUNDIALS

  void IdasInterface::initIterativeLinearSolverB() {
    int flag;
    switch (itsol_g_) {
    case SD_GMRES:
      flag = IDASpgmrB(mem_, whichB_, max_krylovB_);
      if (flag != IDA_SUCCESS) idas_error("IDASpgmrB", flag);
      break;
    case SD_BCGSTAB:
      flag = IDASpbcgB(mem_, whichB_, max_krylovB_);
      if (flag != IDA_SUCCESS) idas_error("IDASpbcgB", flag);
      break;
    case SD_TFQMR:
      flag = IDASptfqmrB(mem_, whichB_, max_krylovB_);
      if (flag != IDA_SUCCESS) idas_error("IDASptfqmrB", flag);
      break;
    default: casadi_error("Uncaught switch");
    }

    // Attach functions for jacobian information
    if (exact_jacobianB_) {
      // Form the Jacobian-times-vector function
      g_fwd_ = g_.derivative(1, 0);

#ifdef WITH_SYSTEM_SUNDIALS
      flag = IDASpilsSetJacTimesVecFnBPatched(mem_, whichB_, jtimesB_wrapper);
#else
      flag = IDASpilsSetJacTimesVecFnB(mem_, whichB_, jtimesB_wrapper);
#endif
      if (flag != IDA_SUCCESS) idas_error("IDASpilsSetJacTimesVecFnB", flag);
    }

    // Add a preconditioner
    if (use_preconditionerB_) {
      // Make sure that a Jacobian has been provided
      if (jacB_.isNull())
          throw CasadiException("IdasInterface::init(): No backwards Jacobian has been provided.");

      // Make sure that a linear solver has been provided
      if (linsolB_.isNull())
          throw CasadiException("IdasInterface::init(): No backwards user "
                                "defined linear solver has been provided.");

      // Pass to IDA
      flag = IDASpilsSetPreconditionerB(mem_, whichB_, psetupB_wrapper, psolveB_wrapper);
      if (flag != IDA_SUCCESS) idas_error("IDASpilsSetPreconditionerB", flag);
    }

  }

  void IdasInterface::initUserDefinedLinearSolverB() {
    // Make sure that a Jacobian has been provided
    casadi_assert(!jacB_.isNull());

    // Make sure that a linear solver has been provided
    casadi_assert(!linsolB_.isNull());

    //  Set fields in the IDA memory
    IDAMem IDA_mem = IDAMem(mem_);
    IDAadjMem IDAADJ_mem = IDA_mem->ida_adj_mem;
    IDABMem IDAB_mem = IDAADJ_mem->IDAB_mem;
    IDAB_mem->ida_lmem   = this;

    IDAB_mem->IDA_mem->ida_lmem = this;
    IDAB_mem->IDA_mem->ida_lsetup = lsetupB_wrapper;
    IDAB_mem->IDA_mem->ida_lsolve = lsolveB_wrapper;
    IDAB_mem->IDA_mem->ida_setupNonNull = TRUE;
  }

  template<typename FunctionType>
  FunctionType IdasInterface::getJacGen() {
    FunctionType f = shared_cast<FunctionType>(f_);
    casadi_assert(!f.isNull());

    // Get the Jacobian in the Newton iteration
    typename FunctionType::MatType cj = FunctionType::MatType::sym("cj");
    typename FunctionType::MatType jac = f.jac(DAE_X, DAE_ODE) - cj*FunctionType::MatType::eye(nx_);
    if (nz_>0) {
      jac = horzcat(vertcat(jac, f.jac(DAE_X, DAE_ALG)),
                    vertcat(f.jac(DAE_Z, DAE_ODE), f.jac(DAE_Z, DAE_ALG)));
    }

    // Jacobian function
    std::vector<typename FunctionType::MatType> jac_in = f.inputExpr();
    jac_in.push_back(cj);

    // Return generated function
    return FunctionType(jac_in, jac);
  }

  template<typename FunctionType>
  FunctionType IdasInterface::getJacGenB() {
    FunctionType g = shared_cast<FunctionType>(g_);
    casadi_assert(!g.isNull());

    // Get the Jacobian in the Newton iteration
    typename FunctionType::MatType cj = FunctionType::MatType::sym("cj");
    typename FunctionType::MatType jac =
        g.jac(RDAE_RX, RDAE_ODE) + cj*FunctionType::MatType::eye(nrx_);
    if (nrz_>0) {
      jac = horzcat(vertcat(jac, g.jac(RDAE_RX, RDAE_ALG)),
                    vertcat(g.jac(RDAE_RZ, RDAE_ODE), g.jac(RDAE_RZ, RDAE_ALG)));
    }

    // Jacobian function
    std::vector<typename FunctionType::MatType> jac_in = g.inputExpr();
    jac_in.push_back(cj);

    // return generated function
    return FunctionType(jac_in, jac);
  }

  Function IdasInterface::getJacB() {
    if (is_a<SXFunction>(g_)) {
      return getJacGenB<SXFunction>();
    } else if (is_a<MXFunction>(g_)) {
      return getJacGenB<MXFunction>();
    } else {
      throw CasadiException("IdasInterface::getJacB(): Not an SXFunction or MXFunction");
    }
  }

  Function IdasInterface::getJac() {
    if (is_a<SXFunction>(f_)) {
      return getJacGen<SXFunction>();
    } else if (is_a<MXFunction>(f_)) {
      return getJacGen<MXFunction>();
    } else {
      throw CasadiException("IdasInterface::getJac(): Not an SXFunction or MXFunction");
    }
  }


} // namespace casadi

