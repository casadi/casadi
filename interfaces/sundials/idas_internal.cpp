/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010 by Joel Andersson, Moritz Diehl, K.U.Leuven. All rights reserved.
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

#include "idas_internal.hpp"
#include "casadi/stl_vector_tools.hpp"
#include "casadi/fx/linear_solver_internal.hpp"
#include "casadi/fx/mx_function.hpp"
#include "casadi/sx/sx_tools.hpp"
#include "casadi/mx/mx_tools.hpp"

using namespace std;
namespace CasADi{
namespace Sundials{

IdasInternal* IdasInternal::clone() const{
  // Return a deep copy
  IdasInternal* node = new IdasInternal(f_,g_);
  node->setOption(dictionary());
  node->jac_ = jac_;
  node->linsol_ = linsol_;
  return node;
}

void IdasInternal::deepCopyMembers(std::map<SharedObjectNode*,SharedObject>& already_copied){
  SundialsInternal::deepCopyMembers(already_copied);
}

IdasInternal::IdasInternal(const FX& f, const FX& g) : SundialsInternal(f,g){
  addOption("suppress_algebraic",          OT_BOOLEAN, false, "supress algebraic variables in the error testing");
  addOption("calc_ic",                     OT_BOOLEAN, true,  "use IDACalcIC to get consistent initial conditions. This only works for semi-explicit index-one systems. Else, you must provide consistent initial conditions yourself.");
  addOption("calc_icB",                    OT_BOOLEAN, false, "use IDACalcIC to get consistent initial conditions. This only works for semi-explicit index-one systems. Else, you must provide consistent initial conditions yourself.");
  addOption("abstolv",                     OT_REALVECTOR);
  addOption("fsens_abstolv",               OT_REALVECTOR); 
  addOption("max_step_size",               OT_REAL, 0, "maximim step size");
  addOption("first_time",                  OT_REAL, GenericType(), "first requested time as a fraction of the time interval");
  addOption("cj_scaling",                  OT_BOOLEAN, false, "IDAS scaling on cj for the user-defined linear solver module");
  addOption("extra_fsens_calc_ic",         OT_BOOLEAN, false, "Call calc ic an extra time, with fsens=0");
  addOption("disable_internal_warnings",   OT_BOOLEAN,false, "Disable IDAS internal warning messages");
  addOption("monitor",                     OT_STRINGVECTOR, GenericType(),  "", "correctInitialConditions|res|resS", true);
  addOption("init_xdot",                   OT_REALVECTOR, GenericType(), "Initial values for the state derivatives");
  addOption("init_z",                      OT_REALVECTOR, GenericType(), "Initial values for the algebraic states");
  
  mem_ = 0;
  
  xz_  = 0; 
  xzdot_ = 0, 
  q_ = 0;

  rxz_ = 0;
  rxzdot_ = 0;
  rq_ = 0;

  is_init = false;
  isInitAdj_ = false;
  isInitTaping_ = false;
  disable_internal_warnings_ = false;
}

IdasInternal::~IdasInternal(){ 
  freeIDAS();
}

void IdasInternal::freeIDAS(){
  if(mem_) IDAFree(&mem_);

  // Forward integration
  if(xz_) N_VDestroy_Serial(xz_);
  if(xzdot_) N_VDestroy_Serial(xzdot_);
  if(q_) N_VDestroy_Serial(q_);
  
  // Backward integration
  if(rxz_) N_VDestroy_Serial(rxz_);
  if(rxzdot_) N_VDestroy_Serial(rxzdot_);
  if(rq_) N_VDestroy_Serial(rq_);
  
    // Forward problem
  for(vector<N_Vector>::iterator it=xzF_.begin(); it != xzF_.end(); ++it)   if(*it) N_VDestroy_Serial(*it);
  for(vector<N_Vector>::iterator it=xzdotF_.begin(); it != xzdotF_.end(); ++it)   if(*it) N_VDestroy_Serial(*it);
  for(vector<N_Vector>::iterator it=qF_.begin(); it != qF_.end(); ++it)   if(*it) N_VDestroy_Serial(*it);
}

void IdasInternal::updateNumSens(bool recursive){
  // Not supported re-initalization needed
  init();
}

void IdasInternal::init(){
  log("IdasInternal::init","begin");

  // Call the base class init
  SundialsInternal::init();

  // Free memory if already initialized
  if(is_init) freeIDAS();
  
  if(hasSetOption("linear_solver_creator")){
    // Make sure that a Jacobian has been provided
    if(jac_.isNull()) getJacobian();
    if(!jac_.isInit()) jac_.init();
    
    // Create a linear solver
    linearSolverCreator creator = getOption("linear_solver_creator");
    linsol_ = creator(jac_.output().sparsity());
    linsol_.setSparsity(jac_.output().sparsity());
    linsol_.init();
  } else {
    if(!jac_.isNull()){
      jac_.init();
      if(!linsol_.isNull()){
        linsol_.setSparsity(jac_.output().sparsity());
        linsol_.init();
      }
    }
  }
  
  // Get initial conditions for the state derivatives
  if(hasSetOption("init_xdot")){
    init_xdot_ = getOption("init_xdot").toDoubleVector();
    casadi_assert_message(init_xdot_.size()==nx_,"Option \"init_xdot\" has incorrect length.");
  } else {
    init_xdot_.resize(nx_);
    fill(init_xdot_.begin(),init_xdot_.end(),0);
  }
  
  // Get initial conditions for the algebraic state
  if(hasSetOption("init_z")){
    init_z_ = getOption("init_z").toDoubleVector();
    casadi_assert_message(init_z_.size()==nz_,"Option \"init_z\" has incorrect length.");
  } else {
    init_z_.resize(nz_);
    fill(init_z_.begin(),init_z_.end(),0);
  }

  // Get the number of forward and adjoint directions
  nfdir_f_ = f_.getOption("number_of_fwd_dir");
  nadir_f_ = f_.getOption("number_of_adj_dir");

  cj_scaling_ = getOption("cj_scaling");
  
  // Sundials return flag
  int flag;

  // Create IDAS memory block
  mem_ = IDACreate();
  if(mem_==0) throw CasadiException("IDACreate(): Creation failed");

  // Allocate n-vectors for ivp
  xz_ = N_VNew_Serial(nx_+nz_);
  xzdot_ = N_VNew_Serial(nx_+nz_);

  // Initialize Idas
  double t0 = 0;
  N_VConst(0.0, xz_);
  N_VConst(0.0, xzdot_);
  IDAInit(mem_, res_wrapper, t0, xz_, xzdot_);
  log("IdasInternal::init","IDA initialized");

  // Disable internal warning messages?
  disable_internal_warnings_ = getOption("disable_internal_warnings");
  
  // Set error handler function
  flag = IDASetErrHandlerFn(mem_, ehfun_wrapper, this);
  casadi_assert_message(flag == IDA_SUCCESS,"IDASetErrHandlerFn");

  // Include algebraic variables in error testing
  flag = IDASetSuppressAlg(mem_, getOption("suppress_algebraic").toInt());
  casadi_assert_message(flag == IDA_SUCCESS,"IDASetSuppressAlg");

  // Maxinum order for the multistep method
  flag = IDASetMaxOrd(mem_, getOption("max_multistep_order").toInt());
  casadi_assert_message(flag == IDA_SUCCESS,"IDASetMaxOrd");

  // Set user data
  flag = IDASetUserData(mem_,this);
  casadi_assert_message(flag == IDA_SUCCESS,"IDASetUserData");

  // Set maximum step size
  flag = IDASetMaxStep(mem_, getOption("max_step_size").toDouble());
  casadi_assert_message(flag == IDA_SUCCESS,"IDASetMaxStep");
  
  if(hasSetOption("abstolv")){
    // Vector absolute tolerances
    vector<double> abstolv = getOption("abstolv").toDoubleVector();
    N_Vector nv_abstol = N_VMake_Serial(abstolv.size(),getPtr(abstolv));
    flag = IDASVtolerances(mem_, reltol_, nv_abstol);
    casadi_assert_message(flag == IDA_SUCCESS,"IDASVtolerances");
    N_VDestroy_Serial(nv_abstol);
  } else {
    // Scalar absolute tolerances
    flag = IDASStolerances(mem_, reltol_, abstol_); 
    casadi_assert_message(flag == IDA_SUCCESS,"IDASStolerances");
  }
  
  // Maximum number of steps
  IDASetMaxNumSteps(mem_, getOption("max_num_steps").toInt());
  if(flag != IDA_SUCCESS) idas_error("IDASetMaxNumSteps",flag);

  // Set algebraic components
  N_Vector id = N_VNew_Serial(nx_+nz_);
  fill_n(NV_DATA_S(id),nx_,1);
  fill_n(NV_DATA_S(id)+nx_,nz_,0);
  
  // Pass this information to IDAS
  flag = IDASetId(mem_, id);
  if(flag != IDA_SUCCESS) idas_error("IDASetId",flag);

  // Delete the allocated memory
  N_VDestroy_Serial(id);
  
  // attach a linear solver
  switch(linsol_f_){
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
  }
  
  // Quadrature equations
  if(nq_>0){

    // Allocate n-vectors for quadratures
    q_ = N_VMake_Serial(nq_,&output(INTEGRATOR_QF).front());

    // Initialize quadratures in IDAS
    N_VConst(0.0, q_);
    flag = IDAQuadInit(mem_, rhsQ_wrapper, q_);
    if(flag != IDA_SUCCESS) idas_error("IDAQuadInit",flag);
    
    // Should the quadrature errors be used for step size control?
    if(getOption("quad_err_con").toInt()){
      flag = IDASetQuadErrCon(mem_, true);
      casadi_assert_message(flag == IDA_SUCCESS, "IDASetQuadErrCon");
      
      // Quadrature error tolerances
      flag = IDAQuadSStolerances(mem_, reltol_, abstol_); // TODO: vector absolute tolerances
      if(flag != IDA_SUCCESS) idas_error("IDAQuadSStolerances",flag);
    }
  }
  
  log("IdasInternal::init","attached linear solver");
    
    // Forward sensitivity problem
    if(nfdir_>0){
      // Allocate n-vectors
      xzF_.resize(nfdir_,0);
      xzdotF_.resize(nfdir_,0);
      qF_.resize(nfdir_,0);
      for(int i=0; i<nfdir_; ++i){
        xzF_[i] = N_VNew_Serial(nx_+nz_);
        xzdotF_[i] = N_VNew_Serial(nx_+nz_);
	if(nq_>0){
	  qF_[i] = N_VMake_Serial(nq_,&fwdSens(INTEGRATOR_QF,i).front());
	}
      }
      
    // Get the sensitivity method
    if(getOption("sensitivity_method")== "simultaneous") ism_ = IDA_SIMULTANEOUS;
    else if(getOption("sensitivity_method")=="staggered") ism_ = IDA_STAGGERED;
    else throw CasadiException("IDAS: Unknown sensitivity method");

    // Copy the forward seeds
    for(int i=0; i<nfdir_; ++i){
      const Matrix<double>& x = fwdSeed(INTEGRATOR_X0,i);
      copy(x.begin(), x.begin()+nx_, NV_DATA_S(xzF_[i]));
      N_VConst(0.0,xzdotF_[i]);
    }

    // Initialize forward sensitivities
    if(finite_difference_fsens_){
      
      // Use finite differences to calculate the residual in the forward sensitivity equations
      flag = IDASensInit(mem_,nfdir_,ism_,0,getPtr(xzF_),getPtr(xzdotF_));
      if(flag != IDA_SUCCESS) idas_error("IDASensInit",flag);

      // Scaling factors
      double* pbar = 0;
      vector<double> pbar_vec;
      if(hasSetOption("fsens_scaling_factors")){
        pbar_vec = getOption("fsens_scaling_factors").toDoubleVector();
        pbar = getPtr(pbar_vec);
      }
      
      // Which parameters should be used to estimate the sensitivity equations
      int * plist = 0;
      vector<int> plist_vec;
      if(hasSetOption("fsens_sensitivity_parameters")){
        plist_vec = getOption("fsens_sensitiviy_parameters").toIntVector();
        plist = getPtr(plist_vec);
      }

      // Specify parameters
      flag = IDASetSensParams(mem_,&input(INTEGRATOR_P).front(),pbar,plist);
      if(flag != IDA_SUCCESS) idas_error("IDASetSensParams",flag);

      //  IDASetSensDQMethod
      // ?

    } else {
      // Use AD to calculate the residual in the forward sensitivity equations
      flag = IDASensInit(mem_,nfdir_,ism_,resS_wrapper,&xzF_.front(),&xzdotF_.front());
      if(flag != IDA_SUCCESS) idas_error("IDASensInit",flag);
    }

    // IDAS bugfix
    IDAMem IDA_mem = IDAMem(mem_);
    int max_multistep_order = getOption("max_multistep_order");
    for(int i=0; i<=max_multistep_order; ++i){
      for(int iS=0; iS<nfdir_; ++iS){
        N_VConst(0.0, IDA_mem->ida_phiS[i][iS]);
      }
    }

    // Vector absolute tolerances
    vector<double> fsens_abstol(nfdir_,fsens_abstol_);
    
    // Set tolerances
    if(hasSetOption("fsens_abstolv")){
      // quick hack
      vector<double> fsens_abstolv = getOption("fsens_abstolv");
      N_Vector nv_abstol = N_VMake_Serial(fsens_abstolv.size(),getPtr(fsens_abstolv));
      vector<N_Vector> nv_abstol_all(nfdir_,nv_abstol);
      flag = IDASensSVtolerances(mem_,fsens_reltol_,getPtr(nv_abstol_all));
      if(flag != IDA_SUCCESS) idas_error("IDASensSVtolerances",flag);
      N_VDestroy_Serial(nv_abstol);
    } else {
      flag = IDASensSStolerances(mem_,fsens_reltol_,getPtr(fsens_abstol));
      if(flag != IDA_SUCCESS) idas_error("IDASensSStolerances",flag);
    }

    // Set optional inputs
    bool errconS = getOption("fsens_err_con");
    flag = IDASetSensErrCon(mem_, errconS);
    if(flag != IDA_SUCCESS) idas_error("IDASetSensErrCon",flag);

    // Quadrature equations
    if(nq_>0){
      for(vector<N_Vector>::iterator it=qF_.begin(); it!=qF_.end(); ++it) N_VConst(0.0,*it);
      flag = IDAQuadSensInit(mem_, rhsQS_wrapper, getPtr(qF_));
      if(flag != IDA_SUCCESS) idas_error("IDAQuadSensInit",flag);
      
      // Set tolerances
      flag = IDAQuadSensSStolerances(mem_,fsens_reltol_,getPtr(fsens_abstol));
      if(flag != IDA_SUCCESS) idas_error("IDAQuadSensSStolerances",flag);
    }
    
    log("IdasInternal::init","initialized forward sensitivities");
  } // enable fsens

  // Adjoint sensitivity problem
  if(!g_.isNull()){
    
    // Allocate n-vectors
    rxz_ = N_VNew_Serial(nrx_+nrz_);
    rxzdot_ = N_VNew_Serial(nrx_+nrz_);
    N_VConst(0.0, rxz_);
    N_VConst(0.0, rxzdot_);

    // Allocate n-vectors for quadratures
    rq_ = N_VMake_Serial(nrq_,output(INTEGRATOR_RQF).ptr());
  }
  log("IdasInternal::init","initialized adjoint sensitivities");
 
 is_init = true;
 isInitTaping_ = false;
 isInitAdj_ = false;
 log("IdasInternal::init","end");
}

void IdasInternal::initTaping(){
  casadi_assert(!isInitTaping_);
  int flag;
  
  // Get the number of steos per checkpoint
  int Nd = getOption("steps_per_checkpoint");

  // Get the interpolation type
  int interpType;
  if(getOption("interpolation_type")=="hermite")
    interpType = IDA_HERMITE;
  else if(getOption("interpolation_type")=="polynomial")
    interpType = IDA_POLYNOMIAL;
  else throw CasadiException("\"interpolation_type\" must be \"hermite\" or \"polynomial\"");
      
  // Initialize adjoint sensitivities
  flag = IDAAdjInit(mem_, Nd, interpType);
  if(flag != IDA_SUCCESS) idas_error("IDAAdjInit",flag);
  
  isInitTaping_ = true;
}

void IdasInternal::initAdj(){
  casadi_assert(!isInitAdj_);
  int flag;
  
  // Create backward problem
  flag = IDACreateB(mem_, &whichB_);
  if(flag != IDA_SUCCESS) idas_error("IDACreateB",flag);

  // Initialize the backward problem
  double tB0 = tf_;
  flag = IDAInitB(mem_, whichB_, resB_wrapper, tB0, rxz_, rxzdot_);
  if(flag != IDA_SUCCESS) idas_error("IDAInitB",flag);

  // Set tolerances
  flag = IDASStolerancesB(mem_, whichB_, asens_reltol_, asens_abstol_);
  if(flag!=IDA_SUCCESS) idas_error("IDASStolerancesB",flag);

  // User data
  flag = IDASetUserDataB(mem_, whichB_, this);
  if(flag != IDA_SUCCESS) idas_error("IDASetUserDataB",flag);

  // Maximum number of steps
  IDASetMaxNumStepsB(mem_, whichB_, getOption("max_num_steps").toInt());
  if(flag != IDA_SUCCESS) idas_error("IDASetMaxNumStepsB",flag);

  // Set algebraic components
  N_Vector id = N_VNew_Serial(nrx_+nrz_);
  fill_n(NV_DATA_S(id),nrx_,1);
  fill_n(NV_DATA_S(id)+nrx_,nrz_,0);
  
  // Pass this information to IDAS
  flag = IDASetIdB(mem_, whichB_, id);
  if(flag != IDA_SUCCESS) idas_error("IDASetIdB",flag);

  // Delete the allocated memory
  N_VDestroy_Serial(id);
    
  // attach linear solver
  switch(linsol_g_){
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
  }
  
  // Quadratures for the adjoint problem
  N_VConst(0.0,rq_);
  flag = IDAQuadInitB(mem_,whichB_,rhsQB_wrapper,rq_);
  if(flag!=IDA_SUCCESS) idas_error("IDAQuadInitB",flag);

  // Quadrature error control
  if(getOption("quad_err_con").toInt()){
    flag = IDASetQuadErrConB(mem_, whichB_,true);
    if(flag != IDA_SUCCESS) idas_error("IDASetQuadErrConB",flag);
    
    flag = IDAQuadSStolerancesB(mem_, whichB_, asens_reltol_, asens_abstol_);
    if(flag != IDA_SUCCESS) idas_error("IDAQuadSStolerancesB",flag);
  }
  
  // Mark initialized
  isInitAdj_ = true;
}


void IdasInternal::res(double t, const double* xz, const double* xzdot, double* r){
  log("IdasInternal::res","begin");
  
  // Get time
  time1 = clock();
  
  // Pass input
  f_.setInput(&t,DAE_T);
  f_.setInput(xz,DAE_X);
  f_.setInput(xz+nx_,DAE_Z);
  f_.setInput(xzdot,DAE_XDOT);
  f_.setInput(input(INTEGRATOR_P),DAE_P);

  // Evaluate
  f_.evaluate();
  
  // Get results
  f_.getOutput(r,DAE_ODE);
  f_.getOutput(r+nx_,DAE_ALG);
  
  if(monitored("res")){
    cout << "DAE_T    = " << t << endl;
    cout << "DAE_X    = " << f_.input(DAE_X) << endl;
    cout << "DAE_Z    = " << f_.input(DAE_Z) << endl;
    cout << "DAE_XDOT = " << f_.input(DAE_XDOT) << endl;
    cout << "DAE_P    = " << f_.input(DAE_P) << endl;
    cout << "residual = " << f_.output() << endl;
  }

  // Check the result for consistency
  for(int i=0; i<nx_+nz_; ++i){
    if(isnan(r[i]) || isinf(r[i])){
      if(verbose_){
        stringstream ss;
        ss << "Warning: The " << i << "-th component of the DAE residual is " << r[i] << " at time t=" << t << ".";
        log("IdasInternal::res",ss.str());
        if(monitored("res")){
          cout << "DAE_T    = " << t << endl;
          cout << "DAE_X    = " << f_.input(DAE_X) << endl;
          cout << "DAE_Z    = " << f_.input(DAE_Z) << endl;
          cout << "DAE_XDOT = " << f_.input(DAE_XDOT) << endl;
          cout << "DAE_P    = " << f_.input(DAE_P) << endl;
          cout << "ODE residual = " << f_.output(DAE_ODE) << endl;
          cout << "ALG residual = " << f_.output(DAE_ALG) << endl;
        }
      }
      throw 1;
    }
  }
  
  
  time2 = clock();
  t_res += double(time2-time1)/CLOCKS_PER_SEC;
  log("IdasInternal::res","end");
}

int IdasInternal::res_wrapper(double t, N_Vector xz, N_Vector xzdot, N_Vector rr, void *user_data){
 try{
    IdasInternal *this_ = (IdasInternal*)user_data;
    this_->res(t,NV_DATA_S(xz),NV_DATA_S(xzdot),NV_DATA_S(rr));
    return 0;
  } catch(int flag){ // recoverable error
    return flag;
  } catch(exception& e){ // non-recoverable error
    cerr << "res failed: " << e.what() << endl;
    return -1;
  }
}

void IdasInternal::ehfun_wrapper(int error_code, const char *module, const char *function, char *msg, void *eh_data){
 try{
    IdasInternal *this_ = (IdasInternal*)eh_data;
    this_->ehfun(error_code,module,function,msg);        
  } catch(exception& e){
    cerr << "ehfun failed: " << e.what() << endl;
  }
}
  
void IdasInternal::ehfun(int error_code, const char *module, const char *function, char *msg){
  cerr << msg << endl;
}

void IdasInternal::jtimes(double t, const double *xz, const double *xzdot, const double *rr, const double *v, double *Jv, double cj, double *tmp1, double *tmp2){
  log("IdasInternal::jtimes","begin");
  // Get time
  time1 = clock();
  
   // Pass input
   f_.setInput(&t,DAE_T);
   f_.setInput(xz,DAE_X);
   f_.setInput(xz+nx_,DAE_Z);
   f_.setInput(xzdot,DAE_XDOT);
   f_.setInput(input(INTEGRATOR_P),DAE_P);
     
   // Pass seeds of the state vectors
   f_.setFwdSeed(v,DAE_X);
   f_.setFwdSeed(v+nx_,DAE_Z);
   
   // Pass seeds of the state derivative
   for(int i=0; i<nx_; ++i) tmp1[i] = cj*v[i];
   f_.setFwdSeed(tmp1,DAE_XDOT);
   
   // Evaluate the AD forward algorithm
   f_.evaluate(1,0);
   
   // Get the output seeds
   f_.getFwdSens(Jv);

  // Log time duration
  time2 = clock();
  t_jac += double(time2-time1)/CLOCKS_PER_SEC;
  log("IdasInternal::jtimes","end");
}

int IdasInternal::jtimes_wrapper(double t, N_Vector xz, N_Vector xzdot, N_Vector rr, N_Vector v, N_Vector Jv, double cj, void *user_data, N_Vector tmp1, N_Vector tmp2){
 try{
    IdasInternal *this_ = (IdasInternal*)user_data;
    this_->jtimes(t,NV_DATA_S(xz),NV_DATA_S(xzdot),NV_DATA_S(rr),NV_DATA_S(v),NV_DATA_S(Jv),cj,NV_DATA_S(tmp1),NV_DATA_S(tmp2));
    return 0;
  } catch(exception& e){
    cerr << "jtimes failed: " << e.what() << endl;
    return 1;
  }
}

void IdasInternal::resS(int Ns, double t, const double* xz, const double* xzdot, const double *resval, N_Vector *xzF, N_Vector* xzdotF, N_Vector *rrF, double *tmp1, double *tmp2, double *tmp3){
  log("IdasInternal::resS","begin");
  casadi_assert(Ns==nfdir_);

  // Record the current cpu time
  time1 = clock();
  
  // Pass input
  f_.setInput(&t,DAE_T);
  f_.setInput(xz,DAE_X);
  f_.setInput(xz+nx_,DAE_Z);
  f_.setInput(xzdot,DAE_XDOT);
  f_.setInput(input(INTEGRATOR_P),DAE_P);
  
  // Calculate the forward sensitivities, nfdir_f_ directions at a time
  for(int offset=0; offset<nfdir_; offset += nfdir_f_){
    // Number of directions in this batch
    int nfdir_batch = std::min(nfdir_-offset, nfdir_f_);
    for(int dir=0; dir<nfdir_batch; ++dir){
      // Pass forward seeds
      f_.fwdSeed(DAE_T,dir).setZero();
      f_.setFwdSeed(NV_DATA_S(xzF[offset+dir]),DAE_X,dir);
      f_.setFwdSeed(NV_DATA_S(xzF[offset+dir])+nx_,DAE_Z,dir);
      f_.setFwdSeed(NV_DATA_S(xzdotF[offset+dir]),DAE_XDOT,dir);
      f_.setFwdSeed(fwdSeed(INTEGRATOR_P,offset+dir),DAE_P,dir);
    }
    
    // Evaluate the AD forward algorithm
    f_.evaluate(nfdir_batch,0);
    
    // Get the output seeds
    for(int dir=0; dir<nfdir_batch; ++dir){
      f_.getFwdSens(NV_DATA_S(rrF[offset+dir]),DAE_ODE,dir);
      f_.getFwdSens(NV_DATA_S(rrF[offset+dir])+nx_,DAE_ALG,dir);
    }
  }
  
  // Record timings
  time2 = clock();
  t_fres += double(time2-time1)/CLOCKS_PER_SEC;
  log("IdasInternal::resS","end");
}

int IdasInternal::resS_wrapper(int Ns, double t, N_Vector xz, N_Vector xzdot, N_Vector resval, N_Vector *xzF, N_Vector *xzdotF, N_Vector *rrF, void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3){
 try{
    IdasInternal *this_ = (IdasInternal*)user_data;
    this_->resS(Ns,t,NV_DATA_S(xz),NV_DATA_S(xzdot),NV_DATA_S(resval),xzF,xzdotF,rrF,NV_DATA_S(tmp1),NV_DATA_S(tmp2),NV_DATA_S(tmp3));
    return 0;
  } catch(exception& e){
    cerr << "resS failed: " << e.what() << endl;
    return 1;
  }
}

void IdasInternal::reset(int nfdir){
  log("IdasInternal::reset","begin");
  
  if(!g_.isNull() && !isInitTaping_)
    initTaping();
  
  // If we have forward sensitivities, rest one extra time without forward sensitivities to get a consistent initial guess
  if(nfdir>0 && getOption("extra_fsens_calc_ic").toInt())
    reset(0);

  // Reset timers
  t_res = t_fres = t_jac = t_lsolve = t_lsetup_jac = t_lsetup_fac = 0;
    
  fsens_order_ = nfdir>0;
  asens_order_ = !g_.isNull();
  
  // Get the time horizon
  t_ = t0_;
  
  // Return flag
  int flag;
  
  // Copy to N_Vectors
  const Matrix<double>& x0 = input(INTEGRATOR_X0);
  copy(x0.begin(), x0.begin()+nx_, NV_DATA_S(xz_));
  copy(init_z_.begin(), init_z_.end(), NV_DATA_S(xz_)+nx_);
  copy(init_xdot_.begin(), init_xdot_.end(), NV_DATA_S(xzdot_));
  
  // Re-initialize
  flag = IDAReInit(mem_, t0_, xz_, xzdot_);
  if(flag != IDA_SUCCESS) idas_error("IDAReInit",flag);
  log("IdasInternal::reset","re-initialized IVP solution");


  // Re-initialize quadratures
  if(nq_>0){
    N_VConst(0.0,q_);
    flag = IDAQuadReInit(mem_, q_);
    if(flag != IDA_SUCCESS) idas_error("IDAQuadReInit",flag);
    log("IdasInternal::reset","re-initialized quadratures");
  }
  
  if(fsens_order_>0){
    // Get the forward seeds
    for(int i=0; i<nfdir_; ++i){
      const Matrix<double>& x0_seed = fwdSeed(INTEGRATOR_X0,i);
      copy(x0_seed.begin(), x0_seed.begin()+nx_, NV_DATA_S(xzF_[i]));
      N_VConst(0.0,xzdotF_[i]);
    }
    
    // Re-initialize sensitivities
    flag = IDASensReInit(mem_,ism_,getPtr(xzF_),getPtr(xzdotF_));
    if(flag != IDA_SUCCESS) idas_error("IDASensReInit",flag);
    log("IdasInternal::reset","re-initialized forward sensitivity solution");

    if(nq_>0){
      for(vector<N_Vector>::iterator it=qF_.begin(); it!=qF_.end(); ++it) N_VConst(0.0,*it);
      flag = IDAQuadSensReInit(mem_, getPtr(qF_));
      if(flag != IDA_SUCCESS) idas_error("IDAQuadSensReInit",flag);
      log("IdasInternal::reset","re-initialized forward sensitivity dependent quadratures");
    }
  } else {
    // Turn of sensitivities
    flag = IDASensToggleOff(mem_);
    if(flag != IDA_SUCCESS) idas_error("IDASensToggleOff",flag);
  }

  // Correct initial conditions, if necessary
  int calc_ic = getOption("calc_ic");
  if(calc_ic){
    correctInitialConditions();
  }

  // Set the stop time of the integration -- don't integrate past this point
  if(stop_at_end_) setStopTime(tf_);
    
  log("IdasInternal::reset","end");
}
  
  
void IdasInternal::correctInitialConditions(){
  log("IdasInternal::correctInitialConditions","begin");
  if(monitored("correctInitialConditions")){
    cout << "initial guess: " << endl;
    cout << "p = " << input(INTEGRATOR_P) << endl;
    cout << "x0 = " << input(INTEGRATOR_X0) << endl;
    if(fsens_order_>0){
      for(int dir=0; dir<nfdir_; ++dir){
        cout << "forward seed guess, direction " << dir << ": " << endl;
        cout << "p_seed = " << fwdSeed(INTEGRATOR_P,dir) << endl;
        cout << "x0_seed = " << fwdSeed(INTEGRATOR_X0,dir) << endl;
      }
    }
  }

  int icopt = IDA_YA_YDP_INIT; // calculate z and xdot given x
  // int icopt = IDA_Y_INIT; // calculate z and x given zdot and xdot (e.g. start in stationary)

  double t_first = hasSetOption("first_time") ? double(getOption("first_time")) : tf_;
  int flag = IDACalcIC(mem_, icopt , t_first);
  if(flag != IDA_SUCCESS) idas_error("IDACalcIC",flag);

  // Retrieve the initial values
  flag = IDAGetConsistentIC(mem_, xz_, xzdot_);
  if(flag != IDA_SUCCESS) idas_error("IDAGetConsistentIC",flag);
  
  // Print progress
  log("IdasInternal::correctInitialConditions","found consistent initial values");
  if(monitored("correctInitialConditions")){
    cout << "p = " << input(INTEGRATOR_P) << endl;
    cout << "x0 = " << input(INTEGRATOR_X0) << endl;
    if(fsens_order_>0){
      for(int dir=0; dir<nfdir_; ++dir){
        cout << "forward seed, direction " << dir << ": " << endl;
        cout << "p_seed = " << fwdSeed(INTEGRATOR_P,dir) << endl;
        cout << "x0_seed = " << fwdSeed(INTEGRATOR_X0,dir) << endl;
      }
    }
  }
  log("IdasInternal::correctInitialConditions","end");
}
  
void IdasInternal::integrate(double t_out){
  log("IdasInternal::integrate","begin");
  int flag;
  
  // Check if we are already at the output time
  double ttol = 1e-9;   // tolerance
  if(fabs(t_-t_out)<ttol){
    // No integration necessary
    log("IdasInternal::integrate","already at the end of the horizon end");
    
  } else {
    // Integrate ...
    if(asens_order_>0){
      // ... with taping
      log("IdasInternal::integrate","integration with taping");
      flag = IDASolveF(mem_, t_out, &t_, xz_, xzdot_, IDA_NORMAL, &ncheck_);
      if(flag != IDA_SUCCESS && flag != IDA_TSTOP_RETURN) idas_error("IDASolveF",flag);
    } else {
      // ... without taping
      log("IdasInternal::integrate","integration without taping");
      flag = IDASolve(mem_, t_out, &t_, xz_, xzdot_, IDA_NORMAL);
      if(flag != IDA_SUCCESS && flag != IDA_TSTOP_RETURN) idas_error("IDASolve",flag);
    }
    log("IdasInternal::integrate","integration complete");
    
    // Get quadrature states
    if(nq_>0){
      double tret;
      flag = IDAGetQuad(mem_, &tret, q_);
      if(flag != IDA_SUCCESS) idas_error("IDAGetQuad",flag);
    }
    
    if(fsens_order_>0){
      // Get the sensitivities
      flag = IDAGetSens(mem_,&t_, getPtr(xzF_));
      if(flag != IDA_SUCCESS) idas_error("IDAGetSens",flag);
    
      if(nq_>0){
        double tret;
        flag = IDAGetQuadSens(mem_, &tret, getPtr(qF_));
        if(flag != IDA_SUCCESS) idas_error("IDAGetQuadSens",flag);
      }
    }
  }
  
  // Save the final state
  copy(NV_DATA_S(xz_),NV_DATA_S(xz_)+nx_,output(INTEGRATOR_XF).begin());
  if(fsens_order_>0){
    for(int i=0; i<nfdir_; ++i){
      copy(NV_DATA_S(xzF_[i]),NV_DATA_S(xzF_[i])+nx_,fwdSens(INTEGRATOR_XF,i).begin());
    }
  }
    
  // Print statistics
  if(getOption("print_stats")) printStats(std::cout);
  
  log("IdasInternal::integrate","end");
}

void IdasInternal::resetAdj(){
  int flag;
  
  // Reset adjoint sensitivities for the parameters
  N_VConst(0.0, rq_);
  
  // Get the adjoint seeds
  const Matrix<double> &xf_aseed = input(INTEGRATOR_RX0);
  copy(xf_aseed.begin(),xf_aseed.end(),NV_DATA_S(rxz_));
    
  if(isInitAdj_){
    flag = IDAReInitB(mem_, whichB_, tf_, rxz_, rxzdot_);
    if(flag != IDA_SUCCESS) idas_error("IDAReInitB",flag);
    
    if(np_>0){
      N_VConst(0.0,rq_);
      flag = IDAQuadReInit(IDAGetAdjIDABmem(mem_, whichB_),rq_);
      // flag = IDAQuadReInitB(mem_,whichB_[dir],rq_[dir]); // BUG in Sundials - do not use this!
    }
    if(flag!=IDA_SUCCESS) idas_error("IDAQuadReInitB",flag);
  } else {
    // Initialize the adjoint integration
    initAdj();
  }

  // Correct initial values for the integration if necessary
  int calc_icB = getOption("calc_icB");
  if(calc_icB){
    flag = IDACalcICB(mem_, whichB_, t0_, rxz_, rxzdot_);
    if(flag != IDA_SUCCESS) idas_error("IDACalcICB",flag);

    N_VPrint_Serial(rxz_);
    N_VPrint_Serial(rxzdot_);

    // Retrieve the initial values
    flag = IDAGetConsistentICB(mem_, whichB_, rxz_, rxzdot_);
    if(flag != IDA_SUCCESS) idas_error("IDAGetConsistentICB",flag);
    
    N_VPrint_Serial(rxz_);
    N_VPrint_Serial(rxzdot_);
  }
}

void IdasInternal::integrateAdj(double t_out){
  int flag;
  // Integrate backwards to t_out
  flag = IDASolveB(mem_, t_out, IDA_NORMAL);
  if(flag<IDA_SUCCESS) idas_error("IDASolveB",flag);

  // Get the sensitivities
  double tret;
  flag = IDAGetB(mem_, whichB_, &tret, rxz_, rxzdot_);
  if(flag!=IDA_SUCCESS) idas_error("IDAGetB",flag);

  flag = IDAGetQuadB(mem_, whichB_, &tret, rq_);
  if(flag!=IDA_SUCCESS) idas_error("IDAGetQuadB",flag);
  
  // Save the adjoint sensitivities
  const double *xz = NV_DATA_S(rxz_);
  copy(xz,xz+nx_,output(INTEGRATOR_RXF).begin());
}

void IdasInternal::printStats(std::ostream &stream) const{
  long nsteps, nfevals, nlinsetups, netfails;
  int qlast, qcur;
  double hinused, hlast, hcur, tcur;
  int flag = IDAGetIntegratorStats(mem_, &nsteps, &nfevals, &nlinsetups,&netfails, &qlast, &qcur, &hinused,&hlast, &hcur, &tcur);
  if(flag!=IDA_SUCCESS) idas_error("IDAGetIntegratorStats",flag);
  
  // Get the number of right hand side evaluations in the linear solver
  long nfevals_linsol=0;
  switch(linsol_f_){
    case SD_DENSE:
    case SD_BANDED:
      flag = IDADlsGetNumResEvals(mem_, &nfevals_linsol);
      if(flag!=IDA_SUCCESS) idas_error("IDADlsGetNumResEvals",flag);
      break;
    case SD_ITERATIVE:
      flag = IDASpilsGetNumResEvals(mem_, &nfevals_linsol);
      if(flag!=IDA_SUCCESS) idas_error("IDASpilsGetNumResEvals",flag);
      break;
    default:
      nfevals_linsol = 0;
  }
  
  stream << "number of steps taken by IDAS:            " << nsteps << std::endl;
  stream << "number of calls to the user's f function: " << (nfevals + nfevals_linsol) << std::endl;
  stream << "   step calculation:                      " << nfevals << std::endl;
  stream << "   linear solver:                         " << nfevals_linsol << std::endl;
  stream << "number of calls made to the linear solver setup function: " << nlinsetups << std::endl;
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
  stream << "Time spent in the jacobian function or jacobian times vector function: " << t_jac << " s." << endl;
  stream << "Time spent in the linear solver solve function: " << t_lsolve << " s." << endl;
  stream << "Time spent to generate the jacobian in the linear solver setup function: " << t_lsetup_jac << " s." << endl;
  stream << "Time spent to factorize the jacobian in the linear solver setup function: " << t_lsetup_fac << " s." << endl;
  stream << std::endl;
}

map<int,string> IdasInternal::calc_flagmap(){
  map<int,string> f;
  f[IDA_TSTOP_RETURN] = "IDA_TSTOP_RETURN";
  f[IDA_ROOT_RETURN] = "IDA_ROOT_RETURN";
  f[IDA_MEM_NULL] = "IDA_MEM_NULL";
  f[IDA_ILL_INPUT] = "IDA_ILL_INPUT";
  f[IDA_TOO_MUCH_WORK] = "IDA_TOO_MUCH_WORK";
  f[IDA_TOO_MUCH_ACC] = "IDA_TOO_MUCH_ACC";
  f[IDA_ERR_FAIL] = "IDA_ERR_FAIL";
  f[IDA_CONV_FAIL] = "IDA_CONV_FAIL";
  f[IDA_LINIT_FAIL] = "IDA_LINIT_FAIL";
  f[IDA_LSETUP_FAIL] = "IDA_LSETUP_FAIL";
  f[IDA_LSOLVE_FAIL] = "IDA_LSOLVE_FAIL";
  f[IDA_CONSTR_FAIL] = "IDA_CONSTR_FAIL";
  f[IDA_REP_RES_ERR] = "IDA_REP_RES_ERR";
  f[IDA_RES_FAIL] = "IDA_RES_FAIL";
  f[IDA_RTFUNC_FAIL] = "IDA_RTFUNC_FAIL";
  f[IDA_SUCCESS] = "IDA_SUCCESS";

  return f;
}
  
map<int,string> IdasInternal::flagmap = IdasInternal::calc_flagmap();

void IdasInternal::idas_error(const string& module, int flag){
  // Find the error
  map<int,string>::const_iterator it = flagmap.find(flag);
  const char* flagname = IDAGetReturnFlagName(flag);
  stringstream ss;
  ss << "Module \"" << module << "\" returned flag " << flag << " (\"" << flagname << "\").";
  ss << " Consult Idas documentation." << std::endl;
  delete flagname;
  
  // Heuristics
  if (
    (module=="IDACalcIC" && (flag==IDA_CONV_FAIL || flag==IDA_NO_RECOVERY || flag==IDA_LINESEARCH_FAIL )) ||
    (module=="IDASolve" && flag ==IDA_ERR_FAIL )
    ) {
    ss << "Some common causes for this error: " << std::endl;
    ss << "  - providing an initial guess for which 0=g(y,z,t) is not invertible wrt y. " << std::endl;
    ss << "  - having a DAE-index higher than 1 such that 0=g(y,z,t) is not invertible wrt y over the whole domain." << std::endl;
    ss << "  - having set abstol or reltol too small." << std::endl;
    ss << "  - using 'calcic'=True for systems that are not semi-explicit index-one. You must provide consistent initial conditions yourself in this case. " << std::endl;
    ss << "  - your problem is too hard for IDAcalcIC to solve. Provide consistent initial conditions yourself." << std::endl;
  }
  
  casadi_error(ss.str());
}

int IdasInternal::rhsQ_wrapper(double t, N_Vector xz, N_Vector xzdot, N_Vector rhsQ, void *user_data){
 try{
    IdasInternal *this_ = (IdasInternal*)user_data;
    this_->rhsQ(t,NV_DATA_S(xz),NV_DATA_S(xzdot),NV_DATA_S(rhsQ));
    return 0;
  } catch(exception& e){
    cerr << "rhsQ failed: " << e.what() << endl;
    return 1;
  }
}

void IdasInternal::rhsQ(double t, const double* xz, const double* xzdot, double* rhsQ){
   log("IdasInternal::rhsQ","begin");
   // Pass input
   f_.setInput(&t,DAE_T);
   f_.setInput(xz,DAE_X);
   f_.setInput(xz+nx_,DAE_Z);
   f_.setInput(xzdot,DAE_XDOT);
   f_.setInput(input(INTEGRATOR_P),DAE_P);

    // Evaluate
   f_.evaluate();
    
    // Get results
   f_.getOutput(rhsQ,DAE_QUAD);
   log("IdasInternal::rhsQ","end");
}
  
void IdasInternal::rhsQS(int Ns, double t, N_Vector xz, N_Vector xzdot, N_Vector *xzF, N_Vector *xzdotF, N_Vector rrQ, N_Vector *qdotF, 
                        N_Vector tmp1, N_Vector tmp2, N_Vector tmp3){

  log("IdasInternal::rhsQS","enter");
  casadi_assert(Ns==nfdir_);
  // Pass input
   f_.setInput(&t,DAE_T);
   f_.setInput(NV_DATA_S(xz),DAE_X);
   f_.setInput(NV_DATA_S(xz)+nx_,DAE_Z);
   f_.setInput(NV_DATA_S(xzdot),DAE_XDOT);
   f_.setInput(input(INTEGRATOR_P),DAE_P);
     
   // Pass forward seeds
  for(int i=0; i<nfdir_; ++i){
    f_.fwdSeed(DAE_T).setZero();
    f_.setFwdSeed(NV_DATA_S(xzF[i]),DAE_X);
    f_.setFwdSeed(NV_DATA_S(xzF[i])+nx_,DAE_Z);
    f_.setFwdSeed(NV_DATA_S(xzdotF[i]),DAE_XDOT);
    f_.setFwdSeed(fwdSeed(INTEGRATOR_P,i),DAE_P);
   
    // Evaluate the AD forward algorithm
    f_.evaluate(1,0);
      
    // Get the output seeds
    f_.getFwdSens(NV_DATA_S(qdotF[i]),DAE_QUAD);
  }
  log("IdasInternal::rhsQS","end");
}

int IdasInternal::rhsQS_wrapper(int Ns, double t, N_Vector xz, N_Vector xzdot, N_Vector *xzF, N_Vector *xzdotF, N_Vector rrQ, N_Vector *qdotF, void *user_data, 
                                N_Vector tmp1, N_Vector tmp2, N_Vector tmp3){
  
 try{
    IdasInternal *this_ = (IdasInternal*)user_data;
    this_->rhsQS(Ns,t,xz,xzdot,xzF,xzdotF,rrQ,qdotF,tmp1,tmp2,tmp3);
    return 0;
  } catch(exception& e){
    cerr << "rhsQS failed: " << e.what() << endl;
    return 1;
  }
}

void IdasInternal::resB(double t, const double* xz, const double* xzdot, const double* xzA, const double* xzdotA, double* rrA){
  log("IdasInternal::resB","begin");
  
  // Pass inputs
  g_.setInput(&t,RDAE_T);
  g_.setInput(xz,RDAE_X);
  g_.setInput(xz+nx_,RDAE_Z);
  g_.setInput(xzdot,RDAE_XDOT);
  g_.setInput(input(INTEGRATOR_P),RDAE_P);
  g_.setInput(xzA,RDAE_RX);
  g_.setInput(xzA+nx_,RDAE_RZ);
  g_.setInput(xzdotA,RDAE_RXDOT);
  
  // Evaluate
  g_.evaluate();

  // Save to output
  g_.getOutput(rrA,RDAE_ODE);
  g_.getOutput(rrA+nx_,RDAE_ALG);
  
  // Negate as we are integrating backwards in time
  for(int i=0; i<nx_+nz_; ++i)
    rrA[i] *= -1;

  log("IdasInternal::resB","end");
}

int IdasInternal::resB_wrapper(double t, N_Vector xz, N_Vector xzdot, N_Vector xzA, N_Vector xzdotA, N_Vector rrA, void *user_dataB){
 try{
    IdasInternal *this_ = (IdasInternal*)user_dataB;
    this_->resB(t,NV_DATA_S(xz),NV_DATA_S(xzdot),NV_DATA_S(xzA),NV_DATA_S(xzdotA),NV_DATA_S(rrA));
    return 0;
  } catch(exception& e){
    cerr << "resB failed: " << e.what() << endl;
    return 1;
  }
}
  
void IdasInternal::rhsQB(double t, const double* xz, const double* xzdot, const double* xzA, const double* xzdotA, double *qdotA){
  log("IdasInternal::rhsQB","begin");

    // Pass inputs
  g_.setInput(&t,RDAE_T);
  g_.setInput(xz,RDAE_X);
  g_.setInput(xz+nx_,RDAE_Z);
  g_.setInput(xzdot,RDAE_XDOT);
  g_.setInput(input(INTEGRATOR_P),RDAE_P);
  g_.setInput(xzA,RDAE_RX);
  g_.setInput(xzA+nx_,RDAE_RZ);
  g_.setInput(xzdotA,RDAE_RXDOT);
  
  // Evaluate
  g_.evaluate();

  // Save to output
  g_.getOutput(qdotA,RDAE_QUAD);
  
  // Negate as we are integrating backwards in time
  for(int i=0; i<nrq_; ++i)
    qdotA[i] *= -1;
  
  log("IdasInternal::rhsQB","end");
}

int IdasInternal::rhsQB_wrapper(double t, N_Vector y, N_Vector xzdot, N_Vector xzA, N_Vector xzdotA, N_Vector qdotA, void *user_dataB){
 try{
    IdasInternal *this_ = (IdasInternal*)user_dataB;
    this_->rhsQB(t,NV_DATA_S(y),NV_DATA_S(xzdot),NV_DATA_S(xzA),NV_DATA_S(xzdotA),NV_DATA_S(qdotA));
    return 0;
  } catch(exception& e){
    cerr << "resQB failed: " << e.what() << endl;
    return 1;
  }
}

void IdasInternal::djac(SUNDIALS_INT Neq, double t, double cj, N_Vector xz, N_Vector xzdot, N_Vector rr, DlsMat Jac, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3){
  log("IdasInternal::djac","begin");

  // Get time
  time1 = clock();
  
  // Pass input to the Jacobian function
  jac_.setInput(&t,JAC_T);
  jac_.setInput(NV_DATA_S(xz),JAC_X);
  jac_.setInput(NV_DATA_S(xz)+nx_,JAC_Z);
  jac_.setInput(NV_DATA_S(xzdot),JAC_XDOT);
  jac_.setInput(input(INTEGRATOR_P),JAC_P);
  jac_.setInput(cj,JAC_CJ);
    
  // Evaluate Jacobian
  jac_.evaluate();

  // Get sparsity and non-zero elements
  const vector<int>& rowind = jac_.output().rowind();
  const vector<int>& col = jac_.output().col();
  const vector<double>& val = jac_.output().data();

  // Dimension of the Jacobian
  //int jdim = jac_.output().size1();
    
  // Loop over rows
  for(int i=0; i<rowind.size()-1; ++i){
    // Loop over non-zero entries
    for(int el=rowind[i]; el<rowind[i+1]; ++el){
      
      // Get column
      int j = col[el];
      
      // Add to the element
      DENSE_ELEM(Jac,i,j) = val[el];
    }
  }
  
  // Log time duration
  time2 = clock();
  t_jac += double(time2-time1)/CLOCKS_PER_SEC;
  log("IdasInternal::djac","end");
}

int IdasInternal::djac_wrapper(SUNDIALS_INT Neq, double t, double cj, N_Vector xz, N_Vector xzdot, N_Vector rr, DlsMat Jac, void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3){
 try{
    IdasInternal *this_ = (IdasInternal*)user_data;
    this_->djac(Neq, t, cj, xz, xzdot, rr, Jac, tmp1, tmp2, tmp3);
    return 0;
  } catch(exception& e){
    cerr << "djac failed: " << e.what() << endl;
    return 1;
  }
}

void IdasInternal::bjac(SUNDIALS_INT Neq, SUNDIALS_INT mupper, SUNDIALS_INT mlower, double tt, double cj, N_Vector xz, N_Vector xzdot, N_Vector rr, DlsMat Jac, N_Vector tmp1, N_Vector tmp2,N_Vector tmp3){
  log("IdasInternal::bjac","begin");
  // Get time
  time1 = clock();

  // Pass input to the jacobian function
  jac_.setInput(tt,JAC_T);
  jac_.setInput(NV_DATA_S(xz),JAC_X);
  jac_.setInput(NV_DATA_S(xz)+nx_,JAC_Z);
  jac_.setInput(NV_DATA_S(xzdot),JAC_XDOT);
  jac_.setInput(input(INTEGRATOR_P),JAC_P);
  jac_.setInput(cj,JAC_CJ);

  // Evaluate jacobian
  jac_.evaluate();

  // Get sparsity and non-zero elements
  const vector<int>& rowind = jac_.output().rowind();
  const vector<int>& col = jac_.output().col();
  const vector<double>& val = jac_.output().data();

  // Loop over rows
  for(int i=0; i<rowind.size()-1; ++i){
    // Loop over non-zero entries
    for(int el=rowind[i]; el<rowind[i+1]; ++el){
      // Get column
      int j = col[el];
      
      // Set the element
      if(i-j>=-mupper && i-j<=mlower)
        BAND_ELEM(Jac,i,j) = val[el];
    }
  }
  
  // Log time duration
  time2 = clock();
  t_jac += double(time2-time1)/CLOCKS_PER_SEC;
  log("IdasInternal::bjac","end");
}

int IdasInternal::bjac_wrapper(SUNDIALS_INT Neq, SUNDIALS_INT mupper, SUNDIALS_INT mlower, double tt, double cj, N_Vector xz, N_Vector xzdot, N_Vector rr, DlsMat Jac, void *user_data, N_Vector tmp1, N_Vector tmp2,N_Vector tmp3){
 try{
    IdasInternal *this_ = (IdasInternal*)user_data;
    this_->bjac(Neq, mupper, mlower, tt, cj, xz, xzdot, rr, Jac, tmp1, tmp2, tmp3);
    return 0;
  } catch(exception& e){
    cerr << "bjac failed: " << e.what() << endl;
    return 1;
  }
}

void IdasInternal::setStopTime(double tf){
  // Set the stop time of the integration -- don't integrate past this point
  int flag = IDASetStopTime(mem_, tf);
  if(flag != IDA_SUCCESS) idas_error("IDASetStopTime",flag);
}

int IdasInternal::psolve_wrapper(double t, N_Vector xz, N_Vector xzdot, N_Vector rr, N_Vector rvec, N_Vector zvec, double cj, double delta, void *user_data, N_Vector tmp){
 try{
    IdasInternal *this_ = (IdasInternal*)user_data;
    casadi_assert(this_);
    this_->psolve(t, xz, xzdot, rr, rvec, zvec, cj, delta, tmp);
    return 0;
  } catch(exception& e){
    cerr << "psolve failed: " << e.what() << endl;
    return 1;
  }
}

int IdasInternal::psetup_wrapper(double t, N_Vector xz, N_Vector xzdot, N_Vector rr, double cj, void* user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3){
 try{
    IdasInternal *this_ = (IdasInternal*)user_data;
    casadi_assert(this_);
    this_->psetup(t, xz, xzdot, rr, cj, tmp1, tmp2, tmp3);
    return 0;
  } catch(exception& e){
    cerr << "psetup failed: " << e.what() << endl;
    return 1;
  }
}

void IdasInternal::psolve(double t, N_Vector xz, N_Vector xzdot, N_Vector rr, N_Vector rvec, N_Vector zvec, double cj, double delta, N_Vector tmp){
  log("IdasInternal::psolve","begin");
  
  // Get time
  time1 = clock();

  // Copy input to output, if necessary
  if(rvec!=zvec){
    N_VScale(1.0, rvec, zvec);
  }
  
  // Solve the (possibly factorized) system 
  casadi_assert(linsol_.output().size()*nrhs_ == NV_LENGTH_S(zvec));
  linsol_.solve(NV_DATA_S(zvec),nrhs_);

  // Log time duration
  time2 = clock();
  t_lsolve += double(time2-time1)/CLOCKS_PER_SEC;
  log("IdasInternal::psolve","end");
}

void IdasInternal::psetup(double t, N_Vector xz, N_Vector xzdot, N_Vector rr, double cj, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3){
  log("IdasInternal::psetup","begin");
  
  // Get time
  time1 = clock();

  // Pass input to the jacobian function
  jac_.setInput(t,JAC_T);
  jac_.setInput(NV_DATA_S(xz),JAC_X);
  jac_.setInput(NV_DATA_S(xz)+nx_,JAC_Z);
  jac_.setInput(NV_DATA_S(xzdot),JAC_XDOT);
  jac_.setInput(input(INTEGRATOR_P),JAC_P);
  jac_.setInput(cj,JAC_CJ);

  // Evaluate jacobian
  jac_.evaluate();

  // Log time duration
  time2 = clock();
  t_lsetup_jac += double(time2-time1)/CLOCKS_PER_SEC;

  // Pass non-zero elements to the linear solver
  linsol_.setInput(jac_.output(),0);

  // Prepare the solution of the linear system (e.g. factorize) -- only if the linear solver inherits from LinearSolver
  linsol_.prepare();

  // Log time duration
  time1 = clock();
  t_lsetup_fac += double(time1-time2)/CLOCKS_PER_SEC;

  log("IdasInternal::psetup","end");
  
}

int IdasInternal::lsetup_wrapper(IDAMem IDA_mem, N_Vector xz, N_Vector xzdot, N_Vector resp, N_Vector vtemp1, N_Vector vtemp2, N_Vector vtemp3){
 try{
    IdasInternal *this_ = (IdasInternal*)(IDA_mem->ida_lmem);
    casadi_assert(this_);
    this_->lsetup(IDA_mem,xz,xzdot,resp,vtemp1,vtemp2,vtemp3);
    return 0;
  } catch(exception& e){
    cerr << "lsetup failed: " << e.what() << endl;
    return -1;
  }
}

int IdasInternal::lsolve_wrapper(IDAMem IDA_mem, N_Vector b, N_Vector weight, N_Vector xz, N_Vector xzdot, N_Vector rr){
 try{
   IdasInternal *this_ = (IdasInternal*)(IDA_mem->ida_lmem);
   casadi_assert(this_);
   this_->lsolve(IDA_mem,b,weight,xz,xzdot,rr);
   return 0;
  } catch(int wrn){
/*    cerr << "warning: " << wrn << endl;*/
    return wrn;
  } catch(exception& e){
    cerr << "lsolve failed: " << e.what() << endl;
    return -1;
  }
}

void IdasInternal::lsetup(IDAMem IDA_mem, N_Vector xz, N_Vector xzdot, N_Vector resp, N_Vector vtemp1, N_Vector vtemp2, N_Vector vtemp3){  
  log("IdasInternal::lsetup","begin");
  
  // Current time
  double t = IDA_mem->ida_tn;

  // Multiple of df_dydot to be added to the matrix
  double cj = IDA_mem->ida_cj;

  // Call the preconditioner setup function (which sets up the linear solver)
  psetup(t, xz, xzdot, 0, cj, vtemp1, vtemp1, vtemp3);
  log("IdasInternal::lsetup","end");
}

void IdasInternal::lsolve(IDAMem IDA_mem, N_Vector b, N_Vector weight, N_Vector xz, N_Vector xzdot, N_Vector rr){
  log("IdasInternal::lsolve","begin");
  // Current time
  double t = IDA_mem->ida_tn;

  // Multiple of df_dydot to be added to the matrix
  double cj = IDA_mem->ida_cj;

  // Accuracy
  double delta = 0.0;
  
  // Call the preconditioner solve function (which solves the linear system)
  psolve(t, xz, xzdot, rr, b, b, cj, delta, 0);
  
  // Scale the correction to account for change in cj
  if(cj_scaling_){
    double cjratio = IDA_mem->ida_cjratio;
    if (cjratio != 1.0) N_VScale(2.0/(1.0 + cjratio), b, b);
  }
  log("IdasInternal::lsolve","end");
}

void IdasInternal::initDenseLinearSolver(){
  // Dense jacobian
  int flag = IDADense(mem_, nx_+nz_);
  if(flag != IDA_SUCCESS) idas_error("IDADense",flag);
  if(exact_jacobian_){
    // Generate jacobians if not already provided
    if(jac_.isNull()){
      getJacobian();
      jac_.init();
    }
    
    // Pass to IDA
    flag = IDADlsSetDenseJacFn(mem_, djac_wrapper);
    if(flag!=IDA_SUCCESS) idas_error("IDADlsSetDenseJacFn",flag);
  }
}
  
void IdasInternal::initBandedLinearSolver(){
  // Banded jacobian
  int flag = IDABand(mem_, nx_+nz_, getOption("upper_bandwidth").toInt(), getOption("lower_bandwidth").toInt());
  if(flag != IDA_SUCCESS) idas_error("IDABand",flag);
  
  // Banded Jacobian information
  if(exact_jacobian_){
    flag = IDADlsSetBandJacFn(mem_, bjac_wrapper);
    if(flag != IDA_SUCCESS) idas_error("IDADlsSetBandJacFn",flag);
  }
}
  
void IdasInternal::initIterativeLinearSolver(){
  // Max dimension of the Krylov space
  int maxl = getOption("max_krylov");
      
  // Attach an iterative solver
  int flag;
  switch(itsol_f_){
    case SD_GMRES:
      flag = IDASpgmr(mem_, maxl);
      if(flag != IDA_SUCCESS) idas_error("IDASpgmr",flag);
      break;
    case SD_BCGSTAB:
      flag = IDASpbcg(mem_, maxl);
      if(flag != IDA_SUCCESS) idas_error("IDASpbcg",flag);
      break;
    case SD_TFQMR:
      flag = IDASptfqmr(mem_, maxl);
      if(flag != IDA_SUCCESS) idas_error("IDASptfqmr",flag);
      break;
  }
        
  // Attach functions for jacobian information
  if(exact_jacobian_){
    flag = IDASpilsSetJacTimesVecFn(mem_, jtimes_wrapper);
    if(flag != IDA_SUCCESS) idas_error("IDASpilsSetJacTimesVecFn",flag);
  }
  
  // Add a preconditioner
  if(use_preconditioner_){
    // Make sure that a Jacobian has been provided
    if(jac_.isNull()) throw CasadiException("IdasInternal::init(): No Jacobian has been provided.");

    // Make sure that a linear solver has been providided
    if(linsol_.isNull()) throw CasadiException("IdasInternal::init(): No user defined linear solver has been provided.");

    // Pass to IDA
    flag = IDASpilsSetPreconditioner(mem_, psetup_wrapper, psolve_wrapper);
    if(flag != IDA_SUCCESS) idas_error("IDASpilsSetPreconditioner",flag);
  }
}

void IdasInternal::initUserDefinedLinearSolver(){
  // Make sure that a Jacobian has been provided
  casadi_assert(!jac_.isNull());

  // Make sure that a linear solver has been providided
  casadi_assert(!linsol_.isNull());

  //  Set fields in the IDA memory
  IDAMem IDA_mem = IDAMem(mem_);
  IDA_mem->ida_lmem   = this;
  IDA_mem->ida_lsetup = lsetup_wrapper;
  IDA_mem->ida_lsolve = lsolve_wrapper;
  IDA_mem->ida_setupNonNull = TRUE;
}

void IdasInternal::initDenseLinearSolverB(){
  int flag = IDADenseB(mem_, whichB_, nrx_+nrz_);
  if(flag != IDA_SUCCESS) idas_error("IDADenseB",flag);
}
  
void IdasInternal::initBandedLinearSolverB(){
  int flag = IDABandB(mem_, whichB_, nrx_+nrz_, getOption("asens_upper_bandwidth").toInt(), getOption("asens_lower_bandwidth").toInt());
  if(flag != IDA_SUCCESS) idas_error("IDABand",flag);
}
  
void IdasInternal::initIterativeLinearSolverB(){
  int maxl = getOption("asens_max_krylov");
  int flag;
  switch(itsol_g_){
    case SD_GMRES:
      flag = IDASpgmrB(mem_, whichB_, maxl);
      if(flag != IDA_SUCCESS) idas_error("IDASpgmrB",flag);
      break;
    case SD_BCGSTAB:
      flag = IDASpbcgB(mem_, whichB_, maxl);
      if(flag != IDA_SUCCESS) idas_error("IDASpbcgB",flag);
      break;
    case SD_TFQMR:
      flag = IDASptfqmrB(mem_, whichB_, maxl);
      if(flag != IDA_SUCCESS) idas_error("IDASptfqmrB",flag);
      break;
  }
}
  
void IdasInternal::initUserDefinedLinearSolverB(){
  casadi_assert_message(false, "Not implemented");
}

void IdasInternal::setLinearSolver(const LinearSolver& linsol, const FX& jac){
  linsol_ = linsol;
  jac_ = jac;

  // Try to generate a jacobian of none provided
  if(jac_.isNull())
    getJacobian();
}

FX IdasInternal::getJacobian(){
  // Quick return if already created
  if(!jac_.isNull())
    return jac_;

  // If SXFunction
  SXFunction f_sx = shared_cast<SXFunction>(f_);
  if(!f_sx.isNull()){
    // Get the Jacobian in the Newton iteration
    SX cj("cj");
    SXMatrix jac_ode_x = f_sx.jac(DAE_X,DAE_ODE);
    SXMatrix jac_ode_xdot = f_sx.jac(DAE_XDOT,DAE_ODE);
    SXMatrix jac = jac_ode_x + cj*jac_ode_xdot;
    if(nz_>0){
      SXMatrix jac_alg_x = f_sx.jac(DAE_X,DAE_ALG);
      SXMatrix jac_alg_z = f_sx.jac(DAE_Z,DAE_ALG);
      SXMatrix jac_ode_z = f_sx.jac(DAE_Z,DAE_ODE);
      jac = horzcat(vertcat(jac,jac_alg_x),vertcat(jac_ode_z,jac_alg_z));
    }
    
    // Jacobian function
    vector<Matrix<SX> > jac_in(JAC_NUM_IN);
    jac_in[JAC_T] = f_sx.inputSX(DAE_T);
    jac_in[JAC_X] = f_sx.inputSX(DAE_X);
    jac_in[JAC_Z] = f_sx.inputSX(DAE_Z);
    jac_in[JAC_XDOT] = f_sx.inputSX(DAE_XDOT);
    jac_in[JAC_P] = f_sx.inputSX(DAE_P);
    jac_in[JAC_CJ] = cj;
    SXFunction J(jac_in,jac);
    
    // Save function
    jac_ = J;
    return J;
  }

  // If SXFunction
  MXFunction f_mx = shared_cast<MXFunction>(f_);
  if(!f_mx.isNull()){
    // Get the Jacobian in the Newton iteration
    MX cj("cj");
    MX jac = f_mx.jac(DAE_X,DAE_ODE) + cj*f_mx.jac(DAE_XDOT,DAE_ODE);
    if(nz_>0){
      jac = horzcat(vertcat(jac,f_mx.jac(DAE_X,DAE_ALG)),vertcat(f_mx.jac(DAE_Z,DAE_ODE),f_mx.jac(DAE_Z,DAE_ALG)));
    }
    
    // Jacobian function
    vector<MX> jac_in(JAC_NUM_IN);
    jac_in[JAC_T] = f_mx.inputMX(DAE_T);
    jac_in[JAC_X] = f_mx.inputMX(DAE_X);
    jac_in[JAC_Z] = f_mx.inputMX(DAE_Z);
    jac_in[JAC_XDOT] = f_mx.inputMX(DAE_XDOT);
    jac_in[JAC_P] = f_mx.inputMX(DAE_P);
    jac_in[JAC_CJ] = cj;
    MXFunction J(jac_in,jac);
    
    // Save function
    jac_ = J;
    return J;
  }
  
  throw CasadiException("IdasInternal::getJacobian(): Not an SXFunction or MXFunction");
}

  
LinearSolver IdasInternal::getLinearSolver(){
  return linsol_;
}

} // namespace Sundials
} // namespace CasADi

