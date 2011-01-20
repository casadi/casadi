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
#include "casadi/fx/sx_function_internal.hpp"
#include "casadi/sx/sx_tools.hpp"
#include <cassert>

using namespace std;
namespace CasADi{
namespace Sundials{

IdasInternal* IdasInternal::clone() const{
  // Copying initialized objects are not allowed since they contain pointers
  if(is_init) throw CasadiException("IdasInternal::clone: cannot clone an initialized object");
  
  // Return a deep copy
  return new IdasInternal(*this);
}

// IdasInternal::IdasInternal(const IdasInternal& integrator): IntegratorInternal(integrator){
//   f_ = integrator_.f_;
//   q_ = integrator_.q_;
//   
// }

IdasInternal::IdasInternal(const FX& f, const FX& q) : f_(f), q_(q){
  // Check dimensions
  if(f.getNumInputs()!=DAE_NUM_IN) throw CasadiException("IdasInternal: f has wrong number of inputs");
  if(f.getNumOutputs()!=DAE_NUM_OUT) throw CasadiException("IdasInternal: f has wrong number of outputs");
  if(!q.isNull()){
    if(q.getNumInputs()!=DAE_NUM_IN) throw CasadiException("IdasInternal: q has wrong number of inputs");
    if(q.getNumOutputs()!=DAE_NUM_OUT) throw CasadiException("IdasInternal: q has wrong number of outputs");
  }

  addOption("suppress_algebraic",          OT_BOOLEAN, false); // supress algebraic variables in the error testing
  addOption("calc_ic",                     OT_BOOLEAN, true);  // use IDACalcIC to get consistent initial conditions
  addOption("calc_icB",                    OT_BOOLEAN, false);  // use IDACalcIC to get consistent initial conditions
  addOption("abstolv",                     OT_REALVECTOR);
  addOption("fsens_abstolv",               OT_REALVECTOR); 
  addOption("max_step_size",               OT_REAL, 0); // maximim step size
  addOption("first_time",                  OT_REAL); // first requested time as a fraction of the time interval
  addOption("cj_scaling",                  OT_BOOLEAN, false);  // IDAS scaling on cj for the user-defined linear solver module
  addOption("extra_fsens_calc_ic",         OT_BOOLEAN, false);  // Call calc ic an extra time, with fsens=0
  
  
  mem_ = 0;
  
  yz_  = 0; 
  yP_ = 0, 
  yQ_ = 0;
  id_ = 0;

  is_init = false;

  ny_ = f.input(DAE_Y).numel();
  nq_ = q.isNull() ? 0 : q.output().numel();
  int np = f.input(DAE_P).numel();
  int nz = f.input(DAE_Z).numel();
  setDimensions(ny_+nq_,np,nz);
  nyz_ = ny_ + nz;
  ncheck_ = 0;
}

IdasInternal::~IdasInternal(){ 
  if(mem_) IDAFree(&mem_);

  // N-vectors for the DAE integration
  if(yz_) N_VDestroy_Serial(yz_);
  if(yP_) N_VDestroy_Serial(yP_);
  if(yQ_) N_VDestroy_Serial(yQ_);
  if(id_) N_VDestroy_Serial(id_);
  
    // Forward problem
  for(vector<N_Vector>::iterator it=yzS_.begin(); it != yzS_.end(); ++it)   if(*it) N_VDestroy_Serial(*it);
  for(vector<N_Vector>::iterator it=yPS_.begin(); it != yPS_.end(); ++it)   if(*it) N_VDestroy_Serial(*it);
  for(vector<N_Vector>::iterator it=yQS_.begin(); it != yQS_.end(); ++it)   if(*it) N_VDestroy_Serial(*it);

  // Adjoint problem
  for(vector<N_Vector>::iterator it=yzB_.begin(); it != yzB_.end(); ++it)   if(*it) N_VDestroy_Serial(*it);
  for(vector<N_Vector>::iterator it=yPB_.begin(); it != yPB_.end(); ++it)   if(*it) N_VDestroy_Serial(*it);
  for(vector<N_Vector>::iterator it=yBB_.begin(); it != yBB_.end(); ++it)   if(*it) N_VDestroy_Serial(*it);

}

void IdasInternal::init(){
  // Call the base class init
  IntegratorInternal::init();
  log("IdasInternal::init","begin");
  
  // Print
  if(verbose()){
    cout << "Initializing IDAS with ny_ = " << ny_ << ", nq_ = " << nq_ << ", np_ = " << np_ << " and nz_ = " << nz_ << endl;
  }
  
  // Init ODE rhs function and quadrature functions, jacobian function
  f_.init();
  if(!q_.isNull()) q_.init();
  
  log("IdasInternal::init","functions initialized");
  
  if(!jac_.isNull()){
    jac_.init();
    vector<int> rowind, col;
    jac_.output().sparsity().getSparsityCRS(rowind,col);
    if(!linsol_.isNull())
      linsol_.setSparsity(rowind,col);
      linsol_.init();
    
    log("IdasInternal::init","user defined linear solver initialized");
  }

  // Get the number of forward and adjoint directions
  nfdir_f_ = f_.getOption("number_of_fwd_dir").toInt();
  nadir_f_ = f_.getOption("number_of_adj_dir").toInt();
  nfdir_q_ = q_.isNull() ? 0 : q_.getOption("number_of_fwd_dir").toInt();
  nadir_q_ = q_.isNull() ? 0 : q_.getOption("number_of_adj_dir").toInt();

  // Quick return if already initialized
  if(is_init){
    reset(ad_order_,ad_order_);
    log("IdasInternal::init","end, Idas already initialized");
    return;
  }

  cj_scaling_ = getOption("cj_scaling").toInt();
  
  // Sundials return flag
  int flag;

  // Create IDAS memory block
  mem_ = IDACreate();
  if(mem_==0) throw CasadiException("IDACreate(): Creation failed");

  // Allocate n-vectors for ivp
  yz_ = N_VNew_Serial(nyz_);
  yP_ = N_VNew_Serial(nyz_);
  id_ = N_VNew_Serial(nyz_);

  // Initialize Idas
  double t0 = 0;
  N_VConst(0.0, yz_);
  N_VConst(0.0, yP_);
  IDAInit(mem_, res_wrapper, t0, yz_, yP_);
  log("IdasInternal::init","IDA initialized");

  // Set error handler function
  flag = IDASetErrHandlerFn(mem_, ehfun_wrapper, this);
  if(flag != IDA_SUCCESS) idas_error("IDASetErrHandlerFn",flag);

  // Include algebraic variables in error testing
  flag = IDASetSuppressAlg(mem_, getOption("suppress_algebraic").toInt());
  if(flag != IDA_SUCCESS) idas_error("IDASetSuppressAlg",flag);

  // Maxinum order for the multistep method
  flag = IDASetMaxOrd(mem_, getOption("max_multistep_order").toInt());
  if(flag != IDA_SUCCESS) idas_error("IDASetMaxOrd",flag);

  // Set user data
  flag = IDASetUserData(mem_,this);
  if(flag != IDA_SUCCESS) idas_error("IDASetUserData",flag);

  // Set maximum step size
  flag = IDASetMaxStep(mem_, getOption("max_step_size").toDouble());
  if(flag != IDA_SUCCESS) idas_error("IDASetMaxStep",flag);
  
  if(hasSetOption("abstolv")){
    // Vector absolute tolerances
    vector<double> abstolv = getOption("abstolv").toDoubleVector();
    N_Vector nv_abstol = N_VMake_Serial(abstolv.size(),&abstolv[0]);
    flag = IDASVtolerances(mem_, reltol_, nv_abstol);
    if(flag != IDA_SUCCESS) idas_error("IDASVtolerances",flag);
    N_VDestroy_Serial(nv_abstol);
  } else {
    // Scalar absolute tolerances
    flag = IDASStolerances(mem_, reltol_, abstol_); 
    if(flag != IDA_SUCCESS) idas_error("IDASStolerances",flag);
  }
  
  // Maximum number of steps
  IDASetMaxNumSteps(mem_, getOption("max_num_steps").toInt());
  if(flag != IDA_SUCCESS) idas_error("IDASetMaxNumSteps",flag);

  // Set algebraic components
  if(hasSetOption("is_differential")){
    // User-specified differential components
    const vector<int>& is_differential = getOption("is_differential").toIntVector();
    copy(is_differential.begin(),is_differential.end(),NV_DATA_S(id_));
  } else {
    // By default, the ny_ first components are differential, the nz_ following algebraic
    fill_n(NV_DATA_S(id_),ny_,1);
    fill_n(NV_DATA_S(id_)+ny_,nz_,0);
  }
  
  // Pass this information to IDAS
  flag = IDASetId(mem_, id_);
  if(flag != IDA_SUCCESS) idas_error("IDASetId",flag);
    
  // attach a linear solver
  if(getOption("linear_solver")=="dense"){
    // Dense jacobian
    flag = IDADense(mem_, nyz_);
    if(flag != IDA_SUCCESS) idas_error("IDADense",flag);
    if(exact_jacobian_){
      // Generate jacobians if not already provided
      if(jac_.isNull()){
        if(jacx_.isNull()) jacx_ = f_.jacobian(DAE_Y,DAE_RES);
        if(jacxdot_.isNull()) jacxdot_ = f_.jacobian(DAE_YDOT,DAE_RES);
        if(jacz_.isNull()) jacz_ = f_.jacobian(DAE_Z,DAE_RES);
        jacx_.init();
        jacxdot_.init();
        jacz_.init();
      }
      
      // Pass to IDA
      flag = IDADlsSetDenseJacFn(mem_, djac_wrapper);
      if(flag!=IDA_SUCCESS) idas_error("IDADlsSetDenseJacFn",flag);
    }
  } else if(getOption("linear_solver")=="banded") {
    // Banded jacobian
    flag = IDABand(mem_, nyz_, getOption("upper_bandwidth").toInt(), getOption("lower_bandwidth").toInt());
    if(flag != IDA_SUCCESS) idas_error("IDABand",flag);
    
    // Banded Jacobian information
    if(exact_jacobian_){
      flag = IDADlsSetBandJacFn(mem_, bjac_wrapper);
      if(flag != IDA_SUCCESS) idas_error("IDADlsSetBandJacFn",flag);
    }
  } else if(getOption("linear_solver")=="iterative") {
    // Max dimension of the Krylov space
    int maxl = getOption("max_krylov").toInt();
    
    // Attach an iterative solver
    if(getOption("iterative_solver")=="gmres"){
      flag = IDASpgmr(mem_, maxl);
      if(flag != IDA_SUCCESS) idas_error("IDASpgmr",flag);
    } else if(getOption("iterative_solver")=="bcgstab"){
      flag = IDASpbcg(mem_, maxl);
      if(flag != IDA_SUCCESS) idas_error("IDASpbcg",flag);
    } else if(getOption("iterative_solver")=="tfqmr") {
      flag = IDASptfqmr(mem_, maxl);
      if(flag != IDA_SUCCESS) idas_error("IDASptfqmr",flag);
    } else throw CasadiException("Unknown sparse solver");
       
    // Attach functions for jacobian information
    if(exact_jacobian_){
      flag = IDASpilsSetJacTimesVecFn(mem_, jtimes_wrapper);
      if(flag != IDA_SUCCESS) idas_error("IDASpilsSetJacTimesVecFn",flag);
    }
    
    // Add a preconditioner
    if(getOption("use_preconditioner")==true){
      // Make sure that a Jacobian has been provided
      if(jac_.isNull()) throw CasadiException("IdasInternal::init(): No Jacobian has been provided.");

      // Make sure that a linear solver has been providided
      if(linsol_.isNull()) throw CasadiException("IdasInternal::init(): No user defined linear solver has been provided.");

      // Pass to IDA
      flag = IDASpilsSetPreconditioner(mem_, psetup_wrapper, psolve_wrapper);
      if(flag != IDA_SUCCESS) idas_error("IDASpilsSetPreconditioner",flag);
    }
  } else if(getOption("linear_solver")=="user_defined") {
    initUserDefinedLinearSolver();
  } else throw CasadiException("IDAS: Unknown linear solver");

  // Quadrature equations
  if(nq_>0){

    // Allocate n-vectors for quadratures
    yQ_ = N_VNew_Serial(nq_);

    // Initialize quadratures in IDAS
    flag = IDAQuadInit(mem_, rhsQ_wrapper, yQ_);
    if(flag != IDA_SUCCESS) idas_error("IDAQuadInit",flag);
    
    // Should the quadrature errors be used for step size control?
    if(getOption("quad_err_con").toInt()){
      flag = IDASetQuadErrCon(mem_, true);
      if(flag != IDA_SUCCESS) idas_error("IDASetQuadErrCon",flag);
      
      // Quadrature error tolerances
      flag = IDAQuadSStolerances(mem_, reltol_, abstol_); // TODO: vector absolute tolerances
      if(flag != IDA_SUCCESS) idas_error("IDAQuadSStolerances",flag);
    }
  }
  
  log("IdasInternal::init","attached linear solver");
    
 // Sensitivities
 if(ad_order_>0){
   
   // Forward sensitivity problem
   if(nfdir_>0){
     // Allocate n-vectors
     yzS_.resize(nfdir_,0);
     yPS_.resize(nfdir_,0);
     yQS_.resize(nfdir_,0);
     for(int i=0; i<nfdir_; ++i){
        yzS_[i] = N_VNew_Serial(nyz_);
        yPS_[i] = N_VNew_Serial(nyz_);
      }

      // Allocate n-vectors for quadratures
      if(nq_>0){
        for(int i=0; i<nfdir_; ++i){
          yQS_[i] = N_VNew_Serial(nq_);
        }
      }
     
    // Get the sensitivity method
    if(getOption("sensitivity_method")== "simultaneous") ism_ = IDA_SIMULTANEOUS;
    else if(getOption("sensitivity_method")=="staggered") ism_ = IDA_STAGGERED;
    else throw CasadiException("IDAS: Unknown sensitivity method");

    // Copy the forward seeds
    for(int i=0; i<nfdir_; ++i){
      copyNV(fwdSeed(INTEGRATOR_X0,i),fwdSeed(INTEGRATOR_XP0,i),fwdSeed(INTEGRATOR_Z0,i),yzS_[i],yPS_[i],yQS_[i]);
    }

    // Initialize forward sensitivities
    if(finite_difference_fsens_){
     
      // Use finite differences to calculate the residual in the forward sensitivity equations
      flag = IDASensInit(mem_,nfdir_,ism_,0,&yzS_[0],&yPS_[0]);
      if(flag != IDA_SUCCESS) idas_error("IDASensInit",flag);

      // Scaling factors
      double* pbar = 0;
      vector<double> pbar_vec;
      if(hasSetOption("fsens_scaling_factors")){
        pbar_vec = getOption("fsens_scaling_factors").toDoubleVector();
        pbar = &pbar_vec[0];
      }
      
      // Which parameters should be used to estimate the sensitivity equations
      int * plist = 0;
      vector<int> plist_vec;
      if(hasSetOption("fsens_sensitiviy_parameters")){
        plist_vec = getOption("fsens_sensitiviy_parameters").toIntVector();
        plist = &plist_vec[0];
      }

      // Specify parameters
      flag = IDASetSensParams(mem_,&input(INTEGRATOR_P)[0],pbar,plist);
      if(flag != IDA_SUCCESS) idas_error("IDASetSensParams",flag);

      //  IDASetSensDQMethod
      // ?

    } else {
      // Use AD to calculate the residual in the forward sensitivity equations
      flag = IDASensInit(mem_,nfdir_,ism_,resS_wrapper,&yzS_[0],&yPS_[0]);
      if(flag != IDA_SUCCESS) idas_error("IDASensInit",flag);
    }

    // Vector absolute tolerances
    vector<double> fsens_abstol(nfdir_,fsens_abstol_);
    
    // Set tolerances
    if(hasSetOption("fsens_abstolv")){
      // quick hack
      vector<double> fsens_abstolv = getOption("fsens_abstolv").toDoubleVector();
      N_Vector nv_abstol = N_VMake_Serial(fsens_abstolv.size(),&fsens_abstolv[0]);
      vector<N_Vector> nv_abstol_all(nfdir_,nv_abstol);
      flag = IDASensSVtolerances(mem_,fsens_reltol_,&nv_abstol_all[0]);
      if(flag != IDA_SUCCESS) idas_error("IDASensSVtolerances",flag);
      N_VDestroy_Serial(nv_abstol);
    } else {
      flag = IDASensSStolerances(mem_,fsens_reltol_,&fsens_abstol[0]);
      if(flag != IDA_SUCCESS) idas_error("IDASensSStolerances",flag);
    }

    // Set optional inputs
    int errconS = getOption("fsens_err_con").toInt();
    flag = IDASetSensErrCon(mem_, errconS);
    if(flag != IDA_SUCCESS) idas_error("IDASetSensErrCon",flag);

    // Quadrature equations
    if(nq_>0){
      flag = IDAQuadSensInit(mem_, rhsQS_wrapper, &yQS_[0]);
      if(flag != IDA_SUCCESS) idas_error("IDAQuadSensInit",flag);
      
      // Set tolerances
      flag = IDAQuadSensSStolerances(mem_,fsens_reltol_,&fsens_abstol[0]);
      if(flag != IDA_SUCCESS) idas_error("IDAQuadSensSStolerances",flag);
    }
    
    log("IdasInternal::init","initialized forward sensitivities");
  } // enable fsens

  if(nadir_>0){
      // Adjoint sensitivity problem
      whichB_.resize(nadir_);

      // Allocate n-vectors
      yzB_.resize(nadir_);
      yPB_.resize(nadir_);
      for(int i=0; i<nadir_; ++i){
        yzB_[i] = N_VNew_Serial(nyz_);
        yPB_[i] = N_VNew_Serial(nyz_);
      }

      // Allocate n-vectors for the adjoint sensitivities of the parameters
      if(np_>0){
        yBB_.resize(nadir_);
        for(int i=0; i<nadir_; ++i){
          yBB_[i] = N_VMake_Serial(np_,&adjSens(INTEGRATOR_P,i)[0]);
        }
      }

      // Get the number of steos per checkpoint
      int Nd = getOption("steps_per_checkpoint").toInt();

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
  }
  log("IdasInternal::init","initialized adjoint sensitivities");
 } // ad_order>0
 
 is_init = true;
 isInitAdj_ = false;
 log("IdasInternal::init","end");
}

void IdasInternal::initAdj(){
  int flag;
  
  for(int dir=0; dir<nadir_; ++dir){

    // Create backward problem
    flag = IDACreateB(mem_, &whichB_[dir]);
    if(flag != IDA_SUCCESS) idas_error("IDACreateB",flag);
  
    // Initialize the backward problem
    double tB0 = input(INTEGRATOR_TF)[0];
    flag = IDAInitB(mem_, whichB_[dir], resB_wrapper, tB0, yzB_[dir], yPB_[dir]);
    if(flag != IDA_SUCCESS) idas_error("IDAInitB",flag);

    // Set tolerances
    flag = IDASStolerancesB(mem_, whichB_[dir], asens_reltol_, asens_abstol_);
    if(flag!=IDA_SUCCESS) idas_error("IDASStolerancesB",flag);

    // User data
    flag = IDASetUserDataB(mem_, whichB_[dir], this);
    if(flag != IDA_SUCCESS) idas_error("IDASetUserDataB",flag);

    // Maximum number of steps
    IDASetMaxNumStepsB(mem_, whichB_[dir], getOption("max_num_steps").toInt());
    if(flag != IDA_SUCCESS) idas_error("IDASetMaxNumStepsB",flag);
  
    // Pass this information to IDAS
    flag = IDASetIdB(mem_, whichB_[dir], id_);
    if(flag != IDA_SUCCESS) idas_error("IDASetIdB",flag);
      
    // attach linear solver
    if(getOption("asens_linear_solver")=="dense"){
      // Dense jacobian
      flag = IDADenseB(mem_, whichB_[dir], nyz_);
      if(flag != IDA_SUCCESS) idas_error("IDADenseB",flag);
    } else if(getOption("asens_linear_solver")=="banded") {
      // Banded jacobian
      flag = IDABandB(mem_, whichB_[dir], nyz_, getOption("asens_upper_bandwidth").toInt(), getOption("asens_lower_bandwidth").toInt());
      if(flag != IDA_SUCCESS) idas_error("IDABand",flag);
    } else if(getOption("asens_linear_solver")=="iterative") {
      // Sparse solver  
      int maxl = getOption("asens_max_krylov").toInt();
      Option is = getOption("asens_iterative_solver");
      if(is=="gmres"){
        flag = IDASpgmrB(mem_, whichB_[dir], maxl);
        if(flag != IDA_SUCCESS) idas_error("IDASpgmrB",flag);
      } else if(is=="bcgstab"){
        flag = IDASpbcgB(mem_, whichB_[dir], maxl);
        if(flag != IDA_SUCCESS) idas_error("IDASpbcgB",flag);
      } else if(is=="tfqmr") {
        flag = IDASptfqmrB(mem_, whichB_[dir], maxl);
        if(flag != IDA_SUCCESS) idas_error("IDASptfqmrB",flag);
      } else throw CasadiException("Unknown sparse solver for backward problem");
    } else throw CasadiException("Unknown linear solver for backward problem");
      
    // Quadratures for the adjoint problem
    flag = IDAQuadInitB(mem_,whichB_[dir],rhsQB_wrapper,yBB_[dir]);
    if(flag!=IDA_SUCCESS) idas_error("IDAQuadInitB",flag);
  
    // Quadrature error control
    if(getOption("quad_err_con").toInt()){
      flag = IDASetQuadErrConB(mem_, whichB_[dir],true);
      if(flag != IDA_SUCCESS) idas_error("IDASetQuadErrConB",flag);
      
      flag = IDAQuadSStolerancesB(mem_, whichB_[dir], asens_reltol_, asens_abstol_);
      if(flag != IDA_SUCCESS) idas_error("IDAQuadSStolerancesB",flag);
    }
  }
  
  // Mark initialized
  isInitAdj_ = true;
}


void IdasInternal::res(double t, const double* yz, const double* yp, double* r){

  // Get time
  time1 = clock();
  
  // Pass input
  f_.setInput(t,DAE_T);
  f_.setInput(yz,DAE_Y);
  f_.setInput(yp,DAE_YDOT);
  f_.setInput(yz+ny_,DAE_Z);
  f_.setInput(input(INTEGRATOR_P),DAE_P);

  // Evaluate
  f_.evaluate();
  
  // Get results
  f_.getOutput(r);

  // Check the result for consistency
  for(int i=0; i<ny_+nz_; ++i){
    if(isnan(r[i]) || isinf(r[i])){
      if(verbose_){
        stringstream ss;
        ss << "Warning: The " << i << "-th component of the DAE residual is " << r[i] << " at time t=" << t << ".";
        log("IdasInternal::res",ss.str());
        if(monitored("IdasInternal::res")){
          cout << "DAE_T    = " << t << endl;
          cout << "DAE_Y    = " << f_.input(DAE_Y) << endl;
          cout << "DAE_YDOT = " << f_.input(DAE_YDOT) << endl;
          cout << "DAE_Z    = " << f_.input(DAE_Z) << endl;
          cout << "DAE_P    = " << f_.input(DAE_P) << endl;
          cout << "residual = " << f_.output() << endl;
        }
      }
      throw 1;
    }
  }
  
  time2 = clock();
  t_res += double(time2-time1)/CLOCKS_PER_SEC;
}

int IdasInternal::res_wrapper(double t, N_Vector yz, N_Vector yp, N_Vector rr, void *user_data){
 try{
    IdasInternal *this_ = (IdasInternal*)user_data;
    this_->res(t,NV_DATA_S(yz),NV_DATA_S(yp),NV_DATA_S(rr));
    return 0;
  } catch(int flag){ // recoverable error
    return flag;
  } catch(exception& e){ // non-recoverable error
    cerr << "res failed: " << e.what() << endl;;
    return -1;
  }
}

void IdasInternal::ehfun_wrapper(int error_code, const char *module, const char *function, char *msg, void *eh_data){
 try{
    IdasInternal *this_ = (IdasInternal*)eh_data;
    this_->ehfun(error_code,module,function,msg);        
  } catch(exception& e){
    cerr << "ehfun failed: " << e.what() << endl;;
  }
}
  
void IdasInternal::ehfun(int error_code, const char *module, const char *function, char *msg){
  cerr << msg << endl;
}

void IdasInternal::jtimes(double t, const double *yz, const double *yp, const double *rr, const double *v, double *Jv, double cj, double *tmp1, double *tmp2){
  // Get time
  time1 = clock();
  
   // Pass input
   f_.setInput(t,DAE_T);
   f_.setInput(yz,DAE_Y);
   f_.setInput(yp,DAE_YDOT);
   f_.setInput(yz+ny_,DAE_Z);
   f_.setInput(input(INTEGRATOR_P),DAE_P);
     
   // Pass seeds of the state vectors
   f_.setFwdSeed(v,DAE_Y);
   f_.setFwdSeed(v+ny_,DAE_Z);
   
   // Pass seeds of the state derivative
   for(int i=0; i<ny_; ++i) tmp1[i] = cj*v[i];
   f_.setFwdSeed(tmp1,DAE_YDOT);
   
   // Evaluate the AD forward algorithm
   f_.evaluate(1,0);
   
   // Get the output seeds
   f_.getFwdSens(Jv);

  // Log time duration
  time2 = clock();
  t_jac += double(time2-time1)/CLOCKS_PER_SEC;
}

int IdasInternal::jtimes_wrapper(double t, N_Vector yz, N_Vector yp, N_Vector rr, N_Vector v, N_Vector Jv, double cj, void *user_data, N_Vector tmp1, N_Vector tmp2){
 try{
    IdasInternal *this_ = (IdasInternal*)user_data;
    this_->jtimes(t,NV_DATA_S(yz),NV_DATA_S(yp),NV_DATA_S(rr),NV_DATA_S(v),NV_DATA_S(Jv),cj,NV_DATA_S(tmp1),NV_DATA_S(tmp2));
    return 0;
  } catch(exception& e){
    cerr << "jtimes failed: " << e.what() << endl;;
    return 1;
  }  
}

void IdasInternal::resS(int Ns, double t, const double* yz, const double* yp, const double *resval, N_Vector *yS, N_Vector* ypS, N_Vector *resvalS, double *tmp1, double *tmp2, double *tmp3){
  assert(Ns==nfdir_);

  // Record the current cpu time
  time1 = clock();
  
   // Pass input
   f_.setInput(t,DAE_T);
   f_.setInput(yz,DAE_Y);
   f_.setInput(yp,DAE_YDOT);
   f_.setInput(yz+ny_,DAE_Z);
   f_.setInput(input(INTEGRATOR_P),DAE_P);

   // Calculate the forward sensitivities, nfdir_f_ directions at a time
   for(int j=0; j<nfdir_; j += nfdir_f_){
     for(int dir=0; dir<nfdir_f_ && j+dir<nfdir_; ++dir){
       // Pass forward seeds 
       f_.setFwdSeed(0.0,DAE_T,dir);
       f_.setFwdSeed(NV_DATA_S(yS[j+dir]),DAE_Y,dir);
       f_.setFwdSeed(NV_DATA_S(ypS[j+dir]),DAE_YDOT,dir);
       f_.setFwdSeed(NV_DATA_S(yS[j+dir])+ny_,DAE_Z,dir);
       f_.setFwdSeed(fwdSeed(INTEGRATOR_P,j+dir),DAE_P,dir);
     }
   
     // Evaluate the AD forward algorithm
     f_.evaluate(1,0);
      
     // Get the output seeds
     for(int dir=0; dir<nfdir_f_ && j+dir<nfdir_; ++dir){
       f_.getFwdSens(NV_DATA_S(resvalS[j+dir]),DAE_RES,dir);
     }
   }
   
   // Record timings
   time2 = clock();
   t_fres += double(time2-time1)/CLOCKS_PER_SEC;
}

int IdasInternal::resS_wrapper(int Ns, double t, N_Vector yz, N_Vector yp, N_Vector resval, N_Vector *yS, N_Vector *ypS, N_Vector *resvalS, void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3){
 try{
    IdasInternal *this_ = (IdasInternal*)user_data;
    this_->resS(Ns,t,NV_DATA_S(yz),NV_DATA_S(yp),NV_DATA_S(resval),yS,ypS,resvalS,NV_DATA_S(tmp1),NV_DATA_S(tmp2),NV_DATA_S(tmp3));
    return 0;
  } catch(exception& e){
    cerr << "resS failed: " << e.what() << endl;;
    return 1;
  }
}

void IdasInternal::reset(int fsens_order, int asens_order){
  log("IdasInternal::reset","begin");
  
  // If we have forward sensitivities, rest one extra time without forward sensitivities to get a consistent initial guess
  if(fsens_order>0 && getOption("extra_fsens_calc_ic").toInt())
    reset(0,asens_order);

  // Reset timers
  t_res = t_fres = t_jac = t_lsolve = t_lsetup_jac = t_lsetup_fac = 0;
    
  fsens_order_ = fsens_order;
  asens_order_ = asens_order;
  
  // Get the time horizon
  double t0 = input(INTEGRATOR_T0)[0];
  double tf = input(INTEGRATOR_TF)[0];
  t_ = t0;
  
  // Return flag
  int flag;
  
  // Copy to N_Vectors
  copyNV(input(INTEGRATOR_X0),input(INTEGRATOR_XP0),input(INTEGRATOR_Z0),yz_,yP_,yQ_);
  
  // Re-initialize
  flag = IDAReInit(mem_, t0, yz_, yP_);
  if(flag != IDA_SUCCESS) idas_error("IDAReInit",flag);

  // Re-initialize quadratures
  if(nq_>0){
    flag = IDAQuadReInit(mem_, yQ_);
    if(flag != IDA_SUCCESS) idas_error("IDAQuadReInit",flag);
  }
  
  if(fsens_order_>0){
    // Get the forward seeds
    for(int i=0; i<nfdir_; ++i){
      copyNV(fwdSeed(INTEGRATOR_X0,i),fwdSeed(INTEGRATOR_XP0,i),fwdSeed(INTEGRATOR_Z0,i),yzS_[i],yPS_[i],yQS_[i]);
    }
    
    // Re-initialize sensitivities
    flag = IDASensReInit(mem_,ism_,&yzS_[0],&yPS_[0]);
    if(flag != IDA_SUCCESS) idas_error("IDASensReInit",flag);

    if(nq_>0){
      flag = IDAQuadSensReInit(mem_, &yQS_[0]);
      if(flag != IDA_SUCCESS) idas_error("IDAQuadSensReInit",flag);
    }
  } else {
    // Turn of sensitivities
    flag = IDASensToggleOff(mem_);
    if(flag != IDA_SUCCESS) idas_error("IDASensToggleOff",flag);
  }

  // Correct initial conditions, if necessary
  int calc_ic = getOption("calc_ic").toInt();
  if(calc_ic){
    correctInitialConditions();
  }
    
  log("IdasInternal::reset","end");
}
  
  
void IdasInternal::correctInitialConditions(){
  double t0 = input(INTEGRATOR_T0)[0];
  double tf = input(INTEGRATOR_TF)[0];
  log("IdasInternal::correctInitialConditions","begin");
  if(monitored("IdasInternal::correctInitialConditions")){
    cout << "initial guess: " << endl;
    cout << "p = " << input(INTEGRATOR_P) << endl;
    cout << "x0 = " << input(INTEGRATOR_X0) << endl;
    cout << "xp0 = " << input(INTEGRATOR_XP0) << endl;
    cout << "z0 = " << input(INTEGRATOR_Z0) << endl;
    if(fsens_order_>0){
      for(int dir=0; dir<nfdir_; ++dir){
        cout << "forward seed guess, direction " << dir << ": " << endl;
        cout << "p_seed = " << fwdSeed(INTEGRATOR_P,dir) << endl;
        cout << "x0_seed = " << fwdSeed(INTEGRATOR_X0,dir) << endl;
        cout << "xp0_seed = " << fwdSeed(INTEGRATOR_XP0,dir) << endl;
        cout << "z0_seed = " << fwdSeed(INTEGRATOR_Z0,dir) << endl;
      }
    }
  }

  int icopt = IDA_YA_YDP_INIT; // calculate z and xdot given x
  // int icopt = IDA_Y_INIT; // calculate z and x given zdot and xdot (e.g. start in stationary)

  double t_first = hasSetOption("first_time") ? getOption("first_time").toDouble() : tf;
  int flag = IDACalcIC(mem_, icopt , t_first);
  if(flag != IDA_SUCCESS) idas_error("IDACalcIC",flag);

  // Retrieve the initial values
  flag = IDAGetConsistentIC(mem_, yz_, yP_);
  if(flag != IDA_SUCCESS) idas_error("IDAGetConsistentIC",flag);

  // Save the corrected input to the function arguments
  copyNV(yz_,yP_,yQ_,input(INTEGRATOR_X0),input(INTEGRATOR_XP0),input(INTEGRATOR_Z0));
  if(fsens_order_>0){
    for(int i=0; i<nfdir_; ++i){
      copyNV(yzS_[i],yPS_[i],yQS_[i],fwdSeed(INTEGRATOR_X0,i),fwdSeed(INTEGRATOR_XP0,i),fwdSeed(INTEGRATOR_Z0,i));
    }
  }
  
  // Print progress
  log("IdasInternal::correctInitialConditions","found consistent initial values");
  if(monitored("IdasInternal::correctInitialConditions")){
    cout << "p = " << input(INTEGRATOR_P) << endl;
    cout << "x0 = " << input(INTEGRATOR_X0) << endl;
    cout << "xp0 = " << input(INTEGRATOR_XP0) << endl;
    cout << "z0 = " << input(INTEGRATOR_Z0) << endl;
    if(fsens_order_>0){
      for(int dir=0; dir<nfdir_; ++dir){
        cout << "forward seed, direction " << dir << ": " << endl;
        cout << "p_seed = " << fwdSeed(INTEGRATOR_P,dir) << endl;
        cout << "x0_seed = " << fwdSeed(INTEGRATOR_X0,dir) << endl;
        cout << "xp0_seed = " << fwdSeed(INTEGRATOR_XP0,dir) << endl;
        cout << "z0_seed = " << fwdSeed(INTEGRATOR_Z0,dir) << endl;
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
      flag = IDASolveF(mem_, t_out, &t_, yz_, yP_, IDA_NORMAL, &ncheck_);
      if(flag != IDA_SUCCESS && flag != IDA_TSTOP_RETURN) idas_error("IDASolveF",flag);
    } else {
      // ... without taping
      log("IdasInternal::integrate","integration without taping");
      flag = IDASolve(mem_, t_out, &t_, yz_, yP_, IDA_NORMAL);
      if(flag != IDA_SUCCESS && flag != IDA_TSTOP_RETURN) idas_error("IDASolve",flag);
    }
    log("IdasInternal::integrate","integration complete");
    
    // Get quadrature states
    if(nq_>0){
      double tret;
      flag = IDAGetQuad(mem_, &tret, yQ_);
      if(flag != IDA_SUCCESS) idas_error("IDAGetQuad",flag);
    }
    
    if(fsens_order_>0){
      // Get the sensitivities
      flag = IDAGetSens(mem_,&t_, &yzS_[0]);
      if(flag != IDA_SUCCESS) idas_error("IDAGetSens",flag);
    
      if(nq_>0){
        double tret;
        flag = IDAGetQuadSens(mem_, &tret, &yQS_[0]);
        if(flag != IDA_SUCCESS) idas_error("IDAGetQuadSens",flag);
      }
    }
  }
  
  // Save the final state
  copyNV(yz_,yP_,yQ_,output(INTEGRATOR_XF),output(INTEGRATOR_XPF),output(INTEGRATOR_ZF));
  if(fsens_order_>0){
    for(int i=0; i<nfdir_; ++i){
      copyNV(yzS_[i],yPS_[i],yQS_[i],fwdSens(INTEGRATOR_XF,i),fwdSens(INTEGRATOR_XPF,i),fwdSens(INTEGRATOR_ZF,i));
    }
  }
  log("IdasInternal::integrate","end");
}

void IdasInternal::resetAdj(){
  double t0 = input(INTEGRATOR_T0)[0];
  double tf = input(INTEGRATOR_TF)[0];
  
  int flag;
  // Reset adjoint sensitivities for the parameters
  if(np_>0)
    for(int dir=0; dir<nadir_; ++dir)
      N_VConst(0.0, yBB_[dir]);
  
  // Get the adjoint seeds
  getAdjointSeeds();
    
  if(isInitAdj_){
    for(int dir=0; dir<nadir_; ++dir){
      flag = IDAReInitB(mem_, whichB_[dir], tf, yzB_[dir], yPB_[dir]);
      if(flag != IDA_SUCCESS) idas_error("IDAReInitB",flag);
      
      flag = IDAQuadReInit(IDAGetAdjIDABmem(mem_, whichB_[dir]),yBB_[dir]);
//      flag = IDAQuadReInitB(mem_,whichB_[dir],yBB_[dir]); // BUG in Sundials - do not use this!
      if(flag!=IDA_SUCCESS) idas_error("IDAQuadReInitB",flag);
    }
  } else {
    // Initialize the adjoint integration
    initAdj();
  }

  // Correct initial values for the integration if necessary
  int calc_icB = getOption("calc_icB").toInt();
  if(calc_icB){
    for(int dir=0; dir<nadir_; ++dir){
      flag = IDACalcICB(mem_, whichB_[dir], t0, yzB_[dir], yPB_[dir]);
      if(flag != IDA_SUCCESS) idas_error("IDACalcICB",flag);

      N_VPrint_Serial(yzB_[dir]);
      N_VPrint_Serial(yPB_[dir]);

      // Retrieve the initial values
      flag = IDAGetConsistentICB(mem_, whichB_[dir], yzB_[dir], yPB_[dir]);
      if(flag != IDA_SUCCESS) idas_error("IDAGetConsistentICB",flag);
      
      N_VPrint_Serial(yzB_[dir]);
      N_VPrint_Serial(yPB_[dir]);
    }
  }
}

void IdasInternal::integrateAdj(double t_out){
  int flag;
  // Integrate backwards to t_out
  flag = IDASolveB(mem_, t_out, IDA_NORMAL);
  if(flag<IDA_SUCCESS) idas_error("IDASolveB",flag);

  // Get the sensitivities
  double tret;
  for(int dir=0; dir<nadir_; ++dir){
    flag = IDAGetB(mem_, whichB_[dir], &tret, yzB_[dir], yPB_[dir]);
    if(flag!=IDA_SUCCESS) idas_error("IDAGetB",flag);

    flag = IDAGetQuadB(mem_, whichB_[dir], &tret, yBB_[dir]);
    if(flag!=IDA_SUCCESS) idas_error("IDAGetQuadB",flag);
  }
  
  // Save the adjoint sensitivities
  setAdjointSensitivities();
}

void IdasInternal::printStats(std::ostream &stream) const{
    long nsteps, nfevals, nlinsetups, netfails;
    int qlast, qcur;
    double hinused, hlast, hcur, tcur;

    int flag = IDAGetIntegratorStats(mem_, &nsteps, &nfevals, &nlinsetups,&netfails, &qlast, &qcur, &hinused,&hlast, &hcur, &tcur);
    if(flag!=IDA_SUCCESS) idas_error("IDAGetIntegratorStats",flag);
    
    stream << "number of steps taken by IDAS: " << nsteps << std::endl;
    stream << "number of calls to the user's f function: " << nfevals << std::endl;
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
  ss << " Consult Idas documentation.";
  throw CasadiException(ss.str());
}

int IdasInternal::rhsQ_wrapper(double t, N_Vector yz, N_Vector yp, N_Vector rhsQ, void *user_data){
 try{
    IdasInternal *this_ = (IdasInternal*)user_data;
    this_->rhsQ(t,NV_DATA_S(yz),NV_DATA_S(yp),NV_DATA_S(rhsQ));
    return 0;
  } catch(exception& e){
    cerr << "rhsQ failed: " << e.what() << endl;;
    return 1;
  }
}

void IdasInternal::rhsQ(double t, const double* yz, const double* yp, double* rhsQ){
   // Pass input
   q_.setInput(t,DAE_T);
   q_.setInput(yz,DAE_Y);
   q_.setInput(yp,DAE_YDOT);
   q_.setInput(yz+ny_,DAE_Z);
   q_.setInput(input(INTEGRATOR_P),DAE_P);

    // Evaluate
   q_.evaluate();
    
    // Get results
   q_.getOutput(rhsQ);
}
  
void IdasInternal::rhsQS(int Ns, double t, N_Vector yz, N_Vector yp, N_Vector *yzS, N_Vector *ypS, N_Vector rrQ, N_Vector *rhsvalQS, 
                        N_Vector tmp1, N_Vector tmp2, N_Vector tmp3){
  assert(Ns==nfdir_);

  // Pass input
   q_.setInput(t,DAE_T);
   q_.setInput(NV_DATA_S(yz),DAE_Y);
   q_.setInput(NV_DATA_S(yp),DAE_YDOT);
   q_.setInput(NV_DATA_S(yz)+ny_,DAE_Z);
   q_.setInput(input(INTEGRATOR_P),DAE_P);
     
   // Pass forward seeds
  for(int i=0; i<nfdir_; ++i){
    q_.setFwdSeed(0.0,DAE_T);
    q_.setFwdSeed(NV_DATA_S(yzS[i]),DAE_Y);
    q_.setFwdSeed(NV_DATA_S(ypS[i]),DAE_YDOT);
    q_.setFwdSeed(NV_DATA_S(yzS[i])+ny_,DAE_Z);
    q_.setFwdSeed(fwdSeed(INTEGRATOR_P,i),DAE_P);
   
    // Evaluate the AD forward algorithm
    q_.evaluate(1,0);
      
    // Get the output seeds
    q_.getFwdSens(NV_DATA_S(rhsvalQS[i]));
  }
}

int IdasInternal::rhsQS_wrapper(int Ns, double t, N_Vector yz, N_Vector yp, N_Vector *yzS, N_Vector *ypS, N_Vector rrQ, N_Vector *rhsvalQS, void *user_data, 
                                N_Vector tmp1, N_Vector tmp2, N_Vector tmp3){
  
 try{
    IdasInternal *this_ = (IdasInternal*)user_data;
    this_->rhsQS(Ns,t,yz,yp,yzS,ypS,rrQ,rhsvalQS,tmp1,tmp2,tmp3);
    return 0;
  } catch(exception& e){
    cerr << "rhsQS failed: " << e.what() << endl;;
    return 1;
  }
}

void IdasInternal::resB(double t, const double* yz, const double* yp, const double* yB, const double* ypB, double* resvalB){
  // Pass input
  f_.setInput(t,DAE_T);
  f_.setInput(yz,DAE_Y);
  f_.setInput(yp,DAE_YDOT);
  f_.setInput(yz+ny_,DAE_Z);
  f_.setInput(input(INTEGRATOR_P),DAE_P);
  
  // Pass adjoint seeds
  f_.setAdjSeed(yB,DAE_RES);

  // Evaluate
  f_.evaluate(0,1);

  // Save to output
  f_.getAdjSens(resvalB,DAE_Y);
  f_.getAdjSens(resvalB+ny_,DAE_Z);

  // Pass adjoint seeds
  f_.setAdjSeed(ypB,DAE_RES);
  
  // Evaluate AD adjoint
  f_.evaluate(0,1);

  // Save to output
  const vector<double>& asens_ydot = f_.adjSens(DAE_YDOT);
  for(int i=0; i<ny_; ++i)
    resvalB[i] -= asens_ydot[i];
  
  // If quadratures are included
  if(nq_>0){
    // Pass input to quadratures
    q_.setInput(t,DAE_T);
    q_.setInput(yz,DAE_Y);
    q_.setInput(yp,DAE_YDOT);
    q_.setInput(yz+ny_,DAE_Z);
    q_.setInput(input(INTEGRATOR_P),DAE_P);

    // Pass adjoint seeds
    q_.setAdjSeed(&output(INTEGRATOR_XF)[ny_],DAE_RES);

    // Evaluate
    q_.evaluate(0,1);
    
    // Get the input seeds
    const vector<double>& asens_y = q_.adjSens(DAE_Y);
    for(int i=0; i<ny_; ++i)
      resvalB[i] += asens_y[i];

    const vector<double>& asens_z = q_.adjSens(DAE_Z);
    for(int i=0; i<nz_; ++i)
      resvalB[i+ny_] += asens_z[i];
  }
}

int IdasInternal::resB_wrapper(double t, N_Vector y, N_Vector yp, N_Vector yB, N_Vector ypB, N_Vector resvalB, void *user_dataB){
 try{
    IdasInternal *this_ = (IdasInternal*)user_dataB;
    this_->resB(t,NV_DATA_S(y),NV_DATA_S(yp),NV_DATA_S(yB),NV_DATA_S(ypB),NV_DATA_S(resvalB));
    return 0;
  } catch(exception& e){
    cerr << "resB failed: " << e.what() << endl;;
    return 1;
  }
}
  
void IdasInternal::rhsQB(double t, const double* yz, const double* yp, const double* yB, const double* ypB, double *rhsvalBQ){
  // Pass input
  f_.setInput(t,DAE_T);
  f_.setInput(yz,DAE_Y);
  f_.setInput(yp,DAE_YDOT);
  f_.setInput(yz+ny_,DAE_Z);
  f_.setInput(input(INTEGRATOR_P),DAE_P);

  // Pass adjoint seeds
  f_.setAdjSeed(yB,DAE_RES);

  // Evaluate
  f_.evaluate(0,1);

  // Save to output
  f_.getAdjSens(rhsvalBQ,DAE_P);
  
  // If quadratures are included
  if(nq_>0){
    // Pass input to quadratures
    q_.setInput(t,DAE_T);
    q_.setInput(yz,DAE_Y);
    q_.setInput(yp,DAE_YDOT);
    q_.setInput(yz+ny_,DAE_Z);
    q_.setInput(input(INTEGRATOR_P),DAE_P);

    // Pass adjoint seeds
    q_.setAdjSeed(&output(INTEGRATOR_XF)[ny_],DAE_RES);

    // Evaluate
    q_.evaluate(0,1);
    
    // Get the input seeds
    const vector<double>& qres = q_.adjSens(DAE_P);
    
    // Copy to result
    for(int i=0; i<np_; ++i){
      rhsvalBQ[i] -= qres[i];
    }
  }
}

int IdasInternal::rhsQB_wrapper(double t, N_Vector y, N_Vector yp, N_Vector yB, N_Vector ypB, N_Vector rhsvalBQ, void *user_dataB){
 try{
    IdasInternal *this_ = (IdasInternal*)user_dataB;
    this_->rhsQB(t,NV_DATA_S(y),NV_DATA_S(yp),NV_DATA_S(yB),NV_DATA_S(ypB),NV_DATA_S(rhsvalBQ));
    return 0;
  } catch(exception& e){
    cerr << "resQB failed: " << e.what() << endl;;
    return 1;
  }
}

void IdasInternal::djac(int Neq, double t, double cj, N_Vector yz, N_Vector yp, N_Vector rr, DlsMat Jac, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3){
  // NOTE: This function is extra complicated due to the fact that we either have a function that calculates df_dx and df_dxdot together or separately
  // TODO: Change this when MXFunction becomes more stable
  
  // Get time
  time1 = clock();

  // Number of jacobian functions
  int njac = jac_.isNull() ? 2 : 1;
  assert(njac!=2);
  
  // Evaluate the functions and add to the result
  for(int ijac = 0; ijac<njac; ++ijac){
    // Jacobian function reference (jac_, jacx_ or jacxdot_)
    FX& jac = njac==1 ? jac_ : ijac==0 ? jacx_ : jacxdot_;
    
    // Pass input to the jacobian function
    if(njac==1){
      // if a function calculating df_dx + cj*df_dxdot has been provided
      jac.setInput(t,JAC_T);
      jac.setInput(NV_DATA_S(yz),JAC_Y);
      jac.setInput(NV_DATA_S(yp),JAC_YDOT);
      jac.setInput(NV_DATA_S(yz)+ny_,JAC_Z);
      jac.setInput(input(INTEGRATOR_P),JAC_P);
      jac.setInput(cj,JAC_CJ);
    } else {
      // if we need to calculate df_dx and df_dxdot separately
      jac.setInput(t,DAE_T);
      jac.setInput(NV_DATA_S(yz),DAE_Y);
      jac.setInput(NV_DATA_S(yp),DAE_YDOT);
      jac.setInput(NV_DATA_S(yz)+ny_,DAE_Z);
      jac.setInput(input(INTEGRATOR_P),DAE_P);
    }
    
    // Evaluate jacobian
    jac.evaluate();

    // Get sparsity and non-zero elements
    const vector<int>& rowind = jac.output().rowind();
    const vector<int>& col = jac.output().col();
    const vector<double>& val = jac.output();

    // Factor
    double c = ijac==0 ? 1 : cj;
    
    // Dimension of the jacobian
    int jdim = jac.output().size1();
    
    // Loop over rows
    for(int offset=0; offset<ny_; offset += jdim){
      for(int i=0; i<rowind.size()-1; ++i){
        // Loop over non-zero entries
        for(int el=rowind[i]; el<rowind[i+1]; ++el){
          // Get column
          int j = col[el];
      
          // Add to the element
          DENSE_ELEM(Jac,i+offset,j+offset) += c*val[el];
        }
      }
    }
  }
  
  // Log time duration
  time2 = clock();
  t_jac += double(time2-time1)/CLOCKS_PER_SEC;
}

int IdasInternal::djac_wrapper(int Neq, double t, double cj, N_Vector yz, N_Vector yp, N_Vector rr, DlsMat Jac, void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3){
 try{
    IdasInternal *this_ = (IdasInternal*)user_data;
    this_->djac(Neq, t, cj, yz, yp, rr, Jac, tmp1, tmp2, tmp3);
    return 0;
  } catch(exception& e){
    cerr << "djac failed: " << e.what() << endl;;
    return 1;
  }
}

void IdasInternal::bjac(int Neq, int mupper, int mlower, double tt, double cj, N_Vector yz, N_Vector yp, N_Vector rr, DlsMat Jac, N_Vector tmp1, N_Vector tmp2,N_Vector tmp3){
  // Get time
  time1 = clock();

  // Pass input to the jacobian function
  jac_.setInput(tt,JAC_T);
  jac_.setInput(NV_DATA_S(yz),JAC_Y);
  jac_.setInput(NV_DATA_S(yp),JAC_YDOT);
  jac_.setInput(NV_DATA_S(yz)+ny_,JAC_Z);
  jac_.setInput(input(INTEGRATOR_P),JAC_P);
  jac_.setInput(cj,JAC_CJ);

  // Evaluate jacobian
  jac_.evaluate();

  // Get sparsity and non-zero elements
  const vector<int>& rowind = jac_.output().rowind();
  const vector<int>& col = jac_.output().col();
  const vector<double>& val = jac_.output();

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
}

int IdasInternal::bjac_wrapper(int Neq, int mupper, int mlower, double tt, double cj, N_Vector yz, N_Vector yp, N_Vector rr, DlsMat Jac, void *user_data, N_Vector tmp1, N_Vector tmp2,N_Vector tmp3){
 try{
    IdasInternal *this_ = (IdasInternal*)user_data;
    this_->bjac(Neq, mupper, mlower, tt, cj, yz, yp, rr, Jac, tmp1, tmp2, tmp3);
    return 0;
  } catch(exception& e){
    cerr << "bjac failed: " << e.what() << endl;;
    return 1;
  }
}

void IdasInternal::setStopTime(double tf){
  // Set the stop time of the integration -- don't integrate past this point
  int flag = IDASetStopTime(mem_, tf);
  if(flag != IDA_SUCCESS) idas_error("IDASetStopTime",flag);
}

int IdasInternal::psolve_wrapper(double t, N_Vector yz, N_Vector yp, N_Vector rr, N_Vector rvec, N_Vector zvec, double cj, double delta, void *user_data, N_Vector tmp){
 try{
    IdasInternal *this_ = (IdasInternal*)user_data;
    assert(this_);
    this_->psolve(t, yz, yp, rr, rvec, zvec, cj, delta, tmp);
    return 0;
  } catch(exception& e){
    cerr << "psolve failed: " << e.what() << endl;;
    return 1;
  }
}

int IdasInternal::psetup_wrapper(double t, N_Vector yz, N_Vector yp, N_Vector rr, double cj, void* user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3){
 try{
    IdasInternal *this_ = (IdasInternal*)user_data;
    assert(this_);
    this_->psetup(t, yz, yp, rr, cj, tmp1, tmp2, tmp3);
    return 0;
  } catch(exception& e){
    cerr << "psetup failed: " << e.what() << endl;;
    return 1;
  }
}

void IdasInternal::psolve(double t, N_Vector yz, N_Vector yp, N_Vector rr, N_Vector rvec, N_Vector zvec, double cj, double delta, N_Vector tmp){
  // Get time
  time1 = clock();

  // Number of right-hand-sides
  int nrhs = nrhs_;
  
  // Number of rows of the linear solver
  int nrow = nyz_/nrhs;
  
  // Pass right hand side to the linear solver (transpose necessary)
  linsol_.setInput(NV_DATA_S(rvec),1);
  
  // Solve the (possibly factorized) system 
  linsol_.solve();
  
  // Get the result (transpose necessary)
  linsol_.getOutput(NV_DATA_S(zvec));

  // Log time duration
  time2 = clock();
  t_lsolve += double(time2-time1)/CLOCKS_PER_SEC;

}

void IdasInternal::psetup(double t, N_Vector yz, N_Vector yp, N_Vector rr, double cj, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3){
  // Get time
  time1 = clock();

  // Pass input to the jacobian function
  jac_.setInput(t,JAC_T);
  jac_.setInput(NV_DATA_S(yz),JAC_Y);
  jac_.setInput(NV_DATA_S(yp),JAC_YDOT);
  jac_.setInput(NV_DATA_S(yz)+ny_,JAC_Z);
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

}

int IdasInternal::lsetup_wrapper(IDAMem IDA_mem, N_Vector yzp, N_Vector ypp, N_Vector resp, N_Vector vtemp1, N_Vector vtemp2, N_Vector vtemp3){
 try{
    IdasInternal *this_ = (IdasInternal*)(IDA_mem->ida_lmem);
    assert(this_);
    this_->lsetup(IDA_mem,yzp,ypp,resp,vtemp1,vtemp2,vtemp3);
    return 0;
  } catch(exception& e){
    cerr << "lsetup failed: " << e.what() << endl;;
    return -1;
  }
}

int IdasInternal::lsolve_wrapper(IDAMem IDA_mem, N_Vector b, N_Vector weight, N_Vector ycur, N_Vector ypcur, N_Vector rescur){
 try{
   IdasInternal *this_ = (IdasInternal*)(IDA_mem->ida_lmem);
   assert(this_);
   this_->lsolve(IDA_mem,b,weight,ycur,ypcur,rescur);
   return 0;
  } catch(int wrn){
    return wrn;
  } catch(exception& e){
    cerr << "lsolve failed: " << e.what() << endl;;
    return -1;
  }
}

void IdasInternal::lsetup(IDAMem IDA_mem, N_Vector yzp, N_Vector ypp, N_Vector resp, N_Vector vtemp1, N_Vector vtemp2, N_Vector vtemp3){  
  // Current time
  double t = IDA_mem->ida_tn;

  // Multiple of df_dydot to be added to the matrix
  double cj = IDA_mem->ida_cj;

  // Call the preconditioner setup function (which sets up the linear solver)
  psetup(t, yzp, ypp, 0, cj, vtemp1, vtemp1, vtemp3);
  
}

void IdasInternal::lsolve(IDAMem IDA_mem, N_Vector b, N_Vector weight, N_Vector ycur, N_Vector ypcur, N_Vector rescur){
  
  // Current time
  double t = IDA_mem->ida_tn;

  // Multiple of df_dydot to be added to the matrix
  double cj = IDA_mem->ida_cj;

  // Accuracy
  double delta = 0.0;
  
  // Call the preconditioner solve function (which solves the linear system)
  psolve(t, ycur, ypcur, rescur, b, b, cj, delta, 0);
  
  // Scale the correction to account for change in cj
  if(cj_scaling_){
    double cjratio = IDA_mem->ida_cjratio;
    if (cjratio != 1.0) N_VScale(2.0/(1.0 + cjratio), b, b);
  }
}

void IdasInternal::initUserDefinedLinearSolver(){
  // Make sure that a Jacobian has been provided
  if(jac_.isNull()) throw CasadiException("IdasInternal::initUserDefinedLinearSolver(): No Jacobian has been provided.");

  // Make sure that a linear solver has been providided
  if(linsol_.isNull()) throw CasadiException("IdasInternal::initUserDefinedLinearSolver(): No user defined linear solver has been provided.");

  //  Set fields in the IDA memory
  IDAMem IDA_mem = IDAMem(mem_);
  IDA_mem->ida_lmem   = this;
  IDA_mem->ida_lsetup = lsetup_wrapper;
  IDA_mem->ida_lsolve = lsolve_wrapper;
  IDA_mem->ida_setupNonNull = TRUE;
}

void IdasInternal::setLinearSolver(const LinearSolver& linsol, const FX& jac){
  linsol_ = linsol;
  jac_ = jac;

  // Try to generate a jacobian of none provided
  if(jac_.isNull())
    getJacobian();
}

Integrator IdasInternal::jac(int iind, int oind){
  // Make sure that the derivative is of the state
  if(oind!=INTEGRATOR_XF)
    throw CasadiException("IdasInternal::jacobian: Not derivative of state");
  
  // Type cast to SXMatrix
  SXFunction f = shared_cast<SXFunction>(f_);
  SXFunction q = shared_cast<SXFunction>(q_);

  // Assert that the function is made up by SXFunctions
  if((f_.isNull() != f.isNull()) || (q_.isNull() != q.isNull()))
    throw CasadiException("IdasInternal::jacobian: Only supported for f and q SXFunction");
    
  // Generate Jacobians with respect to state, state derivative and parameters
  SXMatrix df_dy = f.jac(DAE_Y,DAE_RES);
  SXMatrix df_dydot = f.jac(DAE_YDOT,DAE_RES);
  SXMatrix df_dz = f.jac(DAE_Z,DAE_RES);
  
  // Number of sensitivities
  int ns = iind==INTEGRATOR_P ? np_ : ny_;
  
  // Sensitivities and derivatives of sensitivities
  SXMatrix ysens = symbolic("ysens",ny_,ns);
  SXMatrix ypsens = symbolic("ypsens",ny_,ns);
  SXMatrix zsens = symbolic("zsens",nz_,ns);
  
  // Sensitivity equation
  SXMatrix res_s = prod(df_dy,ysens) + prod(df_dydot,ypsens);
  if(nz_>0)
    res_s += prod(df_dz,zsens);
  
  if(iind==INTEGRATOR_P)
    res_s += f.jac(DAE_P,DAE_RES);

  // Augmented DAE
  SXMatrix faug = vec(horzcat(f->outputv.at(oind),res_s));

  // Input to the augmented DAE (start with old)
  vector<SXMatrix> faug_in = f->inputv;

  // Augment differential state
  faug_in[DAE_Y] = vec(horzcat(faug_in[DAE_Y],ysens));

  // Augment state derivatives
  faug_in[DAE_YDOT] = vec(horzcat(faug_in[DAE_YDOT],ypsens));

  // Augment algebraic state
  faug_in[DAE_Z] = vec(horzcat(faug_in[DAE_Z],zsens));
  
  // Create augmented DAE function
  SXFunction ffcn_aug(faug_in,faug);
  ffcn_aug.setOption("ad_order",1);
  
  // Augmented quadratures
  SXFunction qfcn_aug;
  
  if(!q.isNull()){
    // Now lets do the same for the quadrature states
    SXMatrix dq_dy = q.jac(DAE_Y,DAE_RES);
    SXMatrix dq_dydot = q.jac(DAE_YDOT,DAE_RES);
    SXMatrix dq_dz = q.jac(DAE_Z,DAE_RES);
    
    // Sensitivity quadratures
    SXMatrix q_s = prod(dq_dy,ysens) + prod(dq_dydot,ypsens);
    if(nz_>0)
      q_s += prod(dq_dz,zsens);
      
    if(iind==INTEGRATOR_P)
      q_s += q.jac(DAE_P,DAE_RES);

    // Augmented quadratures
    SXMatrix qaug = vec(horzcat(q->outputv.at(oind),q_s));

    // Input to the augmented DAE (start with old)
    vector<SXMatrix> qaug_in = q->inputv;

    // Augmented differential state
    qaug_in[DAE_Y] = vec(horzcat(qaug_in[DAE_Y],ysens));

    // Augment state derivatives
    qaug_in[DAE_YDOT] = vec(horzcat(qaug_in[DAE_YDOT],ypsens));

    // Augmented algebraic state
    qaug_in[DAE_Z] = vec(horzcat(qaug_in[DAE_Z],zsens));

    // Create augmented DAE function
    qfcn_aug = SXFunction(qaug_in,qaug);
    qfcn_aug.setOption("ad_order",1);
  }
  
  // Create integrator instance
  IdasIntegrator integrator(ffcn_aug,qfcn_aug);

  // Set options
  integrator.copyOptions(shared_from_this<IdasIntegrator>());
  integrator.setOption("nrhs",1+ns);
  
  // Transmit information on derivative states
  if(hasSetOption("is_differential")){
    vector<int> is_diff = getOption("is_differential").toIntVector();
    if(is_diff.size()!=nyz_) throw CasadiException("is_differential has incorrect length");
    vector<int> is_diff_aug(nyz_*(1+ns));
    for(int i=0; i<1+ns; ++i)
      for(int j=0; j<nyz_; ++j)
        is_diff_aug[j+i*nyz_] = is_diff[j];
    integrator.setOption("is_differential",is_diff_aug);
  }
  
  // Mapping between states in the augmented dae and the original dae
  vector<int> jacmap(nx_*(1+ns));
  for(int i=0; i<1+ns; ++i){
    for(int j=0; j<nyz_; ++j)
      jacmap[j+nx_*i] = j+nyz_*i;
    for(int j=0; j<nq_; ++j)
      jacmap[nyz_+j+nx_*i] = nyz_*(1+ns) + j + nq_*i;
  }
  integrator.setOption("jacmap",jacmap);

  // Initial value for the integrators
  vector<double> jacinit(nx_*ns,0.0);
  if(iind==INTEGRATOR_X0){
    for(int i=0; i<1+ns; ++i)
      jacinit[i+nx_*i] = 1;
    integrator.setOption("jacinit",jacinit);
  }

  // Pass linear solver
  if(!linsol_.isNull()){
    LinearSolver linsol_aug = shared_cast<LinearSolver>(linsol_.clone());
    linsol_aug->nrhs_ = 1+ns;
    integrator.setLinearSolver(linsol_aug,jac_);
  }
  
  return integrator;
}

FX IdasInternal::getJacobian(){
  // Quick return if already created
  if(!jac_.isNull())
    return jac_;

  SXFunction f = shared_cast<SXFunction>(f_);
  if(f.isNull())
    throw CasadiException("IdasInternal::getJacobian(): Not an SXFunction");
  
  // Get the Jacobian in the Newton iteration
  SX cj("cj");
  SXMatrix jac = horzcat(f.jac(DAE_Y,DAE_RES) + cj*f.jac(DAE_YDOT,DAE_RES),f.jac(DAE_Z,DAE_RES));

  // Jacobian function
  vector<vector<SX> > jac_in(JAC_NUM_IN);
  jac_in[JAC_T] = f->inputv.at(DAE_T);
  jac_in[JAC_Y] = f->inputv.at(DAE_Y);
  jac_in[JAC_YDOT] = f->inputv.at(DAE_YDOT);
  jac_in[JAC_Z] = f->inputv.at(DAE_Z);
  jac_in[JAC_P] = f->inputv.at(DAE_P);
  jac_in[JAC_CJ] = vector<SX>(1,cj);
  SXFunction J(jac_in,jac);
      
  // Save function
  jac_ = J;
  
  return J;
}

  
LinearSolver IdasInternal::getLinearSolver(){
  return linsol_;
}
  
void IdasInternal::copyNV(const Matrix<double>& x, const Matrix<double>& xp, const Matrix<double>& z, N_Vector& yz, N_Vector& yP, N_Vector& yQ){
  double *yzd = NV_DATA_S(yz);
  double *yPd = NV_DATA_S(yP);
  
  copy(x.begin(), x.begin()+ny_, yzd);
  copy(xp.begin(), xp.begin()+ny_, yPd);
  copy(z.begin(), z.end(), yzd+ny_);

  if(nq_>0){
    double *yQd = NV_DATA_S(yQ);
    copy(x.begin()+ny_, x.end(), yQd);
  }
}
  
void IdasInternal::copyNV(const N_Vector& yz, const N_Vector& yP, const N_Vector& yQ, Matrix<double>& x, Matrix<double>& xp, Matrix<double>& z){
  const double *yzd = NV_DATA_S(yz);
  const double *ypd = NV_DATA_S(yP);
  
  copy(yzd,yzd+ny_,x.begin());
  copy(yzd+ny_,yzd+ny_+nz_,z.begin());
  copy(ypd,ypd+ny_,xp.begin());
  
  if(nq_>0){
    const double *yQd = NV_DATA_S(yQ);
    copy(yQd,yQd+nq_, x.begin()+ny_);
  }
}

void IdasInternal::getAdjointSeeds(){
  for(int i=0; i<nadir_; ++i){
    const double *x0 = &adjSeed(INTEGRATOR_XF,i)[0];
    const double *xp0 = &adjSeed(INTEGRATOR_XPF,i)[0];
    const double *z0 = &adjSeed(INTEGRATOR_ZF,i)[0];

    double *yz = NV_DATA_S(yzB_[i]);
    double *yp = NV_DATA_S(yPB_[i]);

    copy(x0,x0+ny_,yz);
    copy(xp0,xp0+ny_,yp);
    copy(z0,z0+nz_,yz+ny_);
  }
}

void IdasInternal::setAdjointSensitivities(){
  for(int i=0; i<nadir_; ++i){
    double *xf = &adjSens(INTEGRATOR_X0,i)[0];
    double *xpf = &adjSens(INTEGRATOR_XP0,i)[0];
    double *zf = &adjSens(INTEGRATOR_Z0,i)[0];
    
    const double *yz = NV_DATA_S(yzB_[i]);
    const double *yp = NV_DATA_S(yPB_[i]);
  
    copy(yz,yz+ny_,xf);
    copy(yp,yp+ny_,xpf);
    copy(yz+ny_,yz+ny_+nz_,zf);
  }
}

} // namespace Sundials
} // namespace CasADi

