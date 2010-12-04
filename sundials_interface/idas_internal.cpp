/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A minimalistic computer algebra system with automatic differentiation 
 *    and framework for dynamic optimization.
 *    Copyright (C) 2010 by Joel Andersson, Moritz Diehl et al., K.U.Leuven. All rights reserved.
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
#include <cassert>

#ifdef PRECOND_TEST
  // DENSE PRECONDITIONER: REMOVE

// LU-Factorize dense matrix (lapack)
extern "C" void dgetrf_(int *m, int *n, double *a, int *lda, int *ipiv, int *info);

// Solve a system of equation using an LU-factorized matrix (lapack)
extern "C" void dgetrs_(char* trans, int *n, int *nrhs, double *a, int *lda, int *ipiv, double *b, int *ldb, int *info);

// QR-factorize dense matrix (lapack)
extern "C" void dgeqrf_(int *m, int *n, double *a, int *lda, double *tau, double *work, int *lwork, int *info);

// Multiply right hand side with Q-transpose (lapack)
extern "C" void dormqr_(char *side, char *trans, int *n, int *m, int *k, double *a, int *lda, double *tau, double *c, int *ldc, double *work, int *lwork, int *info);

// Solve upper triangular system (lapack)
extern "C" void dtrsm_(char *side, char *uplo, char *transa, char *diag, int *m, int *n, double *alpha, double *a, int *lda, double *b, int *ldb);

#endif // PRECOND_TEST

// use QR factorization instead of LU
bool use_qr = false;

using namespace std;
namespace CasADi{
namespace Sundials{

int IdasInternal::getNX(const FX& f, const FX& q){
  // Number of states
  int nx = f.output().numel();
  
  // Add quadratures, if any_
  if(!q.isNull()) nx += q.output().numel();
  
  return nx;
}

int IdasInternal::getNP(const FX& f){
  return f.input(DAE_P).numel();
}

IdasInternal::IdasInternal(const FX& f, const FX& q) : IntegratorInternal(getNX(f,q), getNP(f)), f_(f), q_(q){
  addOption("suppress_algebraic",          OT_BOOLEAN, false); // supress algebraic variables in the error testing
  addOption("calc_ic",                     OT_BOOLEAN, true);  // use IDACalcIC to get consistent initial conditions
  addOption("abstolv",                     OT_REALVECTOR, Option());
  addOption("fsens_abstolv",               OT_REALVECTOR, Option()); 
  addOption("dense_preconditioner",        OT_BOOLEAN, false); // precondition with a dense preconditioner
  addOption("max_step_size",               OT_REAL, 0); // maximim step size
  addOption("first_time",                  OT_REAL, 1.0); // first requested time (to pass to CalcIC)

  mem_ = 0;

  y0_  = y_  = 0; 
  yp0_ = yp_ = 0, 
  yQ0_ = yQ_ = 0;
  
  is_init = false;

  // Get dimensions
  ny_ = f.output().numel();
  nq_ = q.isNull() ? 0 : q.output().numel();

}

IdasInternal::~IdasInternal(){ 
  if(mem_) IDAFree(&mem_);

  // N-vectors for the DAE integration
  if(y0_) N_VDestroy_Serial(y0_);
  if(y_) N_VDestroy_Serial(y_);
  if(yp0_) N_VDestroy_Serial(yp0_);
  if(yp_) N_VDestroy_Serial(yp_);
  if(yQ0_) N_VDestroy_Serial(yQ0_);
  if(yQ_) N_VDestroy_Serial(yQ_);
  
    // Forward problem
  for(vector<N_Vector>::iterator it=yS0_.begin(); it != yS0_.end(); ++it)   N_VDestroy_Serial(*it);
  for(vector<N_Vector>::iterator it=yS_.begin(); it != yS_.end(); ++it)     N_VDestroy_Serial(*it);
  for(vector<N_Vector>::iterator it=ypS0_.begin(); it != ypS0_.end(); ++it) N_VDestroy_Serial(*it);
  for(vector<N_Vector>::iterator it=ypS_.begin(); it != ypS_.end(); ++it)   N_VDestroy_Serial(*it);
  for(vector<N_Vector>::iterator it=yQS0_.begin(); it != yQS0_.end(); ++it) N_VDestroy_Serial(*it);
  for(vector<N_Vector>::iterator it=yQS_.begin(); it != yQS_.end(); ++it)   N_VDestroy_Serial(*it);
  
  // Adjoint problem
  for(vector<N_Vector>::iterator it=yB0_.begin(); it != yB0_.end(); ++it)   N_VDestroy_Serial(*it);
  for(vector<N_Vector>::iterator it=yB_.begin(); it != yB_.end(); ++it)     N_VDestroy_Serial(*it);
  for(vector<N_Vector>::iterator it=ypB0_.begin(); it != ypB0_.end(); ++it) N_VDestroy_Serial(*it);
  for(vector<N_Vector>::iterator it=ypB_.begin(); it != ypB_.end(); ++it)   N_VDestroy_Serial(*it);
//  for(vector<N_Vector>::iterator it=yQB0_.begin(); it != yQB0_.end(); ++it) N_VDestroy_Serial(*it);
  for(vector<N_Vector>::iterator it=yQB_.begin(); it != yQB_.end(); ++it)   N_VDestroy_Serial(*it);

}

void IdasInternal::init(){
  // Call the base class init
  IntegratorInternal::init();

  // Init ODE rhs function and quadrature functions, jacobian function
  f_.init();
  if(!q_.isNull()) q_.init();
  if(!jacx_.isNull()) jacx_.init();
  
  // Get the number of forward and adjoint directions
  nfdir_f_ = f_.getOption("number_of_fwd_dir").toInt();
  nadir_f_ = f_.getOption("number_of_adj_dir").toInt();
  nfdir_q_ = q_.isNull() ? 0 : q_.getOption("number_of_fwd_dir").toInt();
  nadir_q_ = q_.isNull() ? 0 : q_.getOption("number_of_adj_dir").toInt();

  // Quick return if already initialized
  if(is_init){
    reset(ad_order_,ad_order_);
    return;
  }

  calc_ic_ = getOption("calc_ic").toInt();
  
  // Sundials return flag
  int flag;

  // Create IDAS memory block
  mem_ = IDACreate();
  if(mem_==0) throw CasadiException("IDACreate(): Creation failed");

  // Allocate n-vectors for ivp
  y0_ = N_VMake_Serial(ny_,&input(INTEGRATOR_X0).data()[0]);
  yp0_ = N_VMake_Serial(ny_,&input(INTEGRATOR_XP0).data()[0]);
  y_ = N_VMake_Serial(ny_,&output(INTEGRATOR_XF).data()[0]);
  yp_ = N_VMake_Serial(ny_,&output(INTEGRATOR_XPF).data()[0]);

  // Initialize Idas
  double t0 = 0;
  N_VConst(0.0, y0_);
  N_VConst(0.0, yp0_);
  IDAInit(mem_, res_wrapper, t0, y0_, yp0_);

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
    // Find out which components are differential
    const vector<int>& is_differential = getOption("is_differential").toIntVector();
    copy(is_differential.begin(),is_differential.end(),NV_DATA_S(yp_));
    
    // Pass this information to IDAS
    flag = IDASetId(mem_, yp_);
    if(flag != IDA_SUCCESS) idas_error("IDASetId",flag);
    
    // Reset the vector
    N_VConst(0.0, yp_);
  }

  // attach a linear solver
  if(getOption("linear_solver")=="dense"){
    // Dense jacobian
    flag = IDADense(mem_, ny_);
    if(flag != IDA_SUCCESS) idas_error("IDADense",flag);
    if(exact_jacobian_){
      // Pass to IDA
      flag = IDADlsSetDenseJacFn(mem_, djac_wrapper);
      if(flag!=IDA_SUCCESS) idas_error("IDADlsSetDenseJacFn",flag);
    }
  } else if(getOption("linear_solver")=="banded") {
    // Banded jacobian
    flag = IDABand(mem_, ny_, getOption("upper_bandwidth").toInt(), getOption("lower_bandwidth").toInt());
    if(flag != IDA_SUCCESS) idas_error("IDABand",flag);
    
    // Banded Jacobian information
    if(exact_jacobian_) throw CasadiException("IDAS: Banded Jacobian information not implemented");

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
    if(getOption("exact_jacobian")==true){
      flag = IDASpilsSetJacTimesVecFn(mem_, jtimes_wrapper);
      if(flag != IDA_SUCCESS) idas_error("IDASpilsSetJacTimesVecFn",flag);
    }

    // Attach functions for jacobian information
    if(exact_jacobian_) throw CasadiException("IDAS: Iterataive Jacobian information not implemented");
    
    // Add a preconditioner
    if(getOption("dense_preconditioner")==true){
      // Form and jacobian function df/dy + cj* df/dydot, if it has not already been done
      if(jacx_.isNull()){
        throw CasadiException("IDASInternal: no jacobian function supplied");
      }

#if PRECOND_TEST
#error
      // REMOVE THIS!!
      
      // Allocate jacobian
      pc_.resize(ny_*ny_);
      
      // Allocate info needed for the linear solver
      if(use_qr){
        tau_.resize(ny_);
        work_.resize(10*ny_);
      } else {
        ipiv_.resize(ny_);
      }
#else
      // Make sure that a Jacobian has been provided
      if(jacx_.isNull()) throw CasadiException("IdasInternal::init(): No Jacobian has been provided.");

      // Make sure that a linear solver has been providided
      if(linsol_.isNull()) throw CasadiException("IdasInternal::init(): No user defined linear solver has been provided.");
#endif

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
    yQ0_ = N_VMake_Serial(nq_,&input(INTEGRATOR_X0).data()[ny_]);
    yQ_ = N_VMake_Serial(nq_,&output(INTEGRATOR_XF).data()[ny_]);

    // Initialize quadratures in IDAS
    flag = IDAQuadInit(mem_, rhsQ_wrapper, yQ0_);
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
    
 // Sensitivities
 if(ad_order_>0){
   
   // Forward sensitivity problem
   if(nfdir_>0){
     // Allocate n-vectors
     yS0_.resize(nfdir_);
     yS_.resize(nfdir_);
     ypS0_.resize(nfdir_);
     ypS_.resize(nfdir_);
     for(int i=0; i<nfdir_; ++i){
        yS0_[i] = N_VMake_Serial(ny_,&input(INTEGRATOR_X0).dataF(i)[0]);
        yS_[i] = N_VMake_Serial(ny_,&output(INTEGRATOR_XF).dataF(i)[0]);
        ypS0_[i] = N_VMake_Serial(ny_,&input(INTEGRATOR_XP0).dataF(i)[0]);
        ypS_[i] = N_VMake_Serial(ny_,&output(INTEGRATOR_XPF).dataF(i)[0]);
      }

      // Allocate n-vectors for quadratures
      if(nq_>0){
        yQS0_.resize(nfdir_);
        yQS_.resize(nfdir_);
        for(int i=0; i<nfdir_; ++i){
          yQS0_[i] = N_VMake_Serial(nq_,&input(INTEGRATOR_X0).dataF(i)[ny_]);
          yQS_[i] = N_VMake_Serial(nq_,&output(INTEGRATOR_XF).dataF(i)[ny_]);
        }
      }
     
    // Get the sensitivity method
    if(getOption("sensitivity_method")== "simultaneous") ism_ = IDA_SIMULTANEOUS;
    else if(getOption("sensitivity_method")=="staggered") ism_ = IDA_STAGGERED;
    else throw CasadiException("IDAS: Unknown sensitivity method");

    // Initialize forward sensitivities
    if(finite_difference_fsens_){
      // Use finite differences to calculate the residual in the forward sensitivity equations
      flag = IDASensInit(mem_,nfdir_,ism_,0,&yS0_[0],&ypS0_[0]);
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
      flag = IDASetSensParams(mem_,&input(INTEGRATOR_P).data()[0],pbar,plist);
      if(flag != IDA_SUCCESS) idas_error("IDASetSensParams",flag);


      
      //  IDASetSensDQMethod

    } else {
      // Use AD to calculate the residual in the forward sensitivity equations
      flag = IDASensInit(mem_,nfdir_,ism_,resS_wrapper,&yS0_[0],&ypS0_[0]);
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
      flag = IDAQuadSensInit(mem_, rhsQS_wrapper, &yQS0_[0]);
      if(flag != IDA_SUCCESS) idas_error("IDAQuadSensInit",flag);
      
      // Set tolerances
      flag = IDAQuadSensSStolerances(mem_,fsens_reltol_,&fsens_abstol[0]);
      if(flag != IDA_SUCCESS) idas_error("IDAQuadSensSStolerances",flag);
    }
  } // enable fsens

  // Adjoint sensitivity problem
  whichB_.resize(nadir_);

  // Allocate n-vectors
  yB0_.resize(nadir_);
  yB_.resize(nadir_);
  ypB0_.resize(nadir_);
  ypB_.resize(nadir_);
  for(int i=0; i<nadir_; ++i){
    yB0_[i] = N_VMake_Serial(ny_,&output(INTEGRATOR_XF).dataA(i)[0]);
    yB_[i] = N_VMake_Serial(ny_,&input(INTEGRATOR_X0).dataA(i)[0]);
    ypB0_[i] = N_VMake_Serial(ny_,&output(INTEGRATOR_XPF).dataA(i)[0]);
    ypB_[i] = N_VMake_Serial(ny_,&input(INTEGRATOR_XP0).dataA(i)[0]);
  }

  // Allocate n-vectors for quadratures
  yQB_.resize(nadir_);
  for(int i=0; i<nadir_; ++i){
    //yQB0_[i] = N_VNew_Serial(np_);
    yQB_[i] = N_VMake_Serial(np_,&input(INTEGRATOR_P).dataA(i)[0]);
  }

  if(nadir_>0){
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
 } // ad_order>0
 
 is_init = true;
 isInitAdj_ = false;
}

void IdasInternal::initAdj(){
  int flag;
  
  for(int dir=0; dir<nadir_; ++dir){

    // Create backward problem
    flag = IDACreateB(mem_, &whichB_[dir]);
    if(flag != IDA_SUCCESS) idas_error("IDACreateB",flag);
    
    // Initialize the backward problem
    double tB0 = input(INTEGRATOR_TF).data()[0];
    flag = IDAInitB(mem_, whichB_[dir], resB_wrapper, tB0, yB0_[dir], ypB0_[dir]);
    if(flag != IDA_SUCCESS) idas_error("IDAInitB",flag);

    // Set tolerances
    flag = IDASStolerancesB(mem_, whichB_[dir], asens_reltol_, asens_abstol_);
    if(flag!=IDA_SUCCESS) idas_error("IDASStolerancesB",flag);

    // User data
    flag = IDASetUserDataB(mem_, whichB_[dir], this);
    if(flag != IDA_SUCCESS) idas_error("IDASetUserDataB",flag);
      
    // attach linear solver
    if(getOption("asens_linear_solver")=="dense"){
      // Dense jacobian
      flag = IDADenseB(mem_, whichB_[dir], ny_);
      if(flag != IDA_SUCCESS) idas_error("IDADenseB",flag);
    } else if(getOption("asens_linear_solver")=="banded") {
      // Banded jacobian
      flag = IDABandB(mem_, whichB_[dir], ny_, getOption("asens_upper_bandwidth").toInt(), getOption("asens_lower_bandwidth").toInt());
      if(flag != IDA_SUCCESS) idas_error("IDABand",flag);
      //    if(exact_jacobian) sundialsAssert(CVDlsSetBandJacFn(cvode_mem_[k], bjac_static));  
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
    N_VConst(0.0, yQB_[dir]);
    flag = IDAQuadInitB(mem_,whichB_[dir],rhsQB_wrapper,yQB_[dir]);
    if(flag!=IDA_SUCCESS) idas_error("CVodeQuadInitB",flag);
  
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



void IdasInternal::res(double t, const double* y, const double* yp, double* r){
  // Get time
  time1 = clock();
  
   // Pass input
   f_.input(DAE_T).set(t);
   f_.input(DAE_Y).set(y);
   f_.input(DAE_YDOT).set(yp);
   f_.input(DAE_P).set(input(INTEGRATOR_P).data());

    // Evaluate
   f_.evaluate();
    
    // Get results
   f_.output().get(r);

   time2 = clock();
   t_res += double(time2-time1)/CLOCKS_PER_SEC;
   
   // Check the result for consistency
   for(int i=0; i<ny_; ++i)
     if(isnan(r[i]) || isinf(r[i])){
       throw 1;
     }
}

int IdasInternal::res_wrapper(double t, N_Vector yy, N_Vector yp, N_Vector rr, void *user_data){
 try{
    IdasInternal *this_ = (IdasInternal*)user_data;
    this_->res(t,NV_DATA_S(yy),NV_DATA_S(yp),NV_DATA_S(rr));
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

void IdasInternal::jtimes(double t, const double *yy, const double *yp, const double *rr, const double *v, double *Jv, double cj, double *tmp1, double *tmp2){
  // Get time
  time1 = clock();
  
   // Pass input
   f_.setInput(t,DAE_T);
   f_.setInput(yy,DAE_Y);
   f_.setInput(yp,DAE_YDOT);
     
   // Pass seeds of the state vectors
   f_.setFwdSeed(v,DAE_Y);
   
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

int IdasInternal::jtimes_wrapper(double t, N_Vector yy, N_Vector yp, N_Vector rr, N_Vector v, N_Vector Jv, double cj, void *user_data, N_Vector tmp1, N_Vector tmp2){
 try{
    IdasInternal *this_ = (IdasInternal*)user_data;
    this_->jtimes(t,NV_DATA_S(yy),NV_DATA_S(yp),NV_DATA_S(rr),NV_DATA_S(v),NV_DATA_S(Jv),cj,NV_DATA_S(tmp1),NV_DATA_S(tmp2));
    return 0;
  } catch(exception& e){
    cerr << "jtimes failed: " << e.what() << endl;;
    return 1;
  }  
}

void IdasInternal::resS(int Ns, double t, const double* yy, const double* yp, const double *resval, N_Vector *yS, N_Vector* ypS, N_Vector *resvalS, double *tmp1, double *tmp2, double *tmp3){
  assert(Ns==nfdir_);

  // Record the current cpu time
  time1 = clock();
  
   // Pass input
   f_.input(DAE_T).set(t);
   f_.input(DAE_Y).set(yy);
   f_.input(DAE_YDOT).set(yp);
   f_.input(DAE_P).set(input(INTEGRATOR_P).data());

   // Calculate the forward sensitivities, nfdir_f_ directions at a time
   for(int j=0; j<nfdir_; j += nfdir_f_){
     for(int dir=0; dir<nfdir_f_ && j+dir<nfdir_; ++dir){
       // Pass forward seeds 
       f_.setFwdSeed(0.0,DAE_T,dir);
       f_.setFwdSeed(NV_DATA_S(yS[j+dir]),DAE_Y,dir);
       f_.setFwdSeed(NV_DATA_S(ypS[j+dir]),DAE_YDOT,dir);
       f_.setFwdSeed(input(INTEGRATOR_P).dataF(j+dir),DAE_P,dir);
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

int IdasInternal::resS_wrapper(int Ns, double t, N_Vector yy, N_Vector yp, N_Vector resval, N_Vector *yS, N_Vector *ypS, N_Vector *resvalS, void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3){
 try{
    IdasInternal *this_ = (IdasInternal*)user_data;
    this_->resS(Ns,t,NV_DATA_S(yy),NV_DATA_S(yp),NV_DATA_S(resval),yS,ypS,resvalS,NV_DATA_S(tmp1),NV_DATA_S(tmp2),NV_DATA_S(tmp3));
    return 0;
  } catch(exception& e){
    cerr << "resS failed: " << e.what() << endl;;
    return 1;
  }
}

void IdasInternal::reset(int fsens_order, int asens_order){
  // Reset timers
  t_res = t_fres = t_jac = t_psolve = t_psetup_jac = t_psetup_fac = 0;
    
  fsens_order_ = fsens_order;
  asens_order_ = asens_order;
  
  // Get the time horizon
  double t0 = input(INTEGRATOR_T0).data()[0];
  double tf = input(INTEGRATOR_TF).data()[0];
  t_ = t0;
  
  // Return flag
  int flag;
  
  // Re-initialize
  flag = IDAReInit(mem_, t0, y0_, yp0_);
  if(flag != IDA_SUCCESS) idas_error("IDAReInit",flag);

  // Re-initialize quadratures
  if(nq_>0){
    flag = IDAQuadReInit(mem_, yQ0_);
    if(flag != IDA_SUCCESS) idas_error("IDAQuadReInit",flag);
  }
  
  // Re-initialize sensitivities
  if(fsens_order_>0){
    flag = IDASensReInit(mem_,ism_,&yS0_[0],&ypS0_[0]);
    if(flag != IDA_SUCCESS) idas_error("IDASensReInit",flag);

    if(nq_>0){
      flag = IDAQuadSensReInit(mem_, &yQS0_[0]);
      if(flag != IDA_SUCCESS) idas_error("IDAQuadSensReInit",flag);
    }
  } else {
    // Turn of sensitivities
    flag = IDASensToggleOff(mem_);
    if(flag != IDA_SUCCESS) idas_error("IDASensToggleOff",flag);
  }

  calc_ic_ = getOption("calc_ic").toInt();
  if(calc_ic_){
    // Try to calculate consistent initial values
    cout << "Starting IDACalcIC, t = " << t0 << endl;

    double t_scale = getOption("first_time").toDouble();
    flag = IDACalcIC(mem_, IDA_YA_YDP_INIT, t_scale);
//    flag = IDACalcIC(mem_, IDA_Y_INIT, tf);
    if(flag != IDA_SUCCESS) idas_error("IDACalcIC",flag);

    // Retrieve the initial values
    flag = IDAGetConsistentIC(mem_, y0_, yp0_);
    if(flag != IDA_SUCCESS) idas_error("IDAGetConsistentIC",flag);
  }
  
}
  
void IdasInternal::integrate(double t_out){
  int flag;
  
  // tolerance
  double ttol = 1e-9;
  if(fabs(t_-t_out)<ttol){
    copy(input(INTEGRATOR_X0).data().begin(),input(INTEGRATOR_X0).data().end(),output(INTEGRATOR_XF).data().begin());
    copy(input(INTEGRATOR_XP0).data().begin(),input(INTEGRATOR_XP0).data().end(),output(INTEGRATOR_XPF).data().begin());
    if(fsens_order_>0){
      for(int i=0; i<nfdir_; ++i){
        copy(input(INTEGRATOR_X0).dataF(i).begin(),input(INTEGRATOR_X0).dataF(i).end(),output(INTEGRATOR_XF).dataF(i).begin());
        copy(input(INTEGRATOR_XP0).dataF(i).begin(),input(INTEGRATOR_XP0).dataF(i).end(),output(INTEGRATOR_XPF).dataF(i).begin());
      }
    }
    return;
  }

  if(asens_order_>0){
    int ncheck; // number of checkpoints stored so far
    flag = IDASolveF(mem_, t_out, &t_, y_, yp_, IDA_NORMAL, &ncheck);
    if(flag != IDA_SUCCESS && flag != IDA_TSTOP_RETURN) idas_error("IDASolveF",flag);
    
  } else {
    flag = IDASolve(mem_, t_out, &t_, y_, yp_, IDA_NORMAL);
    if(flag != IDA_SUCCESS && flag != IDA_TSTOP_RETURN) idas_error("IDASolve",flag);
  }

  if(nq_>0){
    double tret;
    flag = IDAGetQuad(mem_, &tret, yQ_);
    if(flag != IDA_SUCCESS) idas_error("IDAGetQuad",flag);
  }
  
  if(fsens_order_>0){
    // Get the sensitivities
    flag = IDAGetSens(mem_,&t_, &yS_[0]);
    if(flag != IDA_SUCCESS) idas_error("IDAGetSens",flag);
    
    if(nq_>0){
      double tret;
      flag = IDAGetQuadSens(mem_, &tret, &yQS_[0]);
      if(flag != IDA_SUCCESS) idas_error("IDAGetQuadSens",flag);
    }
  }
}

void IdasInternal::resetAdj(){
  double tf = input(INTEGRATOR_TF).data()[0];

  int flag;

  for(int dir=0; dir<nadir_; ++dir){
    if(isInitAdj_){
      flag = IDAReInitB(mem_, whichB_[dir], tf, yB0_[dir], ypB0_[dir]);
      if(flag != IDA_SUCCESS) idas_error("IDAReInitB",flag);
      
      N_VConst(0.0,yQB_[dir]);
      flag = IDAQuadReInitB(mem_,whichB_[dir],yQB_[dir]);
      if(flag!=IDA_SUCCESS) idas_error("IDAQuadReInitB",flag);

    } else {
      // Initialize the adjoint integration
      initAdj();
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
    flag = IDAGetB(mem_, whichB_[dir], &tret, yB_[dir], ypB_[dir]);
    if(flag!=IDA_SUCCESS) idas_error("IDAGetB",flag);

    flag = IDAGetQuadB(mem_, whichB_[dir], &tret, yQB_[dir]);
    if(flag!=IDA_SUCCESS) idas_error("IDAGetQuadB",flag);
  }
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

    stream << "Time spent in the DAE residual: " << t_res << " s." << endl;
    stream << "Time spent in the forward sensitivity residual: " << t_fres << " s." << endl;
    stream << "Time spent in the jacobian function or jacobian times vector function: " << t_jac << " s." << endl;
    stream << "Time spent in the preconditioner solve function: " << t_psolve << " s." << endl;
    stream << "Time spent to generate the jacobian in the preconditioner setup function: " << t_psetup_jac << " s." << endl;
    stream << "Time spent to factorize the jacobian in the preconditioner setup function: " << t_psetup_fac << " s." << endl;
    
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

int IdasInternal::rhsQ_wrapper(double t, N_Vector yy, N_Vector yp, N_Vector rhsQ, void *user_data){
 try{
    IdasInternal *this_ = (IdasInternal*)user_data;
    this_->rhsQ(t,NV_DATA_S(yy),NV_DATA_S(yp),NV_DATA_S(rhsQ));
    return 0;
  } catch(exception& e){
    cerr << "rhsQ failed: " << e.what() << endl;;
    return 1;
  }
}

void IdasInternal::rhsQ(double t, const double* yy, const double* yp, double* rhsQ){
   // Pass input
   q_.input(DAE_T).set(t);
   q_.input(DAE_Y).set(yy);
   q_.input(DAE_YDOT).set(yp);
   q_.input(DAE_P).set(input(INTEGRATOR_P).data());

    // Evaluate
   q_.evaluate();
    
    // Get results
   q_.output().get(rhsQ);
}
  
void IdasInternal::rhsQS(int Ns, double t, N_Vector yy, N_Vector yp, N_Vector *yyS, N_Vector *ypS, N_Vector rrQ, N_Vector *rhsvalQS, 
                        N_Vector tmp1, N_Vector tmp2, N_Vector tmp3){
  assert(Ns==nfdir_);

  // Pass input
   q_.input(DAE_T).set(t);
   q_.input(DAE_Y).set(NV_DATA_S(yy));
   q_.input(DAE_YDOT).set(NV_DATA_S(yp));
   q_.input(DAE_P).set(input(INTEGRATOR_P).data());
     
   // Pass forward seeds
  for(int i=0; i<nfdir_; ++i){
    q_.input(DAE_T).setF(0.0);
    q_.input(DAE_Y).setF(NV_DATA_S(yyS[i]));
    q_.input(DAE_YDOT).setF(NV_DATA_S(ypS[i]));
    q_.input(DAE_P).setF(input(INTEGRATOR_P).dataF(i));
   
    // Evaluate the AD forward algorithm
    q_.evaluate(1,0);
      
    // Get the output seeds
    q_.output().getF(NV_DATA_S(rhsvalQS[i]));
  }
}

int IdasInternal::rhsQS_wrapper(int Ns, double t, N_Vector yy, N_Vector yp, N_Vector *yyS, N_Vector *ypS, N_Vector rrQ, N_Vector *rhsvalQS, void *user_data, 
                                N_Vector tmp1, N_Vector tmp2, N_Vector tmp3){
  
 try{
    IdasInternal *this_ = (IdasInternal*)user_data;
    this_->rhsQS(Ns,t,yy,yp,yyS,ypS,rrQ,rhsvalQS,tmp1,tmp2,tmp3);
    return 0;
  } catch(exception& e){
    cerr << "rhsQS failed: " << e.what() << endl;;
    return 1;
  }
}

void IdasInternal::resB(double t, const double* y, const double* yp, const double* yB, const double* ypB, double* resvalB){
  // Pass input
  f_.input(DAE_T).set(t);
  f_.input(DAE_Y).set(y);
  f_.input(DAE_YDOT).set(yp);
  f_.input(DAE_P).set(input(INTEGRATOR_P).data());
  
  // Pass adjoint seeds
  f_.output(DAE_RES).setA(yB);

  // Evaluate
  f_.evaluate(0,1);

  // Save to output
  f_.input(DAE_Y).getA(resvalB);

  // Pass adjoint seeds
  f_.output(DAE_RES).setA(ypB);
  
  // Evaluate AD adjoint
  f_.evaluate(0,1);

  // Save to output
  const vector<double>& asens_ydot = f_.input(DAE_YDOT).dataA();
  for(int i=0; i<ny_; ++i)
    resvalB[i] -= asens_ydot[i];
  
  // If quadratures are included
  if(nq_>0){
    // Pass input to quadratures
    q_.input(DAE_T).set(t);
    q_.input(DAE_Y).set(y);
    q_.input(DAE_YDOT).set(yp);
    q_.input(DAE_P).set(input(INTEGRATOR_P).data());

    // Pass adjoint seeds
    q_.output(DAE_RES).setA(&output(INTEGRATOR_XF).dataA()[ny_]);

    // Evaluate
    q_.evaluate(0,1);
    
    // Get the input seeds
    const vector<double>& asens_q = q_.input(DAE_Y).dataA();
//    cout << "asens_q = " << asens_q  << endl;
    
    // Copy to result
    for(int i=0; i<ny_; ++i)
      resvalB[i] += asens_q[i];
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
  
void IdasInternal::rhsQB(double t, const double* y, const double* yp, const double* yB, const double* ypB, double *rhsvalBQ){
  // Pass input
  f_.input(DAE_T).set(t);
  f_.input(DAE_Y).set(y);
  f_.input(DAE_YDOT).set(yp);
  f_.input(DAE_P).set(input(INTEGRATOR_P).data());

  // Pass adjoint seeds
  f_.output(DAE_RES).setA(yB);

  // Evaluate
  f_.evaluate(0,1);

  // Save to output
  f_.input(DAE_P).getA(rhsvalBQ);
  
  // If quadratures are included
  if(nq_>0){
    // Pass input to quadratures
    q_.input(DAE_T).set(t);
    q_.input(DAE_Y).set(y);
    q_.input(DAE_YDOT).set(yp);
    q_.input(DAE_P).set(input(INTEGRATOR_P).data());

    // Pass adjoint seeds
    q_.output(DAE_RES).setA(&output(INTEGRATOR_XF).dataA()[ny_]);

    // Evaluate
    q_.evaluate(0,1);
    
    // Get the input seeds
    const vector<double>& qres = q_.input(DAE_P).dataA();
    
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

void IdasInternal::djac(int Neq, double t, double cj, N_Vector yy, N_Vector yp, N_Vector rr, DlsMat Jac, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3){
cout << "generating jac at t = " << t << endl;
  
  // Get time
  time1 = clock();

  // Pass input to the jacobian function
  jacx_.setInput(t,JAC_T);
  jacx_.setInput(NV_DATA_S(yy),JAC_Y);
  jacx_.setInput(NV_DATA_S(yp),JAC_YDOT);
  jacx_.setInput(input(INTEGRATOR_P).data(),JAC_P);
  jacx_.setInput(cj,JAC_CJ);

  // Evaluate jacobian
  jacx_.evaluate();
  
  // Get the result
  const vector<double> &res = jacx_.getOutputData(JAC_J);
  for(int i=0; i<ny_; ++i){
    for(int j=0; j<ny_; ++j){
      DENSE_ELEM(Jac, i, j) = res[j + i*ny_];
    }
  }
  
  // Log time duration
  time2 = clock();
  t_jac += double(time2-time1)/CLOCKS_PER_SEC;
}

int IdasInternal::djac_wrapper(int Neq, double t, double cj, N_Vector yy, N_Vector yp, N_Vector rr, DlsMat Jac, void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3){
 try{
    IdasInternal *this_ = (IdasInternal*)user_data;
    this_->djac(Neq, t, cj, yy, yp, rr, Jac, tmp1, tmp2, tmp3);
    return 0;
  } catch(exception& e){
    cerr << "djac failed: " << e.what() << endl;;
    return 1;
  }
}



void IdasInternal::setStopTime(double tf){
  // Set the stop time of the integration -- don't integrate past this point
  int flag = IDASetStopTime(mem_, tf);
  if(flag != IDA_SUCCESS) idas_error("IDASetStopTime",flag);
}

int IdasInternal::psolve_wrapper(double t, N_Vector yy, N_Vector yp, N_Vector rr, N_Vector rvec, N_Vector zvec, double cj, double delta, void *user_data, N_Vector tmp){
 try{
    IdasInternal *this_ = (IdasInternal*)user_data;
    assert(this_);
    this_->psolve(t, yy, yp, rr, rvec, zvec, cj, delta, tmp);
    return 0;
  } catch(exception& e){
    cerr << "psolve failed: " << e.what() << endl;;
    return 1;
  }
}

int IdasInternal::psetup_wrapper(double t, N_Vector yy, N_Vector yp, N_Vector rr, double cj, void* user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3){
 try{
    IdasInternal *this_ = (IdasInternal*)user_data;
    assert(this_);
    this_->psetup(t, yy, yp, rr, cj, tmp1, tmp2, tmp3);
    return 0;
  } catch(exception& e){
    cerr << "psetup failed: " << e.what() << endl;;
    return 1;
  }
}

void IdasInternal::psolve(double t, N_Vector yy, N_Vector yp, N_Vector rr, N_Vector rvec, N_Vector zvec, double cj, double delta, N_Vector tmp){
  // Get time
  time1 = clock();

  #ifdef PRECOND_TEST
  // DENSE PRECONDITIONER: REMOVE

  // Copy input to output vector
  N_VScale(1.0,rvec,zvec);

  int nrhs = 1; // number of right hand sides
  int info = 100;
  
  if(use_qr){
    char trans = 'T';
    char side = 'L';
    int k = 10; // ????
    int lwork = work_.size();
    dormqr_(&side, &trans, &ny_, &nrhs, &k, &pc_[0], &ny_, &tau_[0], NV_DATA_S(zvec), &ny_, &work_[0], &lwork, &info);
    
    char uplo = 'U';
    char transa = 'N';
    char diag = 'N';
    double alpha = 1.;
    dtrsm_(&side, &uplo, &transa, &diag, &ny_, &nrhs, &alpha, &pc_[0], &ny_, NV_DATA_S(zvec), &ny_);

    
  } else {
    char trans = 'N';
    dgetrs_(&trans, &ny_, &nrhs, &pc_[0], &ny_, &ipiv_[0], NV_DATA_S(zvec), &ny_, &info);
  }
  
  if(info != 0)
    throw CasadiException("failed to solve the linear system");

#else // PRECOND_TEST

  // Pass right hand side to the linear solver
  linsol_.setInput(NV_DATA_S(rvec),1);
  
  // Solve the system
  linsol_.evaluate();
  
  // Get the result
  linsol_.getOutput(NV_DATA_S(zvec));

#endif // PRECOND_TEST

  // Log time duration
  time2 = clock();
  t_psolve += double(time2-time1)/CLOCKS_PER_SEC;

}

void IdasInternal::psetup(double t, N_Vector yy, N_Vector yp, N_Vector rr, double cj, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3){
  // Get time
  time1 = clock();

  // Pass input to the jacobian function
  jacx_.setInput(t,JAC_T);
  jacx_.setInput(NV_DATA_S(yy),JAC_Y);
  jacx_.setInput(NV_DATA_S(yp),JAC_YDOT);
  jacx_.setInput(input(INTEGRATOR_P).data(),JAC_P);
  jacx_.setInput(cj,JAC_CJ);

  // Evaluate jacobian
  jacx_.evaluate();
  
#ifdef PRECOND_TEST
  // DENSE PRECONDITIONER: REMOVE
  
  // Get the result
  jacx_.output(JAC_J).get(pc_,DENSE);

  // Log time duration
  time2 = clock();
  t_psetup_jac += double(time2-time1)/CLOCKS_PER_SEC;

  // Factorize the matrix
  int info = -100;
  
  if(use_qr){
    int lwork = work_.size();
    dgeqrf_(&ny_, &ny_, &pc_[0], &ny_, &tau_[0], &work_[0], &lwork, &info);
  } else {
    dgetrf_(&ny_, &ny_, &pc_[0], &ny_, &ipiv_[0], &info);
  }
    
  if(info != 0)
    throw CasadiException("dgetrf_ failed to factorize the jacobian");
#else // PRECOND_TEST


  // Pass non-zero elements to the linear solver
  linsol_.setInput(jacx_.getOutputData(),0);

  #endif // PRECOND_TEST

  // Log time duration
  time1 = clock();
  t_psetup_fac += double(time1-time2)/CLOCKS_PER_SEC;

}

int IdasInternal::lsetup_wrapper(IDAMem IDA_mem, N_Vector yyp, N_Vector ypp, N_Vector resp, N_Vector vtemp1, N_Vector vtemp2, N_Vector vtemp3){
 try{
    IdasInternal *this_ = (IdasInternal*)(IDA_mem->ida_lmem);
    assert(this_);
    this_->lsetup(IDA_mem,yyp,ypp,resp,vtemp1,vtemp2,vtemp3);
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
  } catch(exception& e){
    cerr << "lsolve failed: " << e.what() << endl;;
    return -1;
  }
}

void IdasInternal::lsetup(IDAMem IDA_mem, N_Vector yyp, N_Vector ypp, N_Vector resp, N_Vector vtemp1, N_Vector vtemp2, N_Vector vtemp3){
  // Current time
  double t = IDA_mem->ida_tn;

  // Multiple of df_dydot to be added to the matrix
  double cj = IDA_mem->ida_cj;

  // Pass to the jacobian function
  jacx_.setInput(t,JAC_T);
  jacx_.setInput(NV_DATA_S(yyp),JAC_Y);
  jacx_.setInput(NV_DATA_S(ypp),JAC_YDOT);
  jacx_.setInput(input(INTEGRATOR_P).data(),JAC_P);
  jacx_.setInput(cj,JAC_CJ);
  
  // Evaluate the Jacobian function
  jacx_.evaluate();
  
  // Pass non-zero elements to the linear solver
  linsol_.setInput(jacx_.getOutputData(),0);
  
  // Tell the solver to refactorize at next call
  refactor_ = true;  
}

void IdasInternal::lsolve(IDAMem IDA_mem, N_Vector b, N_Vector weight, N_Vector ycur, N_Vector ypcur, N_Vector rescur){

  // Pass right hand side to the linear solver
  linsol_.setInput(NV_DATA_S(b),1);
  
  // Solve the system
  linsol_.evaluate();
  
  // Get the result
  linsol_.getOutput(NV_DATA_S(b));
  
  // Matrix refactorized
  refactor_ = false;
}

void IdasInternal::initUserDefinedLinearSolver(){
  // Make sure that a Jacobian has been provided
  if(jacx_.isNull()) throw CasadiException("IdasInternal::initUserDefinedLinearSolver(): No Jacobian has been provided.");

  // Make sure that a linear solver has been providided
  if(linsol_.isNull()) throw CasadiException("IdasInternal::initUserDefinedLinearSolver(): No user defined linear solver has been provided.");

  // Initialize the linear solver
  linsol_.init();

  //  Set fields in the IDA memory
  IDAMem IDA_mem = IDAMem(mem_);
  IDA_mem->ida_lmem   = this;
  IDA_mem->ida_lsetup = lsetup_wrapper;
  IDA_mem->ida_lsolve = lsolve_wrapper;
                                 
  // Factorize at the first call
  refactor_ = true;
}



} // namespace Sundials
} // namespace CasADi

