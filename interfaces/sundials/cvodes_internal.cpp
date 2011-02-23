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

#include "cvodes_internal.hpp"
#include "casadi/fx/sx_function_internal.hpp"
#include "casadi/stl_vector_tools.hpp"
#include "casadi/sx/sx_tools.hpp"
#include "casadi/fx/linear_solver_internal.hpp"

using namespace std;
namespace CasADi{
namespace Sundials{

CVodesInternal* CVodesInternal::clone() const{
  // Return a deep copy
  FX f = f_;   f.makeUnique();
  FX q = q_;   q.makeUnique();
  CVodesInternal* node = new CVodesInternal(f,q);
  node->setOption(dictionary());
  node->setLinearSolver(linsol_,M_);
  node->linsol_.makeUnique();
  node->M_.makeUnique();
  if(isInit())
    node->init();
  return node;
}
  
CVodesInternal::CVodesInternal(const FX& f, const FX& q) : f_(f), q_(q){
  addOption("linear_multistep_method",     OT_STRING,  "bdf"); // "bdf" or "adams"
  addOption("nonlinear_solver_iteration",  OT_STRING,  "newton"); // "newton" or "functional"
  addOption("fsens_all_at_once",           OT_BOOLEAN,true); // calculate all right hand sides of the sensitivity equations at once

  mem_ = 0;

  y0_ = y_ = 0;
  yQ0_ = yQ_ = 0;

  is_init = false;
  isInitAdj_ = false;


}

CVodesInternal::~CVodesInternal(){ 
  if(mem_) CVodeFree(&mem_);

  // ODE integration
  if(y0_) N_VDestroy_Serial(y0_);
  if(y_) N_VDestroy_Serial(y_);
  if(yQ0_) N_VDestroy_Serial(yQ0_);
  if(yQ_) N_VDestroy_Serial(yQ_);
  
  // Forward problem
  for(vector<N_Vector>::iterator it=yS0_.begin(); it != yS0_.end(); ++it)   if(*it) N_VDestroy_Serial(*it);
  for(vector<N_Vector>::iterator it=yS_.begin(); it != yS_.end(); ++it)     if(*it) N_VDestroy_Serial(*it);
  for(vector<N_Vector>::iterator it=yQS0_.begin(); it != yQS0_.end(); ++it) if(*it) N_VDestroy_Serial(*it);
  for(vector<N_Vector>::iterator it=yQS_.begin(); it != yQS_.end(); ++it)   if(*it) N_VDestroy_Serial(*it);
  
  // Adjoint problem
  for(vector<N_Vector>::iterator it=yB0_.begin(); it != yB0_.end(); ++it)   if(*it) N_VDestroy_Serial(*it);
  for(vector<N_Vector>::iterator it=yB_.begin(); it != yB_.end(); ++it)     if(*it) N_VDestroy_Serial(*it);
//  for(vector<N_Vector>::iterator it=yQB0_.begin(); it != yQB0_.end(); ++it) if(*it) N_VDestroy_Serial(*it);
  for(vector<N_Vector>::iterator it=yQB_.begin(); it != yQB_.end(); ++it)   if(*it) N_VDestroy_Serial(*it);
}

void CVodesInternal::init(){
  // Init ODE rhs function and quadrature functions
  f_.init();
  casadi_assert(f_.getNumInputs()==ODE_NUM_IN);
  casadi_assert(f_.getNumOutputs()==ODE_NUM_OUT);
  if(!q_.isNull()){
    q_.init();
    casadi_assert(q_.getNumInputs()==ODE_NUM_IN);
    casadi_assert(q_.getNumOutputs()==ODE_NUM_OUT);
  }

  // Number of states
  int nx = f_.output(INTEGRATOR_XF).numel();

  // Add quadratures, if any
  if(!q_.isNull()) nx += q_.output().numel();

  // Number of parameters
  int np = f_.input(ODE_P).numel();
  
  // Set dimensions
  setDimensions(nx,np,0);
  
  ny_ = f_.output().numel();
  nq_ = q_.isNull() ? 0 : q_.output().numel();
  
  IntegratorInternal::init();
  
  // Try to generate a jacobian of none provided
  if(!linsol_.isNull() && M_.isNull()){
    SXFunction f = shared_cast<SXFunction>(f_);
    if(!f.isNull()){
      // Get the Jacobian in the Newton iteration
      SX gamma("gamma");
      SXMatrix jac = eyeSX(ny_) - gamma * f.jac(ODE_Y,ODE_RHS);
      
      // Jacobian function
      vector<vector<SX> > jac_in(Sundials::M_NUM_IN);
      jac_in[M_T] = f->inputv.at(ODE_T);
      jac_in[M_Y] = f->inputv.at(ODE_Y);
      jac_in[M_P] = f->inputv.at(ODE_P);
      jac_in[M_GAMMA] = vector<SX>(1,gamma);
      SXFunction M(jac_in,jac);
      
      // Pass sparsity to linear solver
      linsol_.setSparsity(jac.rowind(),jac.col());
      
      // Save function
      M_ = M;
    }
  }
  
  if(!M_.isNull()) M_.init();
  if(!linsol_.isNull()) linsol_.init();

  // Get the number of forward and adjoint directions
  nfdir_f_ = f_.getOption("number_of_fwd_dir").toInt();
  nadir_f_ = f_.getOption("number_of_adj_dir").toInt();
  nfdir_q_ = q_.isNull() ? 0 : q_.getOption("number_of_fwd_dir").toInt();
  nadir_q_ = q_.isNull() ? 0 : q_.getOption("number_of_adj_dir").toInt();

  // Quick return if already initialized
  if(is_init){
    log("CVodesInternal::init","end, Idas already initialized");
    return;
  }

  // Sundials return flag
  int flag;

  if(getOption("linear_multistep_method")=="adams")  lmm_ = CV_ADAMS;
  else if(getOption("linear_multistep_method")=="bdf") lmm_ = CV_BDF;
  else throw CasadiException("Unknown linear multistep method");

  if(getOption("nonlinear_solver_iteration")=="newton") iter_ = CV_NEWTON;
  else if(getOption("nonlinear_solver_iteration")=="functional") iter_ = CV_FUNCTIONAL;
  else throw CasadiException("Unknown nonlinear solver iteration");

  // Create CVodes memory block
  mem_ = CVodeCreate(lmm_,iter_);
  if(mem_==0) throw CasadiException("CVodeCreate: Creation failed");

  // Allocate n-vectors for ivp
  y0_ = N_VMake_Serial(ny_,&input(INTEGRATOR_X0)[0]);
  y_ = N_VMake_Serial(ny_,&output(INTEGRATOR_XF)[0]);

  // Set error handler function
  flag = CVodeSetErrHandlerFn(mem_, ehfun_wrapper, this);
  if(flag != CV_SUCCESS) cvodes_error("CVodeSetErrHandlerFn",flag);

  // Initialize CVodes
  double t0 = 0;
  flag = CVodeInit(mem_, rhs_wrapper, t0, y0_);
  if(flag!=CV_SUCCESS) cvodes_error("CVodeInit",flag);

  // Set tolerances
  flag = CVodeSStolerances(mem_, reltol_, abstol_);
  if(flag!=CV_SUCCESS) cvodes_error("CVodeInit",flag);

  // Maximum number of steps
  CVodeSetMaxNumSteps(mem_, getOption("max_num_steps").toInt());
  if(flag != CV_SUCCESS) cvodes_error("CVodeSetMaxNumSteps",flag);
  
  // attach a linear solver
  if(getOption("linear_solver")=="dense"){
    // Dense jacobian
    flag = CVDense(mem_, ny_);
    if(flag!=CV_SUCCESS) cvodes_error("CVDense",flag);
    if(exact_jacobian_){
      // Create jacobian if it does not exist
      if(jac_f_.isNull()) jac_f_ = f_.jacobian(ODE_Y,ODE_RHS);
      jac_f_.init();
      
      // Pass to CVodes
      flag = CVDlsSetDenseJacFn(mem_, djac_wrapper);
      if(flag!=CV_SUCCESS) cvodes_error("CVDlsSetDenseJacFn",flag);
    }
  } else if(getOption("linear_solver")=="banded") {
    // Banded jacobian
    flag = CVBand(mem_, ny_, getOption("upper_bandwidth").toInt(), getOption("lower_bandwidth").toInt());
    if(flag!=CV_SUCCESS) cvodes_error("CVBand",flag);
    if(exact_jacobian_){
      flag = CVDlsSetBandJacFn(mem_, bjac_wrapper);
      if(flag!=CV_SUCCESS) cvodes_error("CVDlsSetBandJacFn",flag);
    }
  } else if(getOption("linear_solver")=="iterative") {
    // Sparse (iterative) solver  

    // Preconditioning type
    int pretype;
    if(getOption("pretype")=="none")               pretype = PREC_NONE;
    else if(getOption("pretype")=="left")          pretype = PREC_LEFT;
    else if(getOption("pretype")=="right")         pretype = PREC_RIGHT;
    else if(getOption("pretype")=="both")          pretype = PREC_BOTH;
    else                                           throw CasadiException("Unknown preconditioning type");

    // Max dimension of the Krylov space
    int maxl = getOption("max_krylov").toInt();

    // Attach the sparse solver  
    if(getOption("iterative_solver")=="gmres"){
      flag = CVSpgmr(mem_, pretype, maxl);
      if(flag!=CV_SUCCESS) cvodes_error("CVBand",flag);      
    } else if(getOption("iterative_solver")=="bcgstab") {
      flag = CVSpbcg(mem_, pretype, maxl);
      if(flag!=CV_SUCCESS) cvodes_error("CVSpbcg",flag);
    } else if(getOption("iterative_solver")=="tfqmr") {
      flag = CVSptfqmr(mem_, pretype, maxl);
      if(flag!=CV_SUCCESS) cvodes_error("CVSptfqmr",flag);
    } else throw CasadiException("CVODES: Unknown sparse solver");
      
    // Attach functions for jacobian information
    if(exact_jacobian_){
      flag = CVSpilsSetJacTimesVecFn(mem_, jtimes_wrapper);
      if(flag!=CV_SUCCESS) cvodes_error("CVSpilsSetJacTimesVecFn",flag);      
    }
    
    // Add a preconditioner
    if(getOption("use_preconditioner")==true){
      // Make sure that a Jacobian has been provided
      if(M_.isNull()) throw CasadiException("CVodesInternal::init(): No Jacobian has been provided.");

      // Make sure that a linear solver has been providided
      if(linsol_.isNull()) throw CasadiException("CVodesInternal::init(): No user defined linear solver has been provided.");

      // Pass to IDA
      flag = CVSpilsSetPreconditioner(mem_, psetup_wrapper, psolve_wrapper);
      if(flag != CV_SUCCESS) cvodes_error("CVSpilsSetPreconditioner",flag);
    }    
  } else if(getOption("linear_solver")=="user_defined") {
    initUserDefinedLinearSolver();
  } else throw CasadiException("Unknown linear solver ");

  // Set user data
  flag = CVodeSetUserData(mem_,this);
  if(flag!=CV_SUCCESS) cvodes_error("CVodeSetUserData",flag);

  // Quadrature equations
  if(nq_>0){
    // Allocate n-vectors for quadratures
    yQ0_ = N_VMake_Serial(nq_,&input(INTEGRATOR_X0)[ny_]);
    yQ_ = N_VMake_Serial(nq_,&output(INTEGRATOR_XF)[ny_]);

    // Initialize quadratures in CVodes
    flag = CVodeQuadInit(mem_, rhsQ_wrapper, yQ0_);
    if(flag != CV_SUCCESS) cvodes_error("CVodeQuadInit",flag);
    
    // Should the quadrature errors be used for step size control?
    if(getOption("quad_err_con").toInt()){
      flag = CVodeSetQuadErrCon(mem_, true);
      if(flag != CV_SUCCESS) cvodes_error("IDASetQuadErrCon",flag);
      
      // Quadrature error tolerances
      flag = CVodeQuadSStolerances(mem_, reltol_, abstol_); // TODO: vector absolute tolerances
      if(flag != CV_SUCCESS) cvodes_error("CVodeQuadSStolerances",flag);
    }
  }
  
    // Forward sensitivity problem
    if(nfdir_>0){
      // Allocate n-vectors
      yS0_.resize(nfdir_,0);
      yS_.resize(nfdir_,0);
      for(int i=0; i<nfdir_; ++i){
        yS0_[i] = N_VMake_Serial(ny_,&fwdSeed(INTEGRATOR_X0,i)[0]);
        yS_[i] = N_VMake_Serial(ny_,&fwdSens(INTEGRATOR_XF,i)[0]);
      }

      // Allocate n-vectors for quadratures
      if(nq_>0){
        yQS0_.resize(nfdir_,0);
        yQS_.resize(nfdir_,0);
        for(int i=0; i<nfdir_; ++i){
          yQS0_[i] = N_VMake_Serial(nq_,&fwdSeed(INTEGRATOR_X0,i)[ny_]);
          yQS_[i] = N_VMake_Serial(nq_,&fwdSens(INTEGRATOR_XF,i)[ny_]);
        }
      }
      
    // Calculate all forward sensitivity right hand sides at once?
    bool all_at_once = getOption("fsens_all_at_once").toInt();
      
    // Get the sensitivity method
    if(getOption("sensitivity_method")=="simultaneous") ism_ = CV_SIMULTANEOUS;
    else if(getOption("sensitivity_method")=="staggered") ism_ = all_at_once ? CV_STAGGERED : CV_STAGGERED1;
    else throw CasadiException("CVodes: Unknown sensitivity method");
  
    // Initialize forward sensitivities
    if(finite_difference_fsens_){
      // Use finite differences to calculate the residual in the forward sensitivity equations
      if(all_at_once){
        flag = CVodeSensInit(mem_,nfdir_,ism_,0,&yS0_[0]);
        if(flag != CV_SUCCESS) cvodes_error("CVodeSensInit",flag);
      } else {
        flag = CVodeSensInit1(mem_,nfdir_,ism_,0,&yS0_[0]);
        if(flag != CV_SUCCESS) cvodes_error("CVodeSensInit1",flag);
      }
      
      // Pass pointer to parameters
      flag = CVodeSetSensParams(mem_,&input(INTEGRATOR_P)[0],0,0);
      if(flag != CV_SUCCESS) cvodes_error("CVodeSetSensParams",flag);

      //  CVodeSetSensDQMethod

    } else {
      if(all_at_once){
        // Use AD to calculate the residual in the forward sensitivity equations
        flag = CVodeSensInit(mem_,nfdir_,ism_,rhsS_wrapper,&yS0_[0]);
        if(flag != CV_SUCCESS) cvodes_error("CVodeSensInit",flag);
      } else {
        flag = CVodeSensInit1(mem_,nfdir_,ism_,rhsS1_wrapper,&yS0_[0]);
        if(flag != CV_SUCCESS) cvodes_error("CVodeSensInit",flag);
      }
    }
    
    // Set tolerances
    vector<double> fsens_abstol(nfdir_,fsens_abstol_);
    
    flag = CVodeSensSStolerances(mem_,fsens_reltol_,&fsens_abstol[0]);
    if(flag != CV_SUCCESS) cvodes_error("CVodeSensSStolerances",flag);
    
    // Set optional inputs
    int errconS = getOption("fsens_err_con").toInt();
    flag = CVodeSetSensErrCon(mem_, errconS);
    if(flag != CV_SUCCESS) cvodes_error("CVodeSetSensErrCon",flag);
    
    // Quadrature equations
    if(nq_>0){
      flag = CVodeQuadSensInit(mem_, rhsQS_wrapper, &yQS0_[0]);
      if(flag != CV_SUCCESS) cvodes_error("CVodeQuadSensInit",flag);

      // Set tolerances
      flag = CVodeQuadSensSStolerances(mem_,fsens_reltol_,&fsens_abstol[0]);
      if(flag != CV_SUCCESS) cvodes_error("CVodeQuadSensSStolerances",flag);
    }
  } // enable fsens
    
  // Adjoint sensitivity problem
  whichB_.resize(nadir_);

  // Allocate n-vectors
  yB0_.resize(nadir_,0);
  yB_.resize(nadir_,0);
  for(int i=0; i<nadir_; ++i){
    yB0_[i] = N_VMake_Serial(ny_,&adjSeed(INTEGRATOR_XF,i)[0]);
    yB_[i] = N_VMake_Serial(ny_,&adjSens(INTEGRATOR_X0,i)[0]);
  }

  // Allocate n-vectors for quadratures
  yQB_.resize(nadir_,0);
  for(int i=0; i<nadir_; ++i){
    casadi_assert(adjSens(INTEGRATOR_P,i).size()==np_);
    yQB_[i] = N_VMake_Serial(np_,&adjSens(INTEGRATOR_P,i)[0]);
  }
  
  if(nadir_>0){
    // Get the number of steos per checkpoint
    int Nd = getOption("steps_per_checkpoint").toInt();

    // Get the interpolation type
    int interpType;
    if(getOption("interpolation_type")=="hermite")
      interpType = CV_HERMITE;
    else if(getOption("interpolation_type")=="polynomial")
      interpType = CV_POLYNOMIAL;
    else throw CasadiException("\"interpolation_type\" must be \"hermite\" or \"polynomial\"");
      
    // Initialize adjoint sensitivities
    flag = CVodeAdjInit(mem_, Nd, interpType);
    if(flag != CV_SUCCESS) cvodes_error("CVodeAdjInit",flag);
  }
        
  is_init = true;
  isInitAdj_ = false;
}


void CVodesInternal::initAdj(){
  int flag;

  for(int dir=0; dir<nadir_; ++dir){

      // Create backward problem (use the same lmm and iter)
      flag = CVodeCreateB(mem_, lmm_, iter_, &whichB_[dir]);
      if(flag != CV_SUCCESS) cvodes_error("CVodeCreateB",flag);
      
      // Initialize the backward problem
      double tB0 = input(INTEGRATOR_TF)[0];
      flag = CVodeInitB(mem_, whichB_[dir], rhsB_wrapper, tB0, yB0_[dir]);
      if(flag != CV_SUCCESS) cvodes_error("CVodeInitB",flag);

      // Set tolerances
      flag = CVodeSStolerancesB(mem_, whichB_[dir], asens_reltol_, asens_abstol_);
      if(flag!=CV_SUCCESS) cvodes_error("CVodeSStolerancesB",flag);

      // User data
      flag = CVodeSetUserDataB(mem_, whichB_[dir], this);
      if(flag != CV_SUCCESS) cvodes_error("CVodeSetUserDataB",flag);

      // attach linear solver to backward problem
      if(getOption("asens_linear_solver")=="dense"){
      // Dense jacobian
      flag = CVDenseB(mem_, whichB_[dir], ny_);
      if(flag!=CV_SUCCESS) cvodes_error("CVDenseB",flag);
      } else if(getOption("asens_linear_solver")=="banded") {
      // Banded jacobian
      flag = CVBandB(mem_, whichB_[dir], ny_, getOption("asens_upper_bandwidth").toInt(), getOption("asens_lower_bandwidth").toInt());
      if(flag!=CV_SUCCESS) cvodes_error("CVBandB",flag);
    } else if(getOption("asens_linear_solver")=="iterative") {
      // Sparse jacobian
      // Preconditioning type
      int pretype;
      if(getOption("asens_pretype")=="none")               pretype = PREC_NONE;
      else if(getOption("asens_pretype")=="left")          pretype = PREC_LEFT;
      else if(getOption("asens_pretype")=="right")         pretype = PREC_RIGHT;
      else if(getOption("asens_pretype")=="both")          pretype = PREC_BOTH;
      else                                            throw CasadiException("Unknown preconditioning type for backward problem");

      // Attach the sparse solver  
      int maxl = getOption("asens_max_krylov").toInt();
      if(getOption("asens_iterative_solver")=="gmres"){
        flag = CVSpgmrB(mem_, whichB_[dir], pretype, maxl);
        if(flag!=CV_SUCCESS) cvodes_error("CVSpgmrB",flag);      
      } else if(getOption("asens_iterative_solver")=="bcgstab") {
        flag = CVSpbcgB(mem_, whichB_[dir], pretype, maxl);
        if(flag!=CV_SUCCESS) cvodes_error("CVSpbcgB",flag);
      } else if(getOption("asens_iterative_solver")=="tfqmr") {
        flag = CVSptfqmrB(mem_, whichB_[dir], pretype, maxl);
        if(flag!=CV_SUCCESS) cvodes_error("CVSptfqmrB",flag);
      } else throw CasadiException("Unknown sparse solver for backward problem");
        
      // Attach functions for jacobian information
    } else throw CasadiException("Unknown linear solver for backward problem");

    // Quadratures for the adjoint problem

    N_VConst(0.0, yQB_[dir]);
    flag = CVodeQuadInitB(mem_,whichB_[dir],rhsQB_wrapper,yQB_[dir]);
    if(flag!=CV_SUCCESS) cvodes_error("CVodeQuadInitB",flag);
    
    if(getOption("quad_err_con").toInt()){
      flag = CVodeSetQuadErrConB(mem_, whichB_[dir],true);
      if(flag != CV_SUCCESS) cvodes_error("CVodeSetQuadErrConB",flag);
        
      flag = CVodeQuadSStolerancesB(mem_, whichB_[dir], asens_reltol_, asens_abstol_);
      if(flag != CV_SUCCESS) cvodes_error("CVodeQuadSStolerancesB",flag);
    }
  }
  
  // Mark initialized
  isInitAdj_ = true;
}

void CVodesInternal::rhs(double t, const double* y, double* ydot){
  // Get time
  time1 = clock();

  // Pass input
  f_.setInput(t,ODE_T);
  f_.setInput(y,ODE_Y);
  f_.setInput(input(INTEGRATOR_P),ODE_P);

    // Evaluate
  f_.evaluate();
    
    // Get results
  f_.getOutput(ydot);

  // Log time
  time2 = clock();
  t_res += double(time2-time1)/CLOCKS_PER_SEC;
}

int CVodesInternal::rhs_wrapper(double t, N_Vector y, N_Vector ydot, void *user_data){
try{
    casadi_assert(user_data);
    CVodesInternal *this_ = (CVodesInternal*)user_data;
    this_->rhs(t,NV_DATA_S(y),NV_DATA_S(ydot));
    return 0;
  } catch(exception& e){
    cerr << "rhs failed: " << e.what() << endl;
    return 1;
  }
}
  
void CVodesInternal::reset(int fsens_order, int asens_order){
  if(monitored("CVodesInternal::reset")){
    cout << "initial state: " << endl;
    cout << "p = " << input(INTEGRATOR_P) << endl;
    cout << "x0 = " << input(INTEGRATOR_X0) << endl;
  }

  // Reset timers
  t_res = t_fres = t_jac = t_lsolve = t_lsetup_jac = t_lsetup_fac = 0;
  
  fsens_order_ = fsens_order;
  asens_order_ = asens_order;
  
  // Get the time horizon
  double t0 = input(INTEGRATOR_T0)[0];
  double tf = input(INTEGRATOR_TF)[0];
  t_ = t0;

  // Re-initialize
  int flag = CVodeReInit(mem_, t0, y0_);
  if(flag!=CV_SUCCESS) cvodes_error("CVodeReInit",flag);
  
  // Re-initialize quadratures
  if(nq_>0){
    flag = CVodeQuadReInit(mem_, yQ0_);
    if(flag != CV_SUCCESS) cvodes_error("CVodeQuadReInit",flag);
  }
  
  // Re-initialize sensitivities
  if(fsens_order_>0){
    flag = CVodeSensReInit(mem_,ism_,&yS0_[0]);
    if(flag != CV_SUCCESS) cvodes_error("CVodeSensReInit",flag);
    
    if(nq_>0){
      flag = CVodeQuadSensReInit(mem_, &yQS0_[0]);
      if(flag != CV_SUCCESS) cvodes_error("CVodeQuadSensReInit",flag);
    }
  } else {
    // Turn of sensitivities
    flag = CVodeSensToggleOff(mem_);
    if(flag != CV_SUCCESS) cvodes_error("CVodeSensToggleOff",flag);
  }
}

void CVodesInternal::integrate(double t_out){
  int flag;
  
  // tolerance
  double ttol = 1e-9;
  if(fabs(t_-t_out)<ttol){
    copy(input(INTEGRATOR_X0).begin(),input(INTEGRATOR_X0).end(),output(INTEGRATOR_XF).begin());
    if(fsens_order_>0){
      for(int i=0; i<nfdir_; ++i){
        copy(fwdSeed(INTEGRATOR_X0,i).begin(),fwdSeed(INTEGRATOR_X0,i).end(),fwdSens(INTEGRATOR_XF,i).begin());
      }
    }
    return;
  }
  if(asens_order_>0){
    int ncheck; // number of checkpoints stored so far
    flag = CVodeF(mem_, t_out, y_, &t_, CV_NORMAL,&ncheck);
    if(flag!=CV_SUCCESS && flag!=CV_TSTOP_RETURN) cvodes_error("CVodeF",flag);
    
  } else {
    flag = CVode(mem_, t_out, y_, &t_, CV_NORMAL);
    if(flag!=CV_SUCCESS && flag!=CV_TSTOP_RETURN) cvodes_error("CVode",flag);
  }
  
  if(nq_>0){
    double tret;
    flag = CVodeGetQuad(mem_, &tret, yQ_);
    if(flag!=CV_SUCCESS) cvodes_error("CVodeGetQuad",flag);    
  }
  
  if(fsens_order_>0){
    // Get the sensitivities
    flag = CVodeGetSens(mem_, &t_, &yS_[0]);
    if(flag != CV_SUCCESS) cvodes_error("CVodeGetSens",flag);
    
    if(nq_>0){
      double tret;
      flag = CVodeGetQuadSens(mem_, &tret, &yQS_[0]);
      if(flag != CV_SUCCESS) cvodes_error("CVodeGetQuadSens",flag);
    }
  }
}

void CVodesInternal::resetAdj(){
  double tf = input(INTEGRATOR_TF)[0];
  int flag;
  
  if(isInitAdj_){
    for(int dir=0; dir<nadir_; ++dir){
      flag = CVodeReInitB(mem_, whichB_[dir], tf, yB0_[dir]);
      if(flag != CV_SUCCESS) cvodes_error("CVodeReInitB",flag);

      N_VConst(0.0,yQB_.at(dir));
      flag = CVodeQuadReInitB(mem_,whichB_[dir],yQB_[dir]);
      if(flag!=CV_SUCCESS) cvodes_error("CVodeQuadReInitB",flag);
    }
  } else {
    // Initialize the adjoint integration
    initAdj();
  }
}

void CVodesInternal::integrateAdj(double t_out){
  int flag;
  
  // Integrate backwards to t_out
  flag = CVodeB(mem_, t_out, CV_NORMAL);
  if(flag<CV_SUCCESS) cvodes_error("CVodeB",flag);

  // Get the sensitivities
  double tret;
  for(int dir=0; dir<nadir_; ++dir){
    flag = CVodeGetB(mem_, whichB_[dir], &tret, yB_[dir]);
    if(flag!=CV_SUCCESS) cvodes_error("CVodeGetB",flag);

    flag = CVodeGetQuadB(mem_, whichB_[dir], &tret, yQB_[dir]);
    if(flag!=CV_SUCCESS) cvodes_error("CVodeGetQuadB",flag);
        
  }
}

void CVodesInternal::printStats(std::ostream &stream) const{

    long nsteps, nfevals, nlinsetups, netfails;
    int qlast, qcur;
    double hinused, hlast, hcur, tcur;
    
    int flag;

    flag = CVodeGetIntegratorStats(mem_, &nsteps, &nfevals,&nlinsetups, &netfails, &qlast, &qcur,&hinused, &hlast, &hcur, &tcur);
    if(flag!=CV_SUCCESS) cvodes_error("CVodeGetIntegratorStats",flag);

    stream << "number of steps taken by CVODES: " << nsteps << std::endl;
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

    stream << "Time spent in the ODE residual: " << t_res << " s." << endl;
    stream << "Time spent in the forward sensitivity residual: " << t_fres << " s." << endl;
    stream << "Time spent in the jacobian function or jacobian times vector function: " << t_jac << " s." << endl;
    stream << "Time spent in the linear solver solve function: " << t_lsolve << " s." << endl;
    stream << "Time spent to generate the jacobian in the linear solver setup function: " << t_lsetup_jac << " s." << endl;
    stream << "Time spent to factorize the jacobian in the linear solver setup function: " << t_lsetup_fac << " s." << endl;
    stream << std::endl;

    
    
#if 0
  // Quadrature
  if(ops.quadrature && ocp.hasFunction(LTERM)){
      long nfQevals, nQetfails;
      flag = CVodeGetQuadStats(cvode_mem_[k], &nfQevals, &nQetfails);  
      if(flag != CV_SUCCESS) throw "Error in CVodeGetQuadStats";

      stream << "Quadrature: " << std::endl;
      stream << "number of calls made to the user's quadrature right-hand side function: " << nfQevals << std::endl;
      stream << "number of local error test failures due to quadrature variables: " <<  nQetfails << std::endl;
      stream << std::endl;
}
#endif
}
  
map<int,string> CVodesInternal::calc_flagmap(){
  map<int,string> f;
  f[CV_SUCCESS] = "CV_SUCCESS";
  f[CV_TSTOP_RETURN] = "CV_TSTOP_RETURN";
  f[CV_ROOT_RETURN] = "CV_ROOT_RETURN";
  f[CV_WARNING] = "CV_WARNING";
  f[CV_WARNING] = "CV_WARNING";
  f[CV_TOO_MUCH_WORK] = "CV_TOO_MUCH_WORK";
  f[CV_TOO_MUCH_ACC] = "CV_TOO_MUCH_ACC";
  f[CV_ERR_FAILURE] = "CV_ERR_FAILURE";
  f[CV_CONV_FAILURE] = "CV_CONV_FAILURE";
  f[CV_LINIT_FAIL] = "CV_LINIT_FAIL";
  f[CV_LSETUP_FAIL] = "CV_LSETUP_FAIL";
  f[CV_LSOLVE_FAIL] = "CV_LSOLVE_FAIL";
  f[CV_RHSFUNC_FAIL] = "CV_RHSFUNC_FAIL";
  f[CV_FIRST_RHSFUNC_ERR] = "CV_FIRST_RHSFUNC_ERR";
  f[CV_UNREC_RHSFUNC_ERR] = "CV_UNREC_RHSFUNC_ERR";
  f[CV_RTFUNC_FAIL] = "CV_RTFUNC_FAIL";
  f[CV_MEM_FAIL] = "CV_MEM_FAIL";
  f[CV_ILL_INPUT] = "CV_ILL_INPUT";
  f[CV_NO_MALLOC] = "CV_NO_MALLOC";
  f[CV_BAD_K] = "CV_BAD_K";
  f[CV_BAD_T] = "CV_BAD_T";
  f[CV_BAD_DKY] = "CV_BAD_DKY";
  f[CV_TOO_CLOSE] = "CV_TOO_CLOSE";
  f[CV_QRHSFUNC_FAIL] = "CV_QRHSFUNC_FAIL";
  f[CV_FIRST_QRHSFUNC_ERR] = "CV_FIRST_QRHSFUNC_ERR";
  f[CV_REPTD_QRHSFUNC_ERR] = "CV_REPTD_QRHSFUNC_ERR";
  f[CV_UNREC_QRHSFUNC_ERR] = "CV_UNREC_QRHSFUNC_ERR";
  f[CV_NO_SENS ] = "CV_NO_SENS ";
  f[CV_SRHSFUNC_FAIL] = "CV_SRHSFUNC_FAIL";
  return f;
}
  
map<int,string> CVodesInternal::flagmap = CVodesInternal::calc_flagmap();

void CVodesInternal::cvodes_error(const string& module, int flag){
  // Find the error
  map<int,string>::const_iterator it = flagmap.find(flag);
  
  stringstream ss;
  if(it == flagmap.end()){
    ss << "Unknown error (" << flag << ") from module \"" << module << "\".";
  } else {
    ss << "Module \"" << module << "\" returned flag \"" << it->second << "\".";
  }
  ss << " Consult Cvodes documentation.";
  throw CasadiException(ss.str());
}
  
void CVodesInternal::ehfun_wrapper(int error_code, const char *module, const char *function, char *msg, void *eh_data){
try{
    casadi_assert(eh_data);
    CVodesInternal *this_ = (CVodesInternal*)eh_data;
    this_->ehfun(error_code,module,function,msg);        
  } catch(exception& e){
    cerr << "ehfun failed: " << e.what() << endl;
  }
}
  
void CVodesInternal::ehfun(int error_code, const char *module, const char *function, char *msg){
  cerr << msg << endl;
}

void CVodesInternal::rhsS(int Ns, double t, N_Vector y, N_Vector ydot, N_Vector *yS, N_Vector *ySdot, N_Vector tmp1, N_Vector tmp2){
  casadi_assert(Ns==nfdir_);

  // Record the current cpu time
  time1 = clock();
  
    // Pass input
  f_.setInput(t,ODE_T);
  f_.setInput(NV_DATA_S(y),ODE_Y);
  f_.setInput(input(INTEGRATOR_P),ODE_P);

   // Calculate the forward sensitivities, nfdir_f_ directions at a time
   for(int j=0; j<nfdir_; j += nfdir_f_){
     for(int dir=0; dir<nfdir_f_ && j+dir<nfdir_; ++dir){
       // Pass forward seeds 
       f_.setFwdSeed(0.0,ODE_T,dir);
       f_.setFwdSeed(NV_DATA_S(yS[j+dir]),ODE_Y,dir);
       f_.setFwdSeed(fwdSeed(INTEGRATOR_P,j+dir),ODE_P,dir);
     }

     // Evaluate the AD forward algorithm
     f_.evaluate(1,0);
      
     // Get the output seeds
     for(int dir=0; dir<nfdir_f_ && j+dir<nfdir_; ++dir){
       f_.getFwdSens(NV_DATA_S(ySdot[j+dir]),ODE_RHS,dir);
     }
   }
  
  // Record timings
  time2 = clock();
  t_fres += double(time2-time1)/CLOCKS_PER_SEC;
}

int CVodesInternal::rhsS_wrapper(int Ns, double t, N_Vector y, N_Vector ydot, N_Vector *yS, N_Vector *ySdot, void *user_data, N_Vector tmp1, N_Vector tmp2){
try{
    casadi_assert(user_data);
    CVodesInternal *this_ = (CVodesInternal*)user_data;
    this_->rhsS(Ns,t,y,ydot,yS,ySdot,tmp1,tmp2);
    return 0;
  } catch(exception& e){
    cerr << "fs failed: " << e.what() << endl;
    return 1;
  }
}

void CVodesInternal::rhsS1(int Ns, double t, N_Vector y, N_Vector ydot, int iS, N_Vector yS, N_Vector ySdot, N_Vector tmp1, N_Vector tmp2){
  casadi_assert(Ns==nfdir_);
  
    // Pass input
  f_.setInput(t,ODE_T);
  f_.setInput(NV_DATA_S(y),ODE_Y);
  f_.setInput(input(INTEGRATOR_P),ODE_P);

  // Pass forward seeds
  f_.setFwdSeed(0.0,ODE_T);
  f_.setFwdSeed(NV_DATA_S(yS),ODE_Y);
  f_.setFwdSeed(fwdSeed(INTEGRATOR_P,iS),ODE_P);
    
  // Evaluate the AD forward algorithm
  f_.evaluate(1,0);
  
  // Get the fwd sensitivities
  f_.getFwdSens(NV_DATA_S(ySdot));
}

int CVodesInternal::rhsS1_wrapper(int Ns, double t, N_Vector y, N_Vector ydot, int iS, N_Vector yS, N_Vector ySdot, void *user_data, N_Vector tmp1, N_Vector tmp2){
try{
    casadi_assert(user_data);
    CVodesInternal *this_ = (CVodesInternal*)user_data;
    this_->rhsS1(Ns,t,y,ydot,iS,yS,ySdot,tmp1,tmp2);
    return 0;
  } catch(exception& e){
    cerr << "fs failed: " << e.what() << endl;
    return 1;
  }
}

int CVodesInternal::rhsQ_wrapper(double t, N_Vector yy, N_Vector rhsQ, void *user_data){
try{
    casadi_assert(user_data);
    CVodesInternal *this_ = (CVodesInternal*)user_data;
    this_->rhsQ(t,NV_DATA_S(yy),NV_DATA_S(rhsQ));
    return 0;
  } catch(exception& e){
    cerr << "rhsQ failed: " << e.what() << endl;;
    return 1;
  }
}

void CVodesInternal::rhsQ(double t, const double* yy, double* rhsQ){
// Pass input
  q_.setInput(t,ODE_T);
  q_.setInput(yy,ODE_Y);
  q_.setInput(input(INTEGRATOR_P),ODE_P);

  // Evaluate
  q_.evaluate();
    
  // Get results
  q_.getOutput(rhsQ);
}

void CVodesInternal::rhsQS(int Ns, double t, N_Vector y, N_Vector *yS, N_Vector yQdot, N_Vector *rhsvalQS, N_Vector tmp1, N_Vector tmp2){
  casadi_assert(Ns==nfdir_);
  
  // Pass input
  q_.setInput(t,ODE_T);
  q_.setInput(NV_DATA_S(y),ODE_Y);
  q_.setInput(input(INTEGRATOR_P),ODE_P);

  for(int i=0; i<nfdir_; ++i){
    // Pass forward seeds
    q_.setFwdSeed(0.0,ODE_T);
    q_.setFwdSeed(NV_DATA_S(yS[i]),ODE_Y);
    q_.setFwdSeed(fwdSeed(INTEGRATOR_P,i),ODE_P);

    // Evaluate the AD forward algorithm
    q_.evaluate(1,0);
      
    // Get the forward sensitivities
    q_.getFwdSens(NV_DATA_S(rhsvalQS[i]));
  }
}

int CVodesInternal::rhsQS_wrapper(int Ns, double t, N_Vector y, N_Vector *yS, N_Vector yQdot, N_Vector *rhsvalQS, void *user_data, N_Vector tmp1, N_Vector tmp2){
try{
//    casadi_assert(user_data);
    CVodesInternal *this_ = (CVodesInternal*)user_data;
    if(!this_){
      // SUNDIALS BUG!!!
      for(int i=0; i<Ns; ++i) N_VConst(0.0,rhsvalQS[i]);
      return 0;
    }
    this_->rhsQS(Ns,t,y,yS,yQdot,rhsvalQS,tmp1,tmp2);
    return 0;
  } catch(exception& e){
    cerr << "rhsQS failed: " << e.what() << endl;;
    return 1;
  }
}

void CVodesInternal::rhsB(double t, const double* y, const double *yB, double* yBdot){
  // Pass input
  f_.setInput(t,ODE_T);
  f_.setInput(y,ODE_Y);
  f_.setInput(input(INTEGRATOR_P),ODE_P);

  // Pass adjoint seeds
  f_.setAdjSeed(yB,ODE_RHS);
  
  // Evaluate and tape
  f_.evaluate(0,1);

  // Save to output
  const vector<double>& fres = f_.adjSens(ODE_Y);
  for(int i=0; i<ny_; ++i)
    yBdot[i] = -fres[i];

  // If quadratures are included
  if(nq_>0){
    // Pass input to quadratures
    q_.setInput(t,ODE_T);
    q_.setInput(y,ODE_Y);
    q_.setInput(input(INTEGRATOR_P),ODE_P);

    // Pass adjoint seeds
    q_.setAdjSeed(&adjSeed(INTEGRATOR_XF)[ny_],ODE_RHS);

    // Evaluate
    q_.evaluate(0,1);
    
    // Get the adjoint sensitivities
    const vector<double>& qres = q_.adjSens(ODE_Y);
    
    // Copy to result
    for(int i=0; i<ny_; ++i)
      yBdot[i] -= qres[i];
  }
}

int CVodesInternal::rhsB_wrapper(double t, N_Vector y, N_Vector yB, N_Vector yBdot, void *user_dataB){
try{
    casadi_assert(user_dataB);
    CVodesInternal *this_ = (CVodesInternal*)user_dataB;
    this_->rhsB(t,NV_DATA_S(y),NV_DATA_S(yB),NV_DATA_S(yBdot));
    return 0;
  } catch(exception& e){
    cerr << "rhsB failed: " << e.what() << endl;;
    return 1;
  }
}

int CVodesInternal::rhsQB_wrapper(double t, N_Vector y, N_Vector yB, N_Vector qBdot, void *user_dataB){
  try{
    casadi_assert(user_dataB);
    CVodesInternal *this_ = (CVodesInternal*)user_dataB;
    this_->rhsQB(t,NV_DATA_S(y),NV_DATA_S(yB),NV_DATA_S(qBdot));
    return 0;
  } catch(exception& e){
    cerr << "rhsQB failed: " << e.what() << endl;;
    return 1;
  }
}

void CVodesInternal::rhsQB(double t, const double* y, const double* yB, double* qBdot){
  // Pass input
  f_.setInput(t,ODE_T);
  f_.setInput(y,ODE_Y);
  f_.setInput(input(INTEGRATOR_P),ODE_P);

  // Pass adjoint seeds
  f_.setAdjSeed(yB,ODE_RHS);

  // Evaluate
  f_.evaluate(0,1);

  // Save to output
  f_.getAdjSens(qBdot,ODE_P);
  
  // If quadratures are included
  if(nq_>0){
    // Pass input to quadratures
    q_.setInput(t,ODE_T);
    q_.setInput(y,ODE_Y);
    q_.setInput(input(INTEGRATOR_P),ODE_P);

    // Pass adjoint seeds
    q_.setAdjSeed(&adjSeed(INTEGRATOR_XF)[ny_],ODE_RHS);

    // Evaluate
    q_.evaluate(0,1);
    
    // Get the input seeds
    const vector<double>& qres = q_.adjSens(ODE_P);
    
    // Copy to result
    for(int i=0; i<np_; ++i){
      qBdot[i] += qres[i];
    }
  }
  
  // Negate as we are integrating backwards
  for(int i=0; i<np_; ++i)
    qBdot[i] *= -1;
}

int CVodesInternal::jtimes_wrapper(N_Vector v, N_Vector Jv, double t, N_Vector y, N_Vector fy, void *user_data, N_Vector tmp){
  try{
    casadi_assert(user_data);
    CVodesInternal *this_ = (CVodesInternal*)user_data;
    casadi_assert(this_->f_.fwdSens(ODE_RHS).size() == this_->ny_);
    casadi_assert(NV_LENGTH_S(v) == this_->ny_);
    casadi_assert(NV_LENGTH_S(Jv) == this_->ny_);
    this_->jtimes(NV_DATA_S(v),NV_DATA_S(Jv),t,NV_DATA_S(y),NV_DATA_S(fy),NV_DATA_S(tmp));
    return 0;
  } catch(exception& e){
    cerr << "jtimes failed: " << e.what() << endl;;
    return 1;
  }
}

void CVodesInternal::jtimes(const double *v, double* Jv, double t, const double* y, const double* fy, double* tmp){
  // Get time
  time1 = clock();

  // Pass input
  f_.setInput(t,ODE_T);
  f_.setInput(y,ODE_Y);
  f_.setInput(input(INTEGRATOR_P),ODE_P);

  // Pass input seeds
  f_.setFwdSeed(0.0,ODE_T);
  f_.setFwdSeed(v,ODE_Y);
  fill_n(f_.fwdSeed(ODE_P).begin(),np_,0.0);
  
  // Evaluate
  f_.evaluate(1,0);

  // Get the output seeds
  f_.getFwdSens(Jv,ODE_RHS);
  
  // Log time duration
  time2 = clock();
  t_jac += double(time2-time1)/CLOCKS_PER_SEC;
}

int CVodesInternal::djac_wrapper(int N, double t, N_Vector y, N_Vector fy, DlsMat Jac, void *user_data,N_Vector tmp1, N_Vector tmp2, N_Vector tmp3){
  try{
    casadi_assert(user_data);
    CVodesInternal *this_ = (CVodesInternal*)user_data;
    this_->djac(N, t, y, fy, Jac, tmp1, tmp2, tmp3);
    return 0;
  } catch(exception& e){
    cerr << "djac failed: " << e.what() << endl;;
    return 1;
  }
}

void CVodesInternal::djac(int N, double t, N_Vector y, N_Vector fy, DlsMat Jac, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3){
  // Get time
  time1 = clock();

  // Pass inputs to the jacobian function
  jac_f_.setInput(t,ODE_T);
  jac_f_.setInput(NV_DATA_S(y),ODE_Y);
  jac_f_.setInput(f_.input(ODE_P),ODE_P);

  // Evaluate
  jac_f_.evaluate();
  
  // Get sparsity and non-zero elements
  const vector<int>& rowind = jac_f_.output().rowind();
  const vector<int>& col = jac_f_.output().col();
  const vector<double>& val = jac_f_.output();

  // Loop over rows
  for(int i=0; i<rowind.size()-1; ++i){
    // Loop over non-zero entries
    for(int el=rowind[i]; el<rowind[i+1]; ++el){
      // Get column
      int j = col[el];
      
      // Set the element
      DENSE_ELEM(Jac,i,j) = val[el];
    }
  }
  
  // Log time duration
  time2 = clock();
  t_jac += double(time2-time1)/CLOCKS_PER_SEC;
}

int CVodesInternal::bjac_wrapper(int N, int mupper, int mlower, double t, N_Vector y, N_Vector fy, DlsMat Jac, void *user_data,     
                        N_Vector tmp1, N_Vector tmp2, N_Vector tmp3){
  try{
    casadi_assert(user_data);
    CVodesInternal *this_ = (CVodesInternal*)user_data;
    this_->bjac(N, mupper, mlower, t, y, fy, Jac, tmp1, tmp2, tmp3);
    return 0;
  } catch(exception& e){
    cerr << "bjac failed: " << e.what() << endl;;
    return 1;
  }
}

void CVodesInternal::bjac(int N, int mupper, int mlower, double t, N_Vector y, N_Vector fy, DlsMat Jac, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3){
  // Get time
  time1 = clock();

  // Pass inputs to the jacobian function
  jac_f_.setInput(t,ODE_T);
  jac_f_.setInput(NV_DATA_S(y),ODE_Y);
  jac_f_.setInput(f_.input(ODE_P),ODE_P);

  // Evaluate
  jac_f_.evaluate();
  
  // Get sparsity and non-zero elements
  const vector<int>& rowind = jac_f_.output().rowind();
  const vector<int>& col = jac_f_.output().col();
  const vector<double>& val = jac_f_.output();

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

void CVodesInternal::setStopTime(double tf){
  // Set the stop time of the integration -- don't integrate past this point
  int flag = CVodeSetStopTime(mem_, tf);
  if(flag != CV_SUCCESS) cvodes_error("CVodeSetStopTime",flag);
}

int CVodesInternal::psolve_wrapper(double t, N_Vector y, N_Vector fy, N_Vector r, N_Vector z, double gamma, double delta, int lr, void *user_data, N_Vector tmp){
  try{
    CVodesInternal *this_ = (CVodesInternal*)user_data;
    casadi_assert(this_);
    this_->psolve(t, y, fy, r, z, gamma, delta, lr, tmp);
    return 0;
  } catch(exception& e){
    cerr << "psolve failed: " << e.what() << endl;;
    return 1;
  }
}

int CVodesInternal::psetup_wrapper(double t, N_Vector y, N_Vector fy, booleantype jok, booleantype *jcurPtr, double gamma, void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3){
  try{
    CVodesInternal *this_ = (CVodesInternal*)user_data;
    casadi_assert(this_);
    this_->psetup(t, y, fy, jok, jcurPtr, gamma, tmp1, tmp2, tmp3);
    return 0;
  } catch(exception& e){
    cerr << "psetup failed: " << e.what() << endl;;
    return 1;
  }
}

void CVodesInternal::psolve(double t, N_Vector y, N_Vector fy, N_Vector r, N_Vector z, double gamma, double delta, int lr, N_Vector tmp){
  // Get time
  time1 = clock();

  // Pass right hand side to the linear solver
  linsol_.setInput(NV_DATA_S(r),1);
  
  // Solve the (possibly factorized) system 
  linsol_.solve();
  
  // Get the result
  linsol_.getOutput(NV_DATA_S(z));

  // Log time duration
  time2 = clock();
  t_lsolve += double(time2-time1)/CLOCKS_PER_SEC;
}

void CVodesInternal::psetup(double t, N_Vector y, N_Vector fy, booleantype jok, booleantype *jcurPtr, double gamma, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3){
  // Get time
  time1 = clock();

  // Pass input to the jacobian function
  M_.setInput(t,M_T);
  M_.setInput(NV_DATA_S(y),M_Y);
  M_.setInput(input(INTEGRATOR_P),M_P);
  M_.setInput(gamma,M_GAMMA);

  // Evaluate jacobian
  M_.evaluate();
  
  // Log time duration
  time2 = clock();
  t_lsetup_jac += double(time2-time1)/CLOCKS_PER_SEC;

  // Pass non-zero elements, scaled by -gamma, to the linear solver
  linsol_.setInput(M_.output(),0);

  // Prepare the solution of the linear system (e.g. factorize) -- only if the linear solver inherits from LinearSolver
  linsol_.prepare();

  // Log time duration
  time1 = clock();
  t_lsetup_fac += double(time1-time2)/CLOCKS_PER_SEC;
}

void CVodesInternal::lsetup(CVodeMem cv_mem, int convfail, N_Vector ypred, N_Vector fpred, booleantype *jcurPtr, N_Vector vtemp1, N_Vector vtemp2, N_Vector vtemp3){
  // Current time
  double t = cv_mem->cv_tn;

  // Scaling factor before J
  double gamma = cv_mem->cv_gamma;

  // Call the preconditioner setup function (which sets up the linear solver)
  psetup(t, ypred, fpred, FALSE, jcurPtr, gamma, vtemp1, vtemp2, vtemp3);
}

int CVodesInternal::lsetup_wrapper(CVodeMem cv_mem, int convfail, N_Vector ypred, N_Vector fpred, booleantype *jcurPtr, N_Vector vtemp1, N_Vector vtemp2, N_Vector vtemp3){
  try{
    CVodesInternal *this_ = (CVodesInternal*)(cv_mem->cv_lmem);
    casadi_assert(this_);
    this_->lsetup(cv_mem, convfail, ypred, fpred, jcurPtr, vtemp1, vtemp2, vtemp3);
    return 0;
  } catch(exception& e){
    cerr << "lsetup failed: " << e.what() << endl;;
    return 1;
  }
}

void CVodesInternal::lsolve(CVodeMem cv_mem, N_Vector b, N_Vector weight, N_Vector ycur, N_Vector fcur){
  // Current time
  double t = cv_mem->cv_tn;

  // Scaling factor before J
  double gamma = cv_mem->cv_gamma;

  // Accuracy
  double delta = 0.0;
  
  // Left/right preconditioner
  int lr = 1;
  
  // Call the preconditioner solve function (which solves the linear system)
  psolve(t, ycur, fcur, b, b, gamma, delta, lr, 0);
}

int CVodesInternal::lsolve_wrapper(CVodeMem cv_mem, N_Vector b, N_Vector weight, N_Vector ycur, N_Vector fcur){
  try{
    CVodesInternal *this_ = (CVodesInternal*)(cv_mem->cv_lmem);
    casadi_assert(this_);
    this_->lsolve(cv_mem, b, weight, ycur, fcur);
    return 0;
  } catch(exception& e){
    cerr << "lsolve failed: " << e.what() << endl;;
    return 1;
  }
}

void CVodesInternal::initUserDefinedLinearSolver(){
  // Make sure that a Jacobian has been provided
  if(M_.isNull()) throw CasadiException("CVodesInternal::initUserDefinedLinearSolver(): No Jacobian has been provided.");

  // Make sure that a linear solver has been providided
  if(linsol_.isNull()) throw CasadiException("CVodesInternal::initUserDefinedLinearSolver(): No user defined linear solver has been provided.");

  //  Set fields in the IDA memory
  CVodeMem cv_mem = CVodeMem(mem_);
  cv_mem->cv_lmem   = this;
  cv_mem->cv_lsetup = lsetup_wrapper;
  cv_mem->cv_lsolve = lsolve_wrapper;
  cv_mem->cv_setupNonNull = TRUE;
}

void CVodesInternal::setLinearSolver(const LinearSolver& linsol, const FX& jac){
  linsol_ = linsol;
  M_ = jac;
}


bool CVodesInternal::symbjac(){
  SXFunction f = shared_cast<SXFunction>(f_);
  SXFunction q = shared_cast<SXFunction>(q_);

  return !f.isNull() && q_.isNull() == q.isNull();
}

Integrator CVodesInternal::jac(int iind, int oind){
  casadi_assert_message(oind==INTEGRATOR_XF,"CVodesInternal::jacobian: Not derivative of state");
  
  // Type cast to SXMatrix (currently supported)
  SXFunction f = shared_cast<SXFunction>(f_);
  casadi_assert(!f.isNull());
  
  SXFunction q = shared_cast<SXFunction>(q_);
  casadi_assert(q_.isNull() == q.isNull());
    
  // Generate Jacobians with respect to state, state derivative and parameters
  SXMatrix df_dy = f.jac(ODE_Y,ODE_RHS);
  
  // Number of sensitivities
  int ns;
  switch(iind){
    case INTEGRATOR_P: ns = np_; break;
    case INTEGRATOR_X0: ns = nx_; break;
    default: casadi_assert_message(false,"iind must be INTEGRATOR_P or INTEGRATOR_X0");
  }
  
  // Sensitivities and derivatives of sensitivities
  SXMatrix ysens = symbolic("ysens",ny_,ns);
  
  // Sensitivity ODE
  SXMatrix rhs_s = prod(df_dy,ysens);
  if(iind==INTEGRATOR_P) rhs_s += f.jac(ODE_P,ODE_RHS);

  // Augmented ODE
  SXMatrix faug = vec(horzcat(f.outputSX(oind),rhs_s));

  // Input arguments for the augmented DAE
  vector<SXMatrix> faug_in(ODE_NUM_IN);
  faug_in[ODE_T] = f.inputSX(ODE_T);
  faug_in[ODE_Y] = vec(horzcat(f.inputSX(ODE_Y),ysens));
  faug_in[ODE_P] = f.inputSX(ODE_P);
  
  // Create augmented DAE function
  SXFunction ffcn_aug(faug_in,faug);
  
  // Augmented quadratures
  SXFunction qfcn_aug;
  
  if(!q.isNull()){
    // Now lets do the same for the quadrature states
    SXMatrix dq_dy = q.jac(ODE_Y,ODE_RHS);
    
    // Sensitivity quadratures
    SXMatrix q_s = prod(dq_dy,ysens);
    if(iind==INTEGRATOR_P) q_s += q.jac(ODE_P,ODE_RHS);

    // Augmented quadratures
    SXMatrix qaug = vec(horzcat(q.outputSX(oind),q_s));

    // Input to the augmented quadratures
    vector<SXMatrix> qaug_in(ODE_NUM_IN);
    qaug_in[ODE_T] = q.inputSX(ODE_T);
    qaug_in[ODE_Y] = vec(horzcat(q.inputSX(ODE_Y),ysens));
    qaug_in[ODE_P] = q.inputSX(ODE_P);

    // Create augmented DAE function
    qfcn_aug = SXFunction(qaug_in,qaug);
  }
  
  // Create integrator instance
  CVodesIntegrator integrator(ffcn_aug,qfcn_aug);

  // Set options
  integrator.setOption(dictionary());
  integrator.setOption("nrhs",1+ns);
  
  // Pass linear solver
  if(!linsol_.isNull()){
    LinearSolver linsol_aug = shared_cast<LinearSolver>(linsol_.clone());
    linsol_aug->nrhs_ = 1+ns;
    integrator.setLinearSolver(linsol_aug,jac_f_);
  }
  
  return integrator;
}

FX CVodesInternal::getJacobian(){
  return M_;
}
  
LinearSolver CVodesInternal::getLinearSolver(){
  return linsol_;
}

vector<int> CVodesInternal::jacmap(int ns){
  vector<int> jmap(nx_*(1+ns));
  for(int i=0; i<1+ns; ++i){
    for(int j=0; j<ny_; ++j)
      jmap[j+nx_*i] = j+ny_*i;
    for(int j=0; j<nq_; ++j)
      jmap[ny_+j+nx_*i] = ny_*(1+ns) + j + nq_*i;
  }
  return jmap;
}

} // namespace Sundials
} // namespace CasADi

