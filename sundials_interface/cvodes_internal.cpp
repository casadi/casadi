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

#include "cvodes_internal.hpp"
#include "casadi/stl_vector_tools.hpp"
#include <cassert>

using namespace std;
namespace CasADi{
namespace Sundials{

CVodesInternal* CVodesInternal::clone() const{
  // Copying initialized objects are not allowed since they contain pointers
  if(is_init) throw CasadiException("CVodesInternal::clone: cannot clone an initialized object");
  
  // Return a shallow copy
  return new CVodesInternal(*this);
}
  
int CVodesInternal::getNX(const FX& f, const FX& q){
  // Number of states
  int nx = f.output().numel();
  
  // Add quadratures, if any_
  if(!q.isNull()) nx += q.output().numel();
  
  return nx;
}

int CVodesInternal::getNP(const FX& f){
  return f.input(ODE_P).numel();
}
  
CVodesInternal::CVodesInternal(const FX& f, const FX& q) : IntegratorInternal(getNX(f,q), getNP(f)), f_(f), q_(q){
  addOption("linear_multistep_method",     OT_STRING,  "bdf"); // "bdf" or "adams"
  addOption("nonlinear_solver_iteration",  OT_STRING,  "newton"); // "newton" or "functional"
  addOption("fsens_all_at_once",           OT_BOOLEAN,true); // calculate all right hand sides of the sensitivity equations at once

  mem_ = 0;

  y0_ = y_ = 0;
  yQ0_ = yQ_ = 0;

  is_init = false;

  // Get dimensions
  ny_ = f.output().numel();
  nq_ = q.isNull() ? 0 : q.output().numel();

}

CVodesInternal::~CVodesInternal(){ 
  if(mem_) CVodeFree(&mem_);

  // ODE integration
  if(y0_) N_VDestroy_Serial(y0_);
  if(y_) N_VDestroy_Serial(y_);
  if(yQ0_) N_VDestroy_Serial(yQ0_);
  if(yQ_) N_VDestroy_Serial(yQ_);
  
  // Forward problem
  for(vector<N_Vector>::iterator it=yS0_.begin(); it != yS0_.end(); ++it)   N_VDestroy_Serial(*it);
  for(vector<N_Vector>::iterator it=yS_.begin(); it != yS_.end(); ++it)     N_VDestroy_Serial(*it);
  for(vector<N_Vector>::iterator it=yQS0_.begin(); it != yQS0_.end(); ++it) N_VDestroy_Serial(*it);
  for(vector<N_Vector>::iterator it=yQS_.begin(); it != yQS_.end(); ++it)   N_VDestroy_Serial(*it);
  
  // Adjoint problem
  for(vector<N_Vector>::iterator it=yB0_.begin(); it != yB0_.end(); ++it)   N_VDestroy_Serial(*it);
  for(vector<N_Vector>::iterator it=yB_.begin(); it != yB_.end(); ++it)     N_VDestroy_Serial(*it);
//  for(vector<N_Vector>::iterator it=yQB0_.begin(); it != yQB0_.end(); ++it) N_VDestroy_Serial(*it);
  for(vector<N_Vector>::iterator it=yQB_.begin(); it != yQB_.end(); ++it)   N_VDestroy_Serial(*it);
}

void CVodesInternal::init(){
  IntegratorInternal::init();

  // Init ODE rhs function and quadrature functions
  f_.init();
  if(!q_.isNull()) q_.init();

  // Quick return if already initialized
  if(is_init){
    reset(ad_order_,ad_order_);
    return;
  }

  // Sundials return flag
  int flag;

  int lmm; // linear multistep method
  if(getOption("linear_multistep_method")=="adams")  lmm = CV_ADAMS;
  else if(getOption("linear_multistep_method")=="bdf") lmm = CV_BDF;
  else throw CasadiException("Unknown linear multistep method");

  int iter; // nonlinear solver iteration
  if(getOption("nonlinear_solver_iteration")=="newton") iter = CV_NEWTON;
  else if(getOption("nonlinear_solver_iteration")=="functional") iter = CV_FUNCTIONAL;
  else throw CasadiException("Unknown nonlinear solver iteration");

  // Create CVodes memory block
  mem_ = CVodeCreate(lmm,iter);
  if(mem_==0) throw CasadiException("CVodeCreate: Creation failed");

  // Allocate n-vectors for ivp
  y0_ = N_VMake_Serial(ny_,&input(INTEGRATOR_X0).data()[0]);
  y_ = N_VMake_Serial(ny_,&output(INTEGRATOR_XF).data()[0]);

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
      // Create jacobian
      jac_f_ = f_.jacobian(ODE_Y,ODE_RHS);
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
  } else throw CasadiException("Unknown linear solver ");

  // Set user data
  flag = CVodeSetUserData(mem_,this);
  if(flag!=CV_SUCCESS) cvodes_error("CVodeSetUserData",flag);

  // Quadrature equations
  if(nq_>0){
    // Allocate n-vectors for quadratures
    yQ0_ = N_VMake_Serial(nq_,&input(INTEGRATOR_X0).data()[ny_]);
    yQ_ = N_VMake_Serial(nq_,&output(INTEGRATOR_XF).data()[ny_]);

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
  
// Sensitivities
  if(ad_order_>0){

    // Forward sensitivity problem
    if(nfdir_>0){
      // Allocate n-vectors
      yS0_.resize(nfdir_);
      yS_.resize(nfdir_);
      for(int i=0; i<nfdir_; ++i){
        yS0_[i] = N_VMake_Serial(ny_,&input(INTEGRATOR_X0).dataF(i)[0]);
        yS_[i] = N_VMake_Serial(ny_,&output(INTEGRATOR_XF).dataF(i)[0]);
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
      flag = CVodeSetSensParams(mem_,&input(INTEGRATOR_P).data()[0],0,0);
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
  yB0_.resize(nadir_);
  yB_.resize(nadir_);
  for(int i=0; i<nadir_; ++i){
    yB0_[i] = N_VMake_Serial(ny_,&output(INTEGRATOR_XF).dataA(i)[0]);
    yB_[i] = N_VMake_Serial(ny_,&input(INTEGRATOR_X0).dataA(i)[0]);
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
      interpType = CV_HERMITE;
    else if(getOption("interpolation_type")=="polynomial")
      interpType = CV_POLYNOMIAL;
    else throw CasadiException("\"interpolation_type\" must be \"hermite\" or \"polynomial\"");
      
      // Initialize adjoint sensitivities
      flag = CVodeAdjInit(mem_, Nd, interpType);
      if(flag != CV_SUCCESS) cvodes_error("CVodeAdjInit",flag);
  }
  
  for(int dir=0; dir<nadir_; ++dir){

      // Create backward problem (use the same lmm and iter)
      flag = CVodeCreateB(mem_, lmm, iter, &whichB_[dir]);
      if(flag != CV_SUCCESS) cvodes_error("CVodeCreateB",flag);
      
      // Initialize the backward problem
      double tB0 = 1.0; // will be overwritten during re-init
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

  
    } // enable asens
    
  } // ad_order>0
  
  is_init = true;
}

void CVodesInternal::rhs(double t, const double* y, double* ydot){
    // Pass input
  f_.input(ODE_T).set(t);
  f_.input(ODE_Y).set(y);
  f_.input(ODE_P).set(input(INTEGRATOR_P).data());

    // Evaluate
  f_.evaluate();
    
    // Get results
  f_.output().get(ydot);
}

int CVodesInternal::rhs_wrapper(double t, N_Vector y, N_Vector ydot, void *user_data){
try{
    assert(user_data);
    CVodesInternal *this_ = (CVodesInternal*)user_data;
    this_->rhs(t,NV_DATA_S(y),NV_DATA_S(ydot));
    return 0;
  } catch(exception& e){
    cerr << "rhs failed: " << e.what() << endl;
    return 1;
  }
}
  
void CVodesInternal::reset(int fsens_order, int asens_order){
  fsens_order_ = fsens_order;
  asens_order_ = asens_order;
  
/*  if(asens_order_>0) assert(nadir_>0);
  if(fsens_order_>0) assert(nfdir_>0);*/
  
  // Get the time horizon
  double t0 = input(INTEGRATOR_T0).data()[0];
  double tf = input(INTEGRATOR_TF).data()[0];
  t_ = t0;

  // quick-fix
  t0 -= input(INTEGRATOR_T0).data()[0];
  tf -= input(INTEGRATOR_T0).data()[0];
  
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
  // quick-fix
  t_out -= input(INTEGRATOR_T0).data()[0];
  
  int flag;
  
  // tolerance
  double ttol = 1e-9;
  if(fabs(t_-t_out)<ttol){
    copy(input(INTEGRATOR_X0).data().begin(),input(INTEGRATOR_X0).data().end(),output(INTEGRATOR_XF).data().begin());
    if(fsens_order_>0){
      for(int i=0; i<nfdir_; ++i){
        copy(input(INTEGRATOR_X0).dataF(i).begin(),input(INTEGRATOR_X0).dataF(i).end(),output(INTEGRATOR_XF).dataF(i).begin());
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
  double tf = input(INTEGRATOR_TF).data()[0];

  // quick-fix
  tf -= input(INTEGRATOR_T0).data()[0];
  
  int flag;
  
  for(int dir=0; dir<nadir_; ++dir){
    flag = CVodeReInitB(mem_, whichB_[dir], tf, yB0_[dir]);
    if(flag != CV_SUCCESS) cvodes_error("CVodeReInitB",flag);

    N_VConst(0.0,yQB_[dir]);
    flag = CVodeQuadReInitB(mem_,whichB_[dir],yQB_[dir]);
    if(flag!=CV_SUCCESS) cvodes_error("CVodeQuadReInitB",flag);
  }
}

void CVodesInternal::integrateAdj(double t_out){
  int flag;

  // quick-fix
  t_out -= input(INTEGRATOR_T0).data()[0];
  
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
    assert(eh_data);
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
  assert(Ns==nfdir_);
  
    // Pass input
  f_.input(ODE_T).set(t);
  f_.input(ODE_Y).set(NV_DATA_S(y));
  f_.input(ODE_P).set(input(INTEGRATOR_P).data());

  for(int i=0; i<nfdir_; ++i){

    // Pass forward seeds
    f_.input(ODE_T).setF(0.0);
    f_.input(ODE_Y).setF(NV_DATA_S(yS[i]));
    f_.input(ODE_P).setF(input(INTEGRATOR_P).dataF(i));
    
    // Evaluate the AD forward algorithm
    f_.evaluate(1,0);
  
    // Get the output seeds
    f_.output().getF(NV_DATA_S(ySdot[i]));
  }
}

int CVodesInternal::rhsS_wrapper(int Ns, double t, N_Vector y, N_Vector ydot, N_Vector *yS, N_Vector *ySdot, void *user_data, N_Vector tmp1, N_Vector tmp2){
try{
    assert(user_data);
    CVodesInternal *this_ = (CVodesInternal*)user_data;
    this_->rhsS(Ns,t,y,ydot,yS,ySdot,tmp1,tmp2);
    return 0;
  } catch(exception& e){
    cerr << "fs failed: " << e.what() << endl;
    return 1;
  }
}

void CVodesInternal::rhsS1(int Ns, double t, N_Vector y, N_Vector ydot, int iS, N_Vector yS, N_Vector ySdot, N_Vector tmp1, N_Vector tmp2){
  assert(Ns==nfdir_);
  
    // Pass input
  f_.input(ODE_T).set(t);
  f_.input(ODE_Y).set(NV_DATA_S(y));
  f_.input(ODE_P).set(input(INTEGRATOR_P).data());

  // Pass forward seeds
  f_.input(ODE_T).setF(0.0);
  f_.input(ODE_Y).setF(NV_DATA_S(yS));
  f_.input(ODE_P).setF(input(INTEGRATOR_P).dataF(iS));
    
  // Evaluate the AD forward algorithm
  f_.evaluate(1,0);
  
  // Get the output seeds
  f_.output().getF(NV_DATA_S(ySdot));
}

int CVodesInternal::rhsS1_wrapper(int Ns, double t, N_Vector y, N_Vector ydot, int iS, N_Vector yS, N_Vector ySdot, void *user_data, N_Vector tmp1, N_Vector tmp2){
try{
    assert(user_data);
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
    assert(user_data);
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
  q_.input(ODE_T).set(t);
  q_.input(ODE_Y).set(yy);
  q_.input(ODE_P).set(input(INTEGRATOR_P).data());

  // Evaluate
  q_.evaluate();
    
  // Get results
  q_.output().get(rhsQ);
}

void CVodesInternal::rhsQS(int Ns, double t, N_Vector y, N_Vector *yS, N_Vector yQdot, N_Vector *rhsvalQS, N_Vector tmp1, N_Vector tmp2){
  assert(Ns==nfdir_);
  
  // Pass input
  q_.input(ODE_T).set(t);
  q_.input(ODE_Y).set(NV_DATA_S(y));
  q_.input(ODE_P).set(input(INTEGRATOR_P).data());

  // Pass forward seeds
  for(int i=0; i<nfdir_; ++i){
    q_.input(ODE_T).setF(0.0);
    q_.input(ODE_Y).setF(NV_DATA_S(yS[i]));
    q_.input(ODE_P).setF(input(INTEGRATOR_P).dataF(i));

    // Evaluate the AD forward algorithm
    q_.evaluate(1,0);
      
    // Get the output seeds
    q_.output().getF(NV_DATA_S(rhsvalQS[i]));
  }
}

int CVodesInternal::rhsQS_wrapper(int Ns, double t, N_Vector y, N_Vector *yS, N_Vector yQdot, N_Vector *rhsvalQS, void *user_data, N_Vector tmp1, N_Vector tmp2){
try{
//    assert(user_data);
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
  f_.input(ODE_T).set(t);
  f_.input(ODE_Y).set(y);
  f_.input(ODE_P).set(input(INTEGRATOR_P).data());

  // Pass adjoint seeds
  f_.output(ODE_RHS).setA(yB);
  
  // Evaluate and tape
  f_.evaluate(0,1);

  // Save to output
  const vector<double>& fres = f_.input(ODE_Y).dataA();
  for(int i=0; i<ny_; ++i)
    yBdot[i] = -fres[i];

  // If quadratures are included
  if(nq_>0){
    // Pass input to quadratures
    q_.input(ODE_T).set(t);
    q_.input(ODE_Y).set(y);
    q_.input(ODE_P).set(input(INTEGRATOR_P).data());

    // Pass adjoint seeds
    q_.output(ODE_RHS).setA(&output(INTEGRATOR_XF).dataA()[ny_]);

    // Evaluate
    q_.evaluate(0,1);
    
    // Get the input seeds
    const vector<double>& qres = q_.input(ODE_Y).dataA();
    
    // Copy to result
    for(int i=0; i<ny_; ++i)
      yBdot[i] -= qres[i];
  }
}

int CVodesInternal::rhsB_wrapper(double t, N_Vector y, N_Vector yB, N_Vector yBdot, void *user_dataB){
try{
    assert(user_dataB);
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
    assert(user_dataB);
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
  f_.input(ODE_T).set(t);
  f_.input(ODE_Y).set(y);
  f_.input(ODE_P).set(input(INTEGRATOR_P).data());

  // Pass adjoint seeds
  f_.output(ODE_RHS).setA(yB);

  // Evaluate
  f_.evaluate(0,1);

  // Save to output
  f_.input(ODE_P).getA(qBdot);
  
  // If quadratures are included
  if(nq_>0){
    // Pass input to quadratures
    q_.input(ODE_T).set(t);
    q_.input(ODE_Y).set(y);
    q_.input(ODE_P).set(input(INTEGRATOR_P).data());

    // Pass adjoint seeds
    q_.output(ODE_RHS).setA(&output(INTEGRATOR_XF).dataA()[ny_]);

    // Evaluate
    q_.evaluate(0,1);
    
    // Get the input seeds
    const vector<double>& qres = q_.input(ODE_P).dataA();
    
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
    assert(user_data);
    CVodesInternal *this_ = (CVodesInternal*)user_data;
    assert(this_->f_.output(ODE_RHS).dataF().size() == this_->ny_);
    assert(NV_LENGTH_S(v) == this_->ny_);
    assert(NV_LENGTH_S(Jv) == this_->ny_);
    this_->jtimes(NV_DATA_S(v),NV_DATA_S(Jv),t,NV_DATA_S(y),NV_DATA_S(fy),NV_DATA_S(tmp));
    return 0;
  } catch(exception& e){
    cerr << "jtimes failed: " << e.what() << endl;;
    return 1;
  }
}

void CVodesInternal::jtimes(const double *v, double* Jv, double t, const double* y, const double* fy, double* tmp){
  // Pass input
  f_.input(ODE_T).set(t);
  f_.input(ODE_Y).set(y);
  f_.input(ODE_P).set(input(INTEGRATOR_P).data());

  // Pass input seeds
  f_.input(ODE_T).setF(0.0);
  f_.input(ODE_Y).setF(v);
  fill_n(f_.input(ODE_P).dataF().begin(),np_,0.0);
  
  // Evaluate
  f_.evaluate(1,0);

  // Get the output seeds
  f_.output(ODE_RHS).getF(Jv);
}

int CVodesInternal::djac_wrapper(int N, double t, N_Vector y, N_Vector fy, DlsMat Jac, void *user_data,N_Vector tmp1, N_Vector tmp2, N_Vector tmp3){
  try{
    assert(user_data);
    CVodesInternal *this_ = (CVodesInternal*)user_data;
    this_->djac(N, t, y, fy, Jac, tmp1, tmp2, tmp3);
    return 0;
  } catch(exception& e){
    cerr << "djac failed: " << e.what() << endl;;
    return 1;
  }
}

void CVodesInternal::djac(int N, double t, N_Vector y, N_Vector fy, DlsMat Jac, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3){
    // Retrieve the structure that holds the Jacobian
    _DlsMat *Jac_ = (_DlsMat *)Jac;

    // Get a pointer to the first element of the of the jacobian
    double *jdata = Jac->data;

    // Retrieve the leading dimension of the jacobian
    int ld_jac = Jac->ldim;

    // Pass inputs to the jacobian function
    jac_f_.input(ODE_T).set(t);
    jac_f_.input(ODE_Y).set(NV_DATA_S(y));
    jac_f_.input(ODE_P).set(f_.input(ODE_P).data());

    // Evaluate
    jac_f_.evaluate();
  
    // Save the results
    assert(ld_jac == ny_);
    jac_f_.output().get(jdata); // transpose?
}

int CVodesInternal::bjac_wrapper(int N, int mupper, int mlower, double t, N_Vector y, N_Vector fy, DlsMat Jac, void *user_data,     
                        N_Vector tmp1, N_Vector tmp2, N_Vector tmp3){
  try{
    assert(user_data);
    CVodesInternal *this_ = (CVodesInternal*)user_data;
    this_->bjac(N, mupper, mlower, t, y, fy, Jac, tmp1, tmp2, tmp3);
    return 0;
  } catch(exception& e){
    cerr << "djac failed: " << e.what() << endl;;
    return 1;
  }
}

void CVodesInternal::bjac(int N, int mupper, int mlower, double t, N_Vector y, N_Vector fy, DlsMat Jac,     
                        N_Vector tmp1, N_Vector tmp2, N_Vector tmp3){
  
    throw CasadiException("bjac not implemented");
}

void CVodesInternal::setStopTime(double tf){
  // Set the stop time of the integration -- don't integrate past this point
  int flag = CVodeSetStopTime(mem_, tf);
  if(flag != CV_SUCCESS) cvodes_error("CVodeSetStopTime",flag);
}

} // namespace Sundials
} // namespace CasADi

