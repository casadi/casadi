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

#include "acado_integrator_internal.hpp"
#include "symbolic/stl_vector_tools.hpp"
#include "symbolic/fx/linear_solver_internal.hpp"
#include "symbolic/fx/mx_function.hpp"
#include "symbolic/sx/sx_tools.hpp"
#include <acado_integrators.hpp>

using namespace std;
namespace CasADi{

  
AcadoIntegratorInternal::AcadoIntegratorInternal(const FX& f, const FX& g) : IntegratorInternal(f,g){
  addOption("time_dependence",   OT_BOOLEAN,  true, "Explicit depencency of time in the DAE");
  addOption("num_algebraic",     OT_INTEGER,  0,    "Number of algebraic states");
  addOption("num_grid_points",   OT_INTEGER,  2,  "Number of uniformly distributed grid points for obtaining the solution, does not influence the integration steps");
  setNull();
}

AcadoIntegratorInternal::~AcadoIntegratorInternal(){
  freeMem();
}

void AcadoIntegratorInternal::setNull(){
  t_ = 0;
  xd_= 0;
  xa_ = 0;
  p_= 0;
  arg_ = 0;
  diff_eq_ = 0;
  integrator_ = 0;
  interval_ = 0;
  algebraicStates_ = 0;
  differentialStates_ = 0;
  tmp_ = 0;
  x_tmp_ = 0;
  p_tmp_ = 0;
  frozen_grid_ = false;
}

void AcadoIntegratorInternal::freeMem(){
  if(t_!= 0) delete t_;
  if(xd_!= 0) delete[] xd_;
  if(xa_!= 0) delete[] xa_;
  if(p_!= 0) delete[] p_;
  if(arg_!= 0) delete arg_;
  if(diff_eq_!= 0) delete diff_eq_;
  if(integrator_ != 0) delete integrator_;
  if(interval_ !=0) delete interval_;
  if(differentialStates_!=0) delete differentialStates_;
  if(algebraicStates_ !=0) delete algebraicStates_;
  if(tmp_!=0) delete tmp_;
  if(x_tmp_!=0) delete x_tmp_;
  if(p_tmp_!=0) delete p_tmp_;
}

AcadoIntegratorInternal* AcadoIntegratorInternal::clone() const{
  // Return a deep copy
  AcadoIntegratorInternal* node = new AcadoIntegratorInternal(f_,g_);
  node->setOption(dictionary());
  return node;
}

void AcadoIntegratorInternal::init(){
  // Call the base class init
  IntegratorInternal::init();
  
  // Free memory and set pointers to NULL
  freeMem();
  setNull();
  
  // The following will make sure that no data lingers in ACADO
  ACADO::AlgebraicState().clearStaticCounters();
  ACADO::Control().clearStaticCounters();
  ACADO::DifferentialState().clearStaticCounters();
  ACADO::DifferentialStateDerivative().clearStaticCounters();
  ACADO::Disturbance().clearStaticCounters();
  ACADO::IntegerControl().clearStaticCounters();
  ACADO::IntegerParameter().clearStaticCounters();
  ACADO::IntermediateState().clearStaticCounters();
  ACADO::Parameter().clearStaticCounters();
  
  // Get the number of differential and algebraic equations
  nxa_ = getOption("num_algebraic");
  nxd_ = nx_ - nxa_;
  nt_ = bool(getOption("time_dependence")) ? 1 : 0;
  
  // Create wrapper function
  rhs_ = AcadoFunction(f_);
  rhs_.init();

  // Declare ACADO variables
  t_ = new ACADO::TIME();
  xd_ = new ACADO::DifferentialState[nxd_];
  xa_ = new ACADO::AlgebraicState[nxa_];
  p_ = new ACADO::Parameter[np_];

  // Temporary vector
  x_tmp_ = new ACADO::Vector(nxd_);
  if(np_>0) p_tmp_ = new ACADO::Vector(np_);
  
  // Augmented state vector
  arg_ = new ACADO::IntermediateState(nt_+nxd_+nxa_+np_);
  int ind=0;
  for(int i=0; i<nt_; ++i)    (*arg_)(ind++) = t_[i];
  for(int i=0; i<nxd_; ++i)   (*arg_)(ind++) = xd_[i];
  for(int i=0; i<nxa_; ++i)   (*arg_)(ind++) = xa_[i];
  for(int i=0; i<np_; ++i)    (*arg_)(ind++) = p_[i];

  // Differential equation
  diff_eq_ = new ACADO::DifferentialEquation();
  *diff_eq_ << (*rhs_.fcn_)(*arg_);

  // Allocate an integrator
  integrator_ = new ACADO::IntegratorBDF(*diff_eq_);

  // Grid points
  num_grid_points_ = getOption("num_grid_points");
  interval_ = new ACADO::Grid( t0_, tf_, num_grid_points_);
  
  // Variablesgrid for the solution
  if(nxd_>0) differentialStates_ = new ACADO::VariablesGrid();
  if(nxa_>0) algebraicStates_ = new ACADO::VariablesGrid();
  
  // Temporary
  tmp_ = new ACADO::Vector();
}

void AcadoIntegratorInternal::reset(int nsens, int nsensB, int nsensB_store){
  int nfdir = 0; // NOTE: need to update the function below to the new integrator formulation
  if(nfdir>0){
    integrator_->freezeAll();
  }
  has_been_integrated_ = false;
  nfsens_ = nfdir;
  nasens_ = 0;
}

void AcadoIntegratorInternal::integrate(double t_out){
  // Integrate the system if this has not already been done
  if(!has_been_integrated_){
    // Initial conditions
    double *x0 = getPtr(input(INTEGRATOR_X0));
    double *z0 = x0 + nxd_;
    
    // Parameters
    double *pp = getPtr(input(INTEGRATOR_P));

    // Integrate
    integrator_->integrate( *interval_, x0, z0, pp );

    // Get solution
    if(nxd_>0) integrator_->getX (*differentialStates_);
    if(nxa_>0) integrator_->getXA(*algebraicStates_);
    
    // Forward sensitivities
    for(int dir=0; dir<nfsens_; ++dir){
      // Order
      const int order = 1;
      
      // Get pointer
      double *xseed = getPtr(fwdSeed(INTEGRATOR_X0,dir));
      double *pseed = getPtr(fwdSeed(INTEGRATOR_P,dir));
      
      // Pass seeds
      *x_tmp_ << xseed;
      if(np_>0) *p_tmp_ << pseed;
      integrator_->setForwardSeed(order, *x_tmp_, *p_tmp_);
      
      // Calculate sensitivities
      integrator_->integrateSensitivities();
      
      // Get pointer
      x_tmp_->setZero();
      integrator_->getForwardSensitivities(*x_tmp_, order);
      double *xsens = getPtr(fwdSens(INTEGRATOR_XF,dir));
      xsens << *x_tmp_;
    }
    
    // Adjoint sensitivities
    for(int dir=0; dir<nasens_; ++dir){
      // Order
      const int order = 1;
      
      // Get pointer
      double *xseed = getPtr(adjSeed(INTEGRATOR_XF,dir));
      
      // Pass seeds
      *x_tmp_ << xseed;
      integrator_->setBackwardSeed(order, *x_tmp_);

      // Clear forward seeds
      x_tmp_->setZero();
      p_tmp_->setZero();
      integrator_->setForwardSeed(order, *x_tmp_, *p_tmp_);
      
      // Calculate sensitivities
      integrator_->integrateSensitivities();
      
      // Get pointer
      double *xsens = getPtr(adjSens(INTEGRATOR_X0,dir));
      double *psens = getPtr(adjSens(INTEGRATOR_P,dir));
      x_tmp_->setZero();
      p_tmp_->setZero();
/*      integrator_->getBackwardSensitivities(ACADO::emptyVector, ACADO::emptyVector, ACADO::emptyVector, ACADO::emptyVector, order);*/
/*      xsens << *x_tmp_;
      if(np_>0) psens << *p_tmp_;*/
    }
    
    // Unfreeze the grid
    if(!frozen_grid_){
      integrator_->unfreeze();
    }
    
    // Mark as integrated
    has_been_integrated_ = true;
  }
  
  // Get the grid point corresponding to the output time
  double grid_point_cont = (num_grid_points_-1)*(t_out-t0_)/(tf_-t0_);
  
  // Get the time point before and after
  int grid_point_before = std::floor(grid_point_cont);
  int grid_point_after  = std::ceil(grid_point_cont);
  
  // Get data
  for(int k=0; k< (nxa_==0 ? 1 : 2); ++k){
    // Pointer to the data
    double* xf = getPtr(output(INTEGRATOR_XF));
    if(k==1) xf += nxd_;
    
    // Variablesgrid
    ACADO::VariablesGrid* vgrid = k==0 ? differentialStates_ : algebraicStates_;
    
    // Get the value before the given time
    *tmp_ = vgrid->getVector(grid_point_before);
    xf << *tmp_;
    
    // Get the value after the given time and interpolate
    if(grid_point_before!=grid_point_after){
      *tmp_ = vgrid->getVector(grid_point_after);
      
      // Weights
      double w_before = grid_point_after-grid_point_cont;
      double w_after = grid_point_cont-grid_point_before;
      
      // Interpolate
      for(int i=0; i<tmp_->getDim(); ++i){
        xf[i] = w_before*xf[i] + w_after*(*tmp_)(i);
      }
    }
  }
}




} // namespace CasADi

