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
#include "casadi/stl_vector_tools.hpp"
#include "casadi/fx/linear_solver_internal.hpp"
#include "casadi/fx/sx_function_internal.hpp"
#include "casadi/fx/mx_function.hpp"
#include "casadi/sx/sx_tools.hpp"
#include <acado_integrators.hpp>

using namespace std;
namespace CasADi{

  
AcadoIntegratorInternal::AcadoIntegratorInternal(const FX& f, const FX& q) : IntegratorInternal(f,q){
  addOption("time_dependence",   OT_BOOLEAN,  true, "Explicit depencency of time in the DAE");
  addOption("num_algebraic",     OT_INTEGER,  0,    "Number of algebraic states");

  casadi_warning("AcadoIntegrator interface is under development");
  t_ = 0;
  xd_= 0;
  xa_ = 0;
  p_= 0;
  arg_ = 0;
  diff_eq_ = 0;
  integrator_ = 0;
}

AcadoIntegratorInternal::~AcadoIntegratorInternal(){
  if(t_!= 0) delete t_;
  if(xd_!= 0) delete[] xd_;
  if(xa_!= 0) delete[] xa_;
  if(p_!= 0) delete[] p_;
  if(arg_!= 0) delete arg_;
  if(diff_eq_!= 0) delete diff_eq_;
  if(integrator_ != 0) delete integrator_;
}

AcadoIntegratorInternal* AcadoIntegratorInternal::clone() const{
  // Return a deep copy
  AcadoIntegratorInternal* node = new AcadoIntegratorInternal(f_,q_);
  node->setOption(dictionary());
  return node;
}

void AcadoIntegratorInternal::init(){
  // Init ODE rhs function and quadrature functions, jacobian function
  if(!f_.isInit()) f_.init();
  casadi_assert(q_.isNull());
  
  int nx = f_.input(DAE_Y).size();
  int np = f_.input(DAE_P).size();
  setDimensions(nx,np);
  
  // Call the base class init
  IntegratorInternal::init();

  // Get the number of differential and algebraic equations
  nxa_ = getOption("num_algebraic");
  nxd_ = nx_ - nxa_;
  nt_ = bool(getOption("time_dependence")) ? 1 : 0;
  
  // Create wrapper function
  rhs_ = AcadoFunction(f_);
  rhs_.init();

  cout << "nxd_ = " << nxd_ << endl;
  cout << "nxa_ = " << nxa_ << endl;
  cout << "np_ = " << np_ << endl;
  
  
  // Declare ACADO variables
  t_ = new ACADO::TIME();
  xd_ = new ACADO::DifferentialState[nxd_];
  xa_ = new ACADO::AlgebraicState[nxa_];
  p_ = new ACADO::Parameter[np_];

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
}

void AcadoIntegratorInternal::integrate(double t_out){
  
    // DEFINE INITIAL VALUES:
    // ----------------------
    double *x0 = getPtr(input(INTEGRATOR_X0));
    double *z0 = x0 + nxd_;
    double *pp = getPtr(input(INTEGRATOR_P));

    ACADO::Grid interval( 0.0, 1.0, 100 );


    // START THE INTEGRATION:
    // ----------------------
    integrator_->integrate( interval, x0, z0, pp );

    ACADO::VariablesGrid differentialStates;
    ACADO::VariablesGrid algebraicStates   ;
    ACADO::VariablesGrid intermediateStates;

    integrator_->getX ( differentialStates );
    integrator_->getXA( algebraicStates    );
    integrator_->getI ( intermediateStates );

    
    cout << "differentialStates:" << endl;
    differentialStates.print();
    cout << endl;

    cout << "algebraicStates:" << endl;
    algebraicStates.print();
    cout << endl;

    cout << "intermediateStates:" << endl;
    intermediateStates.print();
    cout << endl;
}



} // namespace CasADi

