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
#include <cassert>

using namespace std;
namespace CasADi{

ACADOIntegratorInternal* ACADOIntegratorInternal::clone() const{
  throw CasadiException("ACADOIntegratorInternal::clone: cannot clone");
}

void ACADOIntegratorInternal::ffcn( double *x, double *f){
  // Pass variables
  int ind=0;
  f_.input(ODE_T).set(&x[ind]);    ind += 1;
  f_.input(ODE_Y).set(&x[ind]);    ind += nx_;
  f_.input(ODE_P).set(&x[ind]);    ind += np_;

  // Evaluate function
  f_.evaluate();
  
  // Save result
  f_.output().get(f);
}

void ACADOIntegratorInternal::ffcn_wrapper( double *x, double *f, void *user_data ){
  ACADOIntegratorInternal* this_ = (ACADOIntegratorInternal*)user_data;
  this_->ffcn(x,f);
}

void ACADOIntegratorInternal::ffcn_fwd(int number, double *x, double *seed, double *f, double *df){
  // Pass variables
  int ind=0;
  f_.input(ODE_T).set(&x[ind]);    ind += 1;
  f_.input(ODE_Y).set(&x[ind]);    ind += nx_;
  f_.input(ODE_P).set(&x[ind]);    ind += np_;

  // Pass forward seeds
  ind=0;
  f_.input(ODE_T).setF(&seed[ind]);    ind += 1;
  f_.input(ODE_Y).setF(&seed[ind]);    ind += nx_;
  f_.input(ODE_P).setF(&seed[ind]);    ind += np_;
  
  // Evaluate function
  f_.evaluate(1,0);
  
  // Save the results
  f_.output(ODE_RHS).get(f);
  f_.output(ODE_RHS).getF(df);
}

void ACADOIntegratorInternal::ffcn_fwd_wrapper(int number, double *x, double *seed, double *f, double *df, void *user_data){
  ACADOIntegratorInternal* this_ = (ACADOIntegratorInternal*)user_data;
  this_->ffcn_fwd(number,x,seed,f,df);
}

void ACADOIntegratorInternal::ffcn_adj(int number, double *x, double *seed, double *f, double *df){
  // Pass variables
  int ind=0;
  f_.input(ODE_T).set(&x[ind]);    ind += 1;
  f_.input(ODE_Y).set(&x[ind]);    ind += nx_;
  f_.input(ODE_P).set(&x[ind]);    ind += np_;

  // Pass adjoint seed
  f_.output(ODE_RHS).setF(seed); // TODO: correct this

  // Evaluate function
  f_.evaluate(0,1);
  
  // Save the results
  f_.output(ODE_RHS).get(f);

  ind=0;
  f_.input(ODE_T).getA(&df[ind]);    ind += 1;
  f_.input(ODE_Y).getA(&df[ind]);    ind += nx_;
  f_.input(ODE_P).getA(&df[ind]);    ind += np_;
}

void ACADOIntegratorInternal::ffcn_adj_wrapper(int number, double *x, double *seed, double *f, double *df, void *user_data){
  ACADOIntegratorInternal* this_ = (ACADOIntegratorInternal*)user_data;
  this_->ffcn_adj(number,x,seed,f,df);
}

int ACADOIntegratorInternal::getNX(const FX& f){
  int nx = f.output().numel();
  return nx;
}

int ACADOIntegratorInternal::getNP(const FX& f){
  return f.input(ODE_P).numel();
}
  
ACADOIntegratorInternal::ACADOIntegratorInternal(const FX& f) : IntegratorInternal(getNX(f), getNP(f)), f_(f){
  addOption("print_level",   OT_STRING,  "none"); // "none", "low", "medium", "high", "debug"

  x_ = 0;
  p_ = 0;
  t_ = 0;
  augx_ = 0;
  model_ = 0;
  dae_ = 0;
  integrator_ = 0;
  x_temp_ = 0;
  p_temp_ = 0;
  is_init_ = false;
}
  
ACADOIntegratorInternal::~ACADOIntegratorInternal(){
  if(x_) delete[] x_;
  if(p_) delete[] p_;
  if(t_) delete t_;
  if(augx_) delete augx_;
  if(model_) delete model_;
  if(dae_) delete dae_;
  if(integrator_) delete integrator_;
  if(x_temp_) delete x_temp_;
  if(p_temp_) delete p_temp_;
}


void ACADOIntegratorInternal::init(){
  // call the base class method
  IntegratorInternal::init();
  
  assert(!is_init_);
  
  // initialize the dae function
  f_.init();

  x_ = new ACADO::DifferentialState[nx_];
  p_ = new ACADO::Parameter[np_];
  t_ = new ACADO::TIME();

  model_ = new ACADO::CFunction( nx_, ffcn_wrapper, ffcn_fwd_wrapper, ffcn_adj_wrapper);
  model_->setUserData(this);

  // Define a Right-Hand-Side:
  // -------------------------
  augx_ = new ACADO::IntermediateState(1 + nx_ + np_);
  int iaug = 0;
  (*augx_)(iaug++) = *t_;
  for(int i=0; i<nx_; ++i) (*augx_)(iaug++) = x_[i];
  for(int i=0; i<np_; ++i) (*augx_)(iaug++) = p_[i];
  dae_ = new ACADO::DifferentialEquation();
  *dae_ << (*model_)(*augx_);

  // DEFINE AN INTEGRATOR:
  // ---------------------
//  integrator_ = new ACADO::IntegratorRK45( *dae_ );
  integrator_ = new ACADO::IntegratorBDF( *dae_ );
  
  // set print level
  ACADO::PrintLevel printlevel;
  if(getOption("print_level")=="none")        printlevel = ACADO::NONE;
  else if(getOption("print_level")=="low")    printlevel = ACADO::LOW;
  else if(getOption("print_level")=="medium") printlevel = ACADO::MEDIUM;
  else if(getOption("print_level")=="high")   printlevel = ACADO::HIGH;
  else if(getOption("print_level")=="debug")  printlevel = ACADO::DEBUG;
  else throw CasadiException("Illegal print level. Allowed are \"none\", \"low\", \"medium\", \"high\", \"debug\"");
  integrator_->set(ACADO::INTEGRATOR_PRINTLEVEL, printlevel );
  integrator_->set(ACADO::INTEGRATOR_TOLERANCE, reltol_);
  integrator_->set(ACADO::ABSOLUTE_TOLERANCE, abstol_);
  integrator_->set(ACADO::MAX_NUM_INTEGRATOR_STEPS, max_num_steps_);
//  integrator_->set(ACADO::PRINT_INTEGRATOR_PROFILE, ACADO::YES);
//  integrator_->set(ACADO::MAX_INTEGRATOR_STEPSIZE, 1e-3);
  integrator_->set(ACADO::LINEAR_ALGEBRA_SOLVER , ACADO::SPARSE_LU);
  
  // Temporary acado-vectors
  x_temp_ = new ACADO::Vector(nx_);
  p_temp_ = new ACADO::Vector(np_);
  
  is_init_ = true;
}

void ACADOIntegratorInternal::reset(int fsens_order, int asens_order){
  assert(0);
}

void ACADOIntegratorInternal::resetAdj(){
  assert(0);
}

void ACADOIntegratorInternal::integrate(double t_out){
  assert(0);
}

void ACADOIntegratorInternal::evaluate(int fsens_order, int asens_order){
  double t_start    =  input(INTEGRATOR_T0).data()[0];
  double t_end      =  input(INTEGRATOR_TF).data()[0];
  double *x0  = &input(INTEGRATOR_X0).data()[0];
  double *p   = &input(INTEGRATOR_P).data()[0];

  if(fsens_order>0 || asens_order>0)
    integrator_->freezeAll(); // why here and not later?
  integrator_->integrate( t_start, t_end, x0, 0, p, 0 );

  // Save results
  integrator_->getX( *x_temp_);
  vector<double> &xf = output(INTEGRATOR_XF).data();
  for(int i=0; i<nx_; ++i) xf[i] = (*x_temp_)(i);
  
  // First order sensitivities
  if(fsens_order>=1){
    // Pass forward seeds
    const vector<double> &x_seed1 = input(INTEGRATOR_X0).dataF();
    for(int i=0; i<nx_; ++i) (*x_temp_)(i) = x_seed1[i];
    
    const vector<double> &p_seed1 = input(INTEGRATOR_P).dataF();
    for(int i=0; i<np_; ++i) (*p_temp_)(i) = p_seed1[i];
    
    integrator_->setForwardSeed(1,*x_temp_,*p_temp_);
    
    // AD forward
    integrator_->integrateSensitivities();

    // Get the sensitivities
    integrator_->getForwardSensitivities( *x_temp_, 1);
    
    vector<double> &x_sens1 = output(INTEGRATOR_XF).dataF();
    for(int i=0; i<nx_; ++i) x_sens1[i] = (*x_temp_)(i);
  }

#if 0
  // Second order sensitivities
  if(fsens_order>=2){
    // Pass forward seeds
    const vector<double> &x_seed2 = input(INTEGRATOR_X0).data(2);
    for(int i=0; i<nx_; ++i) (*x_temp_)(i) = x_seed2[i];

    const vector<double> &p_seed2 = input(INTEGRATOR_P).data(2);
    for(int i=0; i<np_; ++i) (*p_temp_)(i) = p_seed2[i];

    integrator_->setForwardSeed(2,*x_temp_,*p_temp_);

    // AD forward
    integrator_->integrateSensitivities();
    integrator_->getForwardSensitivities( *x_temp_, 2);

    // Get the sensitivities
    vector<double> &x_sens2 = output(INTEGRATOR_XF).data(2);
    for(int i=0; i<nx_; ++i) x_sens2[i] = (*x_temp_)(i);
  }
#endif

}

void ACADOIntegratorInternal::integrateAdj(double t_out){
}

void ACADOIntegratorInternal::printStats(ostream &stream) const{
}


} // namespace CasADi

