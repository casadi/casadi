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

#include <acado/utils/acado_utils.hpp>
#include <acado/function/function_.hpp>
#include "acado_integrator_backend.hpp"
#include <cassert>
#include "casadi/casadi_exception.hpp"
#include "casadi/stl_vector_tools.hpp"

using namespace std;
using namespace ACADO;

namespace CasADi{

ACADO::Integrator* AcadoIntegratorBackend::create(void *user_data){
  return new AcadoIntegratorBackend(user_data);
}

AcadoIntegratorBackend::AcadoIntegratorBackend(void *user_data){
   ocp_solver_ = (CasADi::AcadoInternal*)user_data;
   my_ref_ = ocp_solver_->getRef(this);
   assert(ocp_solver_->integrators_.size()>my_ref_);
   integrator_ = ocp_solver_->integrators_[my_ref_];
}


AcadoIntegratorBackend::AcadoIntegratorBackend( const ACADO::DifferentialEquation& rhs__ ){
  assert(0);
}

AcadoIntegratorBackend::AcadoIntegratorBackend( const AcadoIntegratorBackend& arg ){
  ocp_solver_ = arg.ocp_solver_;
  my_ref_ = ocp_solver_->getRef(this);
  assert(ocp_solver_->integrators_.size()>my_ref_);
  integrator_ = ocp_solver_->integrators_[my_ref_];
  
  // The following is necessary for the base class (WHY?)
  if(rhs){
    delete rhs;
    rhs = 0;
  }
  if(arg.rhs){
    rhs = new ACADO::DifferentialEquation( *arg.rhs );
  }
  
  // The following are base class member variables that should have been copied by the base class copy constructor
  m   = arg.m;
  ma  = arg.ma;
  mdx = arg.mdx;
  mn  = arg.mn;
  mu  = arg.mu;
  mui = arg.mui;
  mp  = arg.mp;
  mpi = arg.mpi;
  mw  = arg.mw;
  md  = m-ma;
}



AcadoIntegratorBackend::~AcadoIntegratorBackend( ){
  ocp_solver_->returnRef(my_ref_);
}


AcadoIntegratorBackend& AcadoIntegratorBackend::operator=( const AcadoIntegratorBackend& arg ){
  assert(0);
}

ACADO::returnValue AcadoIntegratorBackend::init( const ACADO::DifferentialEquation &rhs__ ){
  cout << "init " << this << endl;
  
  // The following are base class member variables that should have been initialized by the base class copy constructor
  m   = rhs__.getDim ();
  ma  = rhs__.getNXA ();
  mdx = rhs__.getNDX ();
  mn  = rhs__.getN   ();
  mu  = rhs__.getNU  ();
  mui = rhs__.getNUI ();
  mp  = rhs__.getNP  ();
  mpi = rhs__.getNPI ();
  mw  = rhs__.getNW  ();
  md  = rhs__.getNumDynamicEquations();

  // The following is necessary for the base class (WHY?)
  if(rhs) delete rhs;
  rhs = new ACADO::DifferentialEquation( rhs__ );
  if( m < md+ma ) md = m - ma;    
  return SUCCESSFUL_RETURN;
}


ACADO::Integrator* AcadoIntegratorBackend::clone() const{
  return new AcadoIntegratorBackend(*this);
}


ACADO::returnValue AcadoIntegratorBackend::freezeMesh(){
  soa = SOA_FREEZING_MESH;
  return SUCCESSFUL_RETURN;
}


ACADO::returnValue AcadoIntegratorBackend::freezeAll(){
  soa = SOA_FREEZING_ALL;
  return SUCCESSFUL_RETURN;
}


ACADO::returnValue AcadoIntegratorBackend::unfreeze(){
  soa = SOA_UNFROZEN;
  return SUCCESSFUL_RETURN;
}


ACADO::returnValue AcadoIntegratorBackend::evaluate( const Vector &x0  ,
                                     const Vector &xa  ,
                                     const Vector &p   ,
                                     const Vector &u   ,
                                     const Vector &w   ,
                                     const Grid   &t_    ){
  
  // Set the time horizon
  integrator_.setInput(t_.getFirstTime(),CasADi::INTEGRATOR_T0);
  integrator_.setInput(t_.getLastTime(),CasADi::INTEGRATOR_TF);
  
  // Set the initial conditions
  Matrix<double>& yy = integrator_.input(CasADi::INTEGRATOR_X0);
  &yy[0] << const_cast<Vector &>(x0);
//  &yy[md] << const_cast<Vector &>(xa);
  

  // Parameters and controls
  Matrix<double>& pp = integrator_.input(CasADi::INTEGRATOR_P);
  &pp[0] << const_cast<Vector &>(p);
  &pp[mp] << const_cast<Vector &>(u);
  
  // Integrate
  
  bool with_sens = (soa == SOA_FREEZING_ALL);
  
  if(with_sens){
    for(int i=0; i<md; ++i){
      Matrix<double> &xsens = integrator_.fwdSeed(CasADi::INTEGRATOR_X0,i);
      xsens[i] = 1;
    }
    for(int i=0; i<mu+mp; ++i){
      Matrix<double> &usens = integrator_.fwdSeed(CasADi::INTEGRATOR_P,i+md);
      usens[i] = 1;
    }
  }
  
  integrator_.evaluate(with_sens,0);
    
  
  Integrator::initializeOptions();

  timeInterval = t_;
  xStore.init( md+ma, timeInterval );
  iStore.init( mn   , timeInterval );
  return SUCCESSFUL_RETURN;
}


ACADO::returnValue AcadoIntegratorBackend::setProtectedForwardSeed( const Vector &xSeed , const Vector &pSeed ,
                                                    const Vector &uSeed , const Vector &wSeed ,
                                                    const int    &order   ){
  
  if(xSeed.getDim()>0)
    for(int i=0; i<md; ++i)
      if(xSeed(i)==1)
        ider_ = i;

  if(pSeed.getDim()>0)
    for(int i=0; i<mp; ++i)
      if(pSeed(i)==1)
        ider_ = md+i;

  if(uSeed.getDim()>0)
    for(int i=0; i<mu; ++i)
      if(uSeed(i)==1)
        ider_ = md+mp+i;

    return SUCCESSFUL_RETURN;
}


ACADO::returnValue AcadoIntegratorBackend::setProtectedBackwardSeed( const Vector &seed, const int &order ){
  assert(0);
  return ACADOERROR(RET_NOT_IMPLEMENTED_YET);
}


ACADO::returnValue AcadoIntegratorBackend::setBackwardSeed2( const Vector &seed ){
  assert(0);
  return ACADOERROR(RET_NOT_IMPLEMENTED_YET);
}


ACADO::returnValue AcadoIntegratorBackend::evaluateSensitivities(){
  return SUCCESSFUL_RETURN;
}


ACADO::returnValue AcadoIntegratorBackend::step(int number_){
  assert(0);
  return ACADOERROR(RET_NOT_IMPLEMENTED_YET);
}



ACADO::returnValue AcadoIntegratorBackend::stop(){
  assert(0);
  return ACADOERROR(RET_NOT_IMPLEMENTED_YET);
}



ACADO::returnValue AcadoIntegratorBackend::getProtectedX( Vector *xEnd ) const{
  const Matrix<double>& xf = integrator_.output(CasADi::INTEGRATOR_XF);
  if( (int) xEnd[0].getDim() != xf.size() )
    return RET_INPUT_HAS_WRONG_DIMENSION;

  for(int i=0; i<xf.size(); ++i)
    xEnd[0](i) = xf[i];

  return SUCCESSFUL_RETURN;
}


returnValue AcadoIntegratorBackend::getProtectedForwardSensitivities( ACADO::Matrix *Dx, int order ) const{
  const Matrix<double> &sens = integrator_.fwdSens(CasADi::INTEGRATOR_XF,ider_);

  if( Dx == NULL ){
    return SUCCESSFUL_RETURN;
  }
  
  if( order == 1 && nFDirs2 == 0 ){
    for( int run1 = 0; run1 < m; run1++ ){
            Dx[0](run1,0) = sens[run1];
    }
  }

    return SUCCESSFUL_RETURN;
}


ACADO::returnValue AcadoIntegratorBackend:: getProtectedBackwardSensitivities(Vector &Dx_x0, Vector &Dx_p , Vector &Dx_u , Vector &Dx_w , int order ) const{
  assert(0);
  return SUCCESSFUL_RETURN;
}



ACADO::returnValue AcadoIntegratorBackend::setDxInitialization( double *dx0 ){
  assert(0);
  return SUCCESSFUL_RETURN;
}


int AcadoIntegratorBackend::getNumberOfSteps() const{
    return 999;
}

int AcadoIntegratorBackend::getNumberOfRejectedSteps() const{
    return 0;
}


double AcadoIntegratorBackend::getStepSize() const{
    return h[0];
}

int AcadoIntegratorBackend::getDim() const{
    return md+ma;
}

int AcadoIntegratorBackend::getDimX() const{
    return md;
}


} // namespace ACADO

// end of file.
