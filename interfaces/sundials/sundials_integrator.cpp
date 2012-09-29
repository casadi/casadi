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

#include "sundials_integrator.hpp"
#include "sundials_internal.hpp"
#include <cassert>

using namespace std;
namespace CasADi{

SundialsIntegrator::SundialsIntegrator(){
}

SundialsInternal* SundialsIntegrator::operator->(){
  return (SundialsInternal*)(FX::operator->());
}

const SundialsInternal* SundialsIntegrator::operator->() const{
   return (const SundialsInternal*)(FX::operator->()); 
}
  
bool SundialsIntegrator::checkNode() const{
  return dynamic_cast<const SundialsInternal*>(get())!=0;
}

void SundialsIntegrator::setLinearSolver(const LinearSolver& linsol, const FX& jac){
  casadi_warning(
    "Depreciated function \"Integrator::setLinearSolver\",\n"
    "use setOption(\"linear solver_creator\",SolverName::creator) in C++ \n"
    "or setOption(\"linear solver_creator\",SolverName) in Python/Octave instead.\n"
    "Options to the linear solver are passed with setOption(\"linear solver_options\",...)\n"
    "This function will be removed in the next release"
  );
  (*this)->setLinearSolver(linsol,jac);
}

FX SundialsIntegrator::getJacobian(){
  return (*this)->getJacobian();  
}
  
LinearSolver SundialsIntegrator::getLinearSolver(){
  return (*this)->getLinearSolver();  
}

void SundialsIntegrator::setStopTime(double tf){
  (*this)->setStopTime(tf);
}
  
void SundialsIntegrator::setInitialTime(double t0){
  (*this)->setInitialTime(t0);
}

void SundialsIntegrator::setFinalTime(double tf){
  (*this)->setFinalTime(tf);
}

} // namespace CasADi

