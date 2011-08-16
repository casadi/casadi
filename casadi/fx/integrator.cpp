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

#include "integrator.hpp"
#include "integrator_internal.hpp"
#include <cassert>

using namespace std;
namespace CasADi{

Integrator::Integrator(){
}

Integrator  Integrator::clone() const{
  Integrator ret;
  if(!isNull()) ret.assignNode((*this)->clone());
  return ret;
}

void Integrator::printStats(ostream &stream) const{
  (*this)->printStats(stream);
}
  
IntegratorInternal* Integrator::operator->(){
  return (IntegratorInternal*)(FX::operator->());
}

const IntegratorInternal* Integrator::operator->() const{
   return (const IntegratorInternal*)(FX::operator->()); 
}
  
void Integrator::reset(int fsens_order, int asens_order){
  (*this)->reset(fsens_order, asens_order);
}

void Integrator::integrate(double t_out){
  (*this)->integrate(t_out);
}
  
void Integrator::setStopTime(double tf){
  (*this)->setStopTime(tf);
}
  
void Integrator::setInitialTime(double t0){
  (*this)->setInitialTime(t0);
}

void Integrator::setFinalTime(double tf){
  (*this)->setFinalTime(tf);
}
  
bool Integrator::checkNode() const{
  return dynamic_cast<const IntegratorInternal*>(get())!=0;
}

void Integrator::setLinearSolver(const LinearSolver& linsol, const FX& jac){
  casadi_warning(
    "Depreciated function \"Integrator::setLinearSolver\",\n"
    "use setOption(\"linear solver_creator\",SolverName::creator) in C++ \n"
    "or setOption(\"linear solver_creator\",SolverName) in Python/Octave instead.\n"
    "Options to the linear solver are passed with setOption(\"linear solver_options\",...)\n"
    "This function will be removed in the next release"
  );
  (*this)->setLinearSolver(linsol,jac);
}

void Integrator::resetAdj(){
  (*this)->resetAdj();
}

void Integrator::integrateAdj(double t_out){
  (*this)->integrateAdj(t_out);
}

FX Integrator::getJacobian(){
  return (*this)->getJacobian();  
}
  
LinearSolver Integrator::getLinearSolver(){
  return (*this)->getLinearSolver();  
}

FX Integrator::getDAE(){
  return (*this)->f_;
}
  
 
} // namespace CasADi

