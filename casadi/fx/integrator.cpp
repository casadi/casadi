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

#include "integrator.hpp"
#include "integrator_internal.hpp"
#include <cassert>

namespace CasADi{

  
  
Integrator::Integrator(){
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
  
void Integrator::assertNode() const{
  if(!dynamic_cast<const IntegratorInternal*>(get()))
    throw CasadiException("Integrator::assertNode");
}
  
void Integrator::setLinearSolver(const LinearSolver& prototype){
  (*this)->setLinearSolver(prototype.getCreator());
}

 
} // namespace CasADi

