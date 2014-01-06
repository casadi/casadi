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
    return static_cast<IntegratorInternal*>(GenericIntegrator::operator->());
  }

  const IntegratorInternal* Integrator::operator->() const{
    return static_cast<const IntegratorInternal*>(GenericIntegrator::operator->()); 
  }
  
  void Integrator::reset(){
    (*this)->reset();
  }

  void Integrator::integrate(double t_out){
    (*this)->integrate(t_out);
  }
    
  bool Integrator::checkNode() const{
    return dynamic_cast<const IntegratorInternal*>(get())!=0;
  }

  void Integrator::resetB(){
    (*this)->resetB();
  }

  void Integrator::integrateB(double t_out){
    (*this)->integrateB(t_out);
  }

  FX Integrator::getDAE(){
    return (*this)->f_;
  }

  std::pair<FX,FX> Integrator::getAugmented(int nfwd, int nadj){
    IntegratorInternal::AugOffset offset;
    return (*this)->getAugmented(nfwd,nadj,offset);
  }

  DMatrix& Integrator::z(){
    return (*this)->z_;
  }

  const DMatrix& Integrator::z() const{
    return (*this)->z_;
  }

  DMatrix& Integrator::rz(){
    return (*this)->rz_;
  }

  const DMatrix& Integrator::rz() const{
    return (*this)->rz_;
  }
 
} // namespace CasADi

