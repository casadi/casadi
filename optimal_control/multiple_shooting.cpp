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

#include "multiple_shooting_internal.hpp"

namespace CasADi{
  namespace OptimalControl{
    
MultipleShooting::MultipleShooting(){
}
    
MultipleShooting::MultipleShooting(const FX& ffcn, const FX& mfcn, const FX& cfcn, const FX& rfcn){
  assignNode(new MultipleShootingInternal(ffcn,mfcn,cfcn,rfcn));
}

const MultipleShootingInternal* MultipleShooting::operator->() const{
  return (const MultipleShootingInternal*)FX::operator->();
}

MultipleShootingInternal* MultipleShooting::operator->(){
  return (MultipleShootingInternal*)FX::operator->();
}

void MultipleShooting::getGuess(std::vector<double>& V_init) const{
  (*this)->getGuess(V_init);
}
    
void MultipleShooting::getVariableBounds(std::vector<double>& V_min, std::vector<double>& V_max) const{
  (*this)->getVariableBounds(V_min,V_max);
}
    
void MultipleShooting::getConstraintBounds(std::vector<double>& G_min, std::vector<double>& G_max) const{
  (*this)->getConstraintBounds(G_min,G_max);
}

void MultipleShooting::setOptimalSolution( const std::vector<double> &V_opt ){
  (*this)->setOptimalSolution(V_opt);
}


  } // namespace OptimalControl
} // namespace CasADi

