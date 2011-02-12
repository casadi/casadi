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

FX MultipleShooting::getF() const{
  return (*this)->F_;
}
    
FX MultipleShooting::getG() const{
  return (*this)->G_;
}
    
FX MultipleShooting::getJ() const{
  return (*this)->J_;
}

void MultipleShooting::setNLPSolver(const NLPSolver& nlp_solver){
  (*this)->nlp_solver_ = nlp_solver;
}

const MultipleShootingInternal* MultipleShooting::operator->() const{
  return (const MultipleShootingInternal*)FX::operator->();
}

MultipleShootingInternal* MultipleShooting::operator->(){
  return (MultipleShootingInternal*)FX::operator->();
}

  } // namespace OptimalControl
} // namespace CasADi

