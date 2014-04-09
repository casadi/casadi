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

#include "direct_single_shooting_internal.hpp"

namespace casadi{

DirectSingleShooting::DirectSingleShooting(){
}

DirectSingleShooting::DirectSingleShooting(const Function& ffcn, const Function& mfcn, const Function& cfcn, const Function& rfcn){
  assignNode(new DirectSingleShootingInternal(ffcn,mfcn,cfcn,rfcn));
}

const DirectSingleShootingInternal* DirectSingleShooting::operator->() const{
  return static_cast<const DirectSingleShootingInternal*>(Function::operator->());
}

DirectSingleShootingInternal* DirectSingleShooting::operator->(){
  return static_cast<DirectSingleShootingInternal*>(Function::operator->());
}

void DirectSingleShooting::getGuess(std::vector<double>& V_init) const{
  (*this)->getGuess(V_init);
}

void DirectSingleShooting::getVariableBounds(std::vector<double>& V_min, std::vector<double>& V_max) const{
  (*this)->getVariableBounds(V_min,V_max);
}

void DirectSingleShooting::getConstraintBounds(std::vector<double>& G_min, std::vector<double>& G_max) const{
  (*this)->getConstraintBounds(G_min,G_max);
}

void DirectSingleShooting::setOptimalSolution( const std::vector<double> &V_opt ){
  (*this)->setOptimalSolution(V_opt);
}

  NLPSolver DirectSingleShooting::getNLPSolver() const { return isNull() ? NLPSolver(): (*this)->nlp_solver_; }

void DirectSingleShooting::reportConstraints(std::ostream &stream) {
  (*this)->reportConstraints();
}

} // namespace casadi

