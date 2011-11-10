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

#include "collocation_integrator_internal.hpp"
#include "casadi/stl_vector_tools.hpp"
#include "casadi/matrix/sparsity_tools.hpp"
#include "casadi/matrix/matrix_tools.hpp"

using namespace std;
namespace CasADi{

CollocationIntegratorInternal::CollocationIntegratorInternal(const FX& f, const FX& q) : IntegratorInternal(f,q){
}

void CollocationIntegratorInternal::deepCopyMembers(std::map<SharedObjectNode*,SharedObject>& already_copied){
}

CollocationIntegratorInternal::~CollocationIntegratorInternal(){
}

void CollocationIntegratorInternal::init(){
}
  
void CollocationIntegratorInternal::initAdj(){
}

void CollocationIntegratorInternal::reset(int fsens_order, int asens_order){
}

void CollocationIntegratorInternal::resetAdj(){
}

void CollocationIntegratorInternal::integrate(double t_out){
}

void CollocationIntegratorInternal::integrateAdj(double t_out){
}

FX CollocationIntegratorInternal::getJacobian(){
  return FX();
}
  
LinearSolver CollocationIntegratorInternal::getLinearSolver(){
  return LinearSolver();
}

void CollocationIntegratorInternal::setLinearSolver(const LinearSolver& linsol, const FX& jac){
}

void CollocationIntegratorInternal::printStats(std::ostream &stream) const{
}

void CollocationIntegratorInternal::setStopTime(double tf){
}

} // namespace CasADi
