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

#include "kinsol_solver.hpp"
#include "kinsol_internal.hpp"

using namespace std;
namespace CasADi{

  KinsolSolver::KinsolSolver(){ 
  }

  KinsolSolver::KinsolSolver(const FX& f, const FX& jac, const LinearSolver& linsol){
    assignNode(new KinsolInternal(f,jac,linsol));
  }

  KinsolInternal* KinsolSolver::operator->(){
    return (KinsolInternal*)(FX::operator->());
  }

  const KinsolInternal* KinsolSolver::operator->() const{
    return (const KinsolInternal*)(FX::operator->());
  }

  bool KinsolSolver::checkNode() const{
    return dynamic_cast<const KinsolInternal*>(get());
  }

  void KinsolSolver::setLinearSolver(const LinearSolver& linsol){
    casadi_warning(
		   "Depreciated function \"KinsolSolver::setLinearSolver\",\n"
		   "use setOption(\"linear solver_creator\",SolverName::creator) in C++ \n"
		   "or setOption(\"linear solver_creator\",SolverName) in Python/Octave instead.\n"
		   "Options to the linear solver are passed with setOption(\"linear solver_options\",...)\n"
		   "This function will be removed in the next release"
		   );
    (*this)->setLinearSolver(linsol);
  }

  void KinsolSolver::setJacobian(const FX& jac){
    (*this)->setJacobian(jac);
  }

  FX KinsolSolver::getJacobian(){
    return (*this)->getJacobian();  
  }
  
  LinearSolver KinsolSolver::getLinearSolver(){
    return (*this)->getLinearSolver();  
  }

} // namespace CasADi


