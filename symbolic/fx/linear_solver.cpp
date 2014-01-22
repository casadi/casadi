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

#include "linear_solver_internal.hpp"

using namespace std;
namespace CasADi{

  LinearSolver::LinearSolver(){
  }

  LinearSolver::LinearSolver(const CRSSparsity& sp, int nrhs){
    assignNode(new LinearSolverInternal(sp,nrhs));
  }

  LinearSolverInternal* LinearSolver::operator->(){
    return static_cast<LinearSolverInternal*>(FX::operator->());
  }

  const LinearSolverInternal* LinearSolver::operator->() const{
    return static_cast<const LinearSolverInternal*>(FX::operator->());
  }
 
  void LinearSolver::prepare(){
    assertInit();
    (*this)->prepare();
  }

  void LinearSolver::solve(double* x, int nrhs, bool transpose){
    assertInit();
    (*this)->solve(x,nrhs,transpose);
  }
 
  void LinearSolver::solve(bool transpose){
    assertInit();
    (*this)->solve(transpose);
  }

  MX LinearSolver::solve(const MX& A, const MX& B, bool transpose){
    assertInit();
    return (*this)->solve(A,B,transpose);
  }
 
  bool LinearSolver::prepared() const{
    assertInit();
    return (*this)->prepared_;
  }
 
  bool LinearSolver::checkNode() const{
    return dynamic_cast<const LinearSolverInternal*>(get())!=0;
  }

  void LinearSolver::spSolve(bvec_t* X, bvec_t* B, bool transpose) const{
    (*this)->spSolve(X,B,transpose);
  }


} // namespace CasADi

  


