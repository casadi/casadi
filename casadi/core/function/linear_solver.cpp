/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010-2014 Joel Andersson, Joris Gillis, Moritz Diehl,
 *                            K.U. Leuven. All rights reserved.
 *    Copyright (C) 2011-2014 Greg Horn
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
namespace casadi {

  LinearSolver::LinearSolver() {
  }

  LinearSolver::LinearSolver(const Sparsity& sp, int nrhs) {
    assignNode(new LinearSolverInternal(sp, nrhs));
  }

  LinearSolverInternal* LinearSolver::operator->() {
    return static_cast<LinearSolverInternal*>(Function::operator->());
  }

  const LinearSolverInternal* LinearSolver::operator->() const {
    return static_cast<const LinearSolverInternal*>(Function::operator->());
  }

  void LinearSolver::prepare() {
    assertInit();
    (*this)->prepare();
  }

  void LinearSolver::solve(double* x, int nrhs, bool transpose) {
    assertInit();
    (*this)->solve(x, nrhs, transpose);
  }

  void LinearSolver::solve(bool transpose) {
    assertInit();
    (*this)->solve(transpose);
  }

  MX LinearSolver::solve(const MX& A, const MX& B, bool transpose) {
    assertInit();
    return (*this)->solve(A, B, transpose);
  }

  bool LinearSolver::prepared() const {
    assertInit();
    return (*this)->prepared_;
  }

  bool LinearSolver::testCast(const SharedObjectNode* ptr) {
    return dynamic_cast<const LinearSolverInternal*>(ptr)!=0;
  }

  void LinearSolver::spSolve(bvec_t* X, const bvec_t* B, bool transpose) const {
    (*this)->spSolve(X, B, transpose);
  }

  void LinearSolver::spSolve(DMatrix& X, const DMatrix& B, bool transpose) const {
    (*this)->spSolve(X, B, transpose);
  }

  LinearSolver::LinearSolver(const std::string& name, const Sparsity& sp, int nrhs) {
    assignNode(LinearSolverInternal::getPlugin(name).creator(sp, nrhs));
  }

  bool LinearSolver::hasPlugin(const std::string& name) {
    return LinearSolverInternal::hasPlugin(name);
  }

  void LinearSolver::loadPlugin(const std::string& name) {
    LinearSolverInternal::loadPlugin(name);
  }

  std::string LinearSolver::doc(const std::string& name) {
    return LinearSolverInternal::getPlugin(name).doc;
  }

  Sparsity LinearSolver::getFactorizationSparsity(bool transpose) const {
    return (*this)->getFactorizationSparsity(transpose);
  }

  DMatrix LinearSolver::getFactorization(bool transpose) const {
    return (*this)->getFactorization(transpose);
  }

  void LinearSolver::solveL(double* x, int nrhs, bool transpose) {
    return (*this)->solveL(x, nrhs, transpose);
  }

} // namespace casadi
