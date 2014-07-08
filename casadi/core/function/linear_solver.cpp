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

  bool LinearSolver::checkNode() const {
    return dynamic_cast<const LinearSolverInternal*>(get())!=0;
  }

  void LinearSolver::spSolve(bvec_t* X, const bvec_t* B, bool transpose) const {
    (*this)->spSolve(X, B, transpose);
  }

  void LinearSolver::spSolve(DMatrix& X, const DMatrix& B, bool transpose) const {
    (*this)->spSolve(X, B, transpose);
  }

  LinearSolver::LinearSolver(const std::string& name, const Sparsity& sp, int nrhs) {
    // Check if the solver has been loaded
    std::map<std::string, Plugin>::iterator it=solvers_.find(name);

    // Load the solver if needed
    if (it==solvers_.end()) {
      loadPlugin(name);
      it=solvers_.find(name);
    }
    casadi_assert(it!=solvers_.end());
    assignNode(it->second.creator(sp, nrhs));
  }

  std::map<std::string, LinearSolver::Plugin> LinearSolver::solvers_;

  void LinearSolver::registerPlugin(RegFcn regfcn) {
    // Create a temporary struct
    Plugin plugin;
   
    // Set the fields
    int flag = regfcn(&plugin);
    casadi_assert(flag==0);

    // Check if the solver name is in use
    std::map<std::string, Plugin>::iterator it=solvers_.find(plugin.name);
    casadi_assert_message(it==solvers_.end(), "Solver " << plugin.name << " is already in use");

    // Add to list of solvers
    solvers_[plugin.name] = plugin;
  }

  void LinearSolver::loadPlugin(const std::string& name) {
#ifndef WITH_DL
    casadi_error("WITH_DL option needed for dynamic loading");
#else // WITH_DL
    // Retrieve the registration function
    RegFcn reg = FunctionInternal::loadPlugin<RegFcn>(name,"linearsolver");

    // Register the plugin
    registerPlugin(reg);
#endif // WITH_DL
  }


} // namespace casadi




