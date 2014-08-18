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

#include "lp_internal.hpp"

using namespace std;
namespace casadi {

  LpSolver::LpSolver() {
  }

  LpSolverInternal* LpSolver::operator->() {
    return static_cast<LpSolverInternal*>(Function::operator->());
  }

  const LpSolverInternal* LpSolver::operator->() const {
    return static_cast<const LpSolverInternal*>(Function::operator->());
  }

  bool LpSolver::testCast(const SharedObjectNode* ptr) {
    return dynamic_cast<const LpSolverInternal*>(ptr)!=0;
  }

  LpSolver::LpSolver(const std::string& name, const LPStructure& st) {
    assignNode(LpSolverInternal::getPlugin(name).creator(st));
  }

  void LpSolver::loadPlugin(const std::string& name) {
    LpSolverInternal::loadPlugin(name);
  }

  std::string LpSolver::doc(const std::string& name) {
    return LpSolverInternal::getPlugin(name).doc;
  }

} // namespace casadi
