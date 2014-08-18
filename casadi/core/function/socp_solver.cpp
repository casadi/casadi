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

#include "socp_solver_internal.hpp"

using namespace std;
namespace casadi {

  SocpSolver::SocpSolver() {
  }

  SocpSolverInternal* SocpSolver::operator->() {
    return static_cast<SocpSolverInternal*>(Function::operator->());
  }

  const SocpSolverInternal* SocpSolver::operator->() const {
    return static_cast<const SocpSolverInternal*>(Function::operator->());
  }

  bool SocpSolver::testCast(const SharedObjectNode* ptr) {
    return dynamic_cast<const SocpSolverInternal*>(ptr)!=0;
  }

  SocpSolver::SocpSolver(const std::string& name, const SOCPStructure& st) {
    assignNode(SocpSolverInternal::getPlugin(name).creator(st));
  }

  void SocpSolver::loadPlugin(const std::string& name) {
    SocpSolverInternal::loadPlugin(name);
  }

  std::string SocpSolver::doc(const std::string& name) {
    return SocpSolverInternal::getPlugin(name).doc;
  }

} // namespace casadi
