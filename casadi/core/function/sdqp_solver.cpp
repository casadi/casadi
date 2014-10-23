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


#include "sdqp_solver_internal.hpp"

using namespace std;
namespace casadi {

  SdqpSolver::SdqpSolver() {
  }

  SdqpSolverInternal* SdqpSolver::operator->() {
    return static_cast<SdqpSolverInternal*>(Function::operator->());
  }

  const SdqpSolverInternal* SdqpSolver::operator->() const {
    return static_cast<const SdqpSolverInternal*>(Function::operator->());
  }

  bool SdqpSolver::testCast(const SharedObjectNode* ptr) {
    return dynamic_cast<const SdqpSolverInternal*>(ptr)!=0;
  }

  void SdqpSolver::setSOCQPOptions() {
    (*this)->setSOCQPOptions();
  }

  SdqpSolver::SdqpSolver(const std::string& name, const SDQPStructure& st) {
    assignNode(SdqpSolverInternal::instantiatePlugin(name, st));
  }

  bool SdqpSolver::hasPlugin(const std::string& name) {
    return SdqpSolverInternal::hasPlugin(name);
  }

  void SdqpSolver::loadPlugin(const std::string& name) {
    SdqpSolverInternal::loadPlugin(name);
  }

  std::string SdqpSolver::doc(const std::string& name) {
    return SdqpSolverInternal::getPlugin(name).doc;
  }

} // namespace casadi
