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


#include "sdp_solver_internal.hpp"

using namespace std;
namespace casadi {

  SdpSolver::SdpSolver() {
  }

  SdpSolverInternal* SdpSolver::operator->() {
    return static_cast<SdpSolverInternal*>(Function::operator->());
  }

  const SdpSolverInternal* SdpSolver::operator->() const {
    return static_cast<const SdpSolverInternal*>(Function::operator->());
  }

  bool SdpSolver::testCast(const SharedObjectNode* ptr) {
    return dynamic_cast<const SdpSolverInternal*>(ptr)!=0;
  }

  void SdpSolver::setSOCPOptions() {
    (*this)->setSOCPOptions();
  }

  SdpSolver::SdpSolver(const std::string& name, const SDPStructure& st) {
    assignNode(SdpSolverInternal::instantiatePlugin(name, st));
  }

  bool SdpSolver::hasPlugin(const std::string& name) {
    return SdpSolverInternal::hasPlugin(name);
  }

  void SdpSolver::loadPlugin(const std::string& name) {
    SdpSolverInternal::loadPlugin(name);
  }

  std::string SdpSolver::doc(const std::string& name) {
    return SdpSolverInternal::getPlugin(name).doc;
  }

} // namespace casadi
