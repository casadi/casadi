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


#include "stabilized_qp_solver_internal.hpp"

using namespace std;
namespace casadi {

  StabilizedQpSolver::StabilizedQpSolver() {
  }

  StabilizedQpSolverInternal* StabilizedQpSolver::operator->() {
    return static_cast<StabilizedQpSolverInternal*>(Function::operator->());
  }

  const StabilizedQpSolverInternal* StabilizedQpSolver::operator->() const {
    return static_cast<const StabilizedQpSolverInternal*>(Function::operator->());
  }

  bool StabilizedQpSolver::testCast(const SharedObjectNode* ptr) {
    return dynamic_cast<const StabilizedQpSolverInternal*>(ptr)!=0;
  }

  void StabilizedQpSolver::setLPOptions() {
    (*this)->setLPOptions();
  }

  void StabilizedQpSolver::generateNativeCode(const std::string &filename) const {
    std::ofstream file;
    file.open(filename.c_str());
    (*this)->generateNativeCode(file);
  }

  StabilizedQpSolver::StabilizedQpSolver(const std::string& name, const QPStructure& st) {
    assignNode(StabilizedQpSolverInternal::instantiatePlugin(name, st));
  }

  bool StabilizedQpSolver::hasPlugin(const std::string& name) {
    return StabilizedQpSolverInternal::hasPlugin(name);
  }

  void StabilizedQpSolver::loadPlugin(const std::string& name) {
    StabilizedQpSolverInternal::loadPlugin(name);
  }

  std::string StabilizedQpSolver::doc(const std::string& name) {
    return StabilizedQpSolverInternal::getPlugin(name).doc;
  }

} // namespace casadi
