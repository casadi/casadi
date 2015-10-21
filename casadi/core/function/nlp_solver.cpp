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


#include "nlp_solver_internal.hpp"

using namespace std;
namespace casadi {

  NlpSolver::NlpSolver() {
  }

  NlpSolver::NlpSolver(const std::string& name, const std::string& solver,
                       const Function& f, const Dict& opts) {
    Function tmp = Function::nlp_solver(name, solver, f, opts);
    assignNode(tmp.get());
  }

  NlpSolverInternal* NlpSolver::operator->() {
    return static_cast<NlpSolverInternal*>(Function::operator->());
  }

  const NlpSolverInternal* NlpSolver::operator->() const {
    return static_cast<const NlpSolverInternal*>(Function::operator->());
  }

  bool NlpSolver::test_cast(const SharedObjectNode* ptr) {
    return dynamic_cast<const NlpSolverInternal*>(ptr)!=0;
  }

  void NlpSolver::reportConstraints(std::ostream &stream) {
    (*this)->reportConstraints();
  }

  Function NlpSolver::nlp() {
    return (*this)->nlp_;
  }

  Function NlpSolver::gradF() {
    return (*this)->gradF();
  }

  Function NlpSolver::jacG() {
    return (*this)->jacG();
  }

  Function NlpSolver::hessLag() {
    return (*this)->hessLag();
  }

  bool NlpSolver::hasPlugin(const std::string& name) {
    return NlpSolverInternal::hasPlugin(name);
  }

  void NlpSolver::loadPlugin(const std::string& name) {
    NlpSolverInternal::loadPlugin(name);
  }

  std::string NlpSolver::doc(const std::string& name) {
    return NlpSolverInternal::getPlugin(name).doc;
  }

  DMatrix NlpSolver::getReducedHessian() {
    return (*this)->getReducedHessian();
  }

  void NlpSolver::setOptionsFromFile(const std::string & file) {
    (*this)->setOptionsFromFile(file);
  }

} // namespace casadi
