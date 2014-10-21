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


#include "qp_solver_internal.hpp"

using namespace std;
namespace casadi {


  QpSolver::QpSolver() {
  }

  QpSolverInternal* QpSolver::operator->() {
    return static_cast<QpSolverInternal*>(Function::operator->());
  }

  const QpSolverInternal* QpSolver::operator->() const {
    return static_cast<const QpSolverInternal*>(Function::operator->());
  }

  bool QpSolver::testCast(const SharedObjectNode* ptr) {
    return dynamic_cast<const QpSolverInternal*>(ptr)!=0;
  }

  void QpSolver::setLPOptions() {
    (*this)->setLPOptions();
  }

  void QpSolver::generateNativeCode(const std::string &filename) const {
    std::ofstream file;
    file.open(filename.c_str());
    (*this)->generateNativeCode(file);
  }

  void QpSolver::generateNativeCode(std::ostream &file) const {
    (*this)->generateNativeCode(file);
  }

  QpSolver::QpSolver(const std::string& name, const QPStructure& st) {
    assignNode(QpSolverInternal::instantiatePlugin(name, st));
  }

  bool QpSolver::hasPlugin(const std::string& name) {
    return QpSolverInternal::hasPlugin(name);
  }

  void QpSolver::loadPlugin(const std::string& name) {
    QpSolverInternal::loadPlugin(name);
  }

  std::string QpSolver::doc(const std::string& name) {
    return QpSolverInternal::getPlugin(name).doc;
  }

} // namespace casadi
