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


#include "dle_solver.hpp"
#include "dle_internal.hpp"

using namespace std;
namespace casadi {

  DleSolver::DleSolver() {

  }

  DleSolver  DleSolver::clone() const {
    DleSolver ret;
    if (!isNull()) ret.assignNode((*this)->clone());
    return ret;
  }

  void DleSolver::printStats(ostream &stream) const {
    (*this)->printStats(stream);
  }

  DleInternal* DleSolver::operator->() {
    return static_cast<DleInternal*>(Function::operator->());
  }

  const DleInternal* DleSolver::operator->() const {
    return static_cast<const DleInternal*>(Function::operator->());
  }

  bool DleSolver::testCast(const SharedObjectNode* ptr) {
    return dynamic_cast<const DleInternal*>(ptr)!=0;
  }

  bool DleSolver::hasPlugin(const std::string& name) {
    return DleInternal::hasPlugin(name);
  }

  void DleSolver::loadPlugin(const std::string& name) {
    DleInternal::loadPlugin(name);
  }

  std::string DleSolver::doc(const std::string& name) {
    return DleInternal::getPlugin(name).doc;
  }

  DleSolver::DleSolver(const std::string& name,
                         const DleStructure& st) {
    assignNode(DleInternal::instantiatePlugin(name, st));
  }

  Sparsity DleSolver::getSparsity(const DleStructure& st) {
    return DleInternal::getSparsity(st);
  }

} // namespace casadi

