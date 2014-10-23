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


#include "lr_dle_solver.hpp"
#include "lr_dle_internal.hpp"

using namespace std;
namespace casadi {

  LrDleSolver::LrDleSolver() {

  }

  LrDleSolver  LrDleSolver::clone() const {
    LrDleSolver ret;
    if (!isNull()) ret.assignNode((*this)->clone());
    return ret;
  }

  void LrDleSolver::printStats(ostream &stream) const {
    (*this)->printStats(stream);
  }

  LrDleInternal* LrDleSolver::operator->() {
    return static_cast<LrDleInternal*>(Function::operator->());
  }

  const LrDleInternal* LrDleSolver::operator->() const {
    return static_cast<const LrDleInternal*>(Function::operator->());
  }

  bool LrDleSolver::testCast(const SharedObjectNode* ptr) {
    return dynamic_cast<const LrDleInternal*>(ptr)!=0;
  }

  bool LrDleSolver::hasPlugin(const std::string& name) {
    return LrDleInternal::hasPlugin(name);
  }

  void LrDleSolver::loadPlugin(const std::string& name) {
    LrDleInternal::loadPlugin(name);
  }

  std::string LrDleSolver::doc(const std::string& name) {
    return LrDleInternal::getPlugin(name).doc;
  }

  LrDleSolver::LrDleSolver(const std::string& name,
                         const LrDleStructure& st, const std::vector<int> &Hs) {
    assignNode(LrDleInternal::getPlugin(name).creator(st, Hs));
  }

  Sparsity LrDleSolver::getSparsity(const LrDleStructure& st, const std::vector<int> &Hs) {
    return LrDleInternal::getSparsity(st, Hs);
  }

} // namespace casadi

