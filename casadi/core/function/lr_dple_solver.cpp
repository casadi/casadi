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


#include "lr_dple_solver.hpp"
#include "lr_dple_internal.hpp"

using namespace std;
namespace casadi {

  LrDpleSolver::LrDpleSolver() {

  }

  LrDpleSolver  LrDpleSolver::clone() const {
    LrDpleSolver ret;
    if (!isNull()) ret.assignNode((*this)->clone());
    return ret;
  }

  void LrDpleSolver::printStats(ostream &stream) const {
    (*this)->printStats(stream);
  }

  LrDpleInternal* LrDpleSolver::operator->() {
    return static_cast<LrDpleInternal*>(Function::operator->());
  }

  const LrDpleInternal* LrDpleSolver::operator->() const {
    return static_cast<const LrDpleInternal*>(Function::operator->());
  }

  bool LrDpleSolver::testCast(const SharedObjectNode* ptr) {
    return dynamic_cast<const LrDpleInternal*>(ptr)!=0;
  }

  bool LrDpleSolver::hasPlugin(const std::string& name) {
    return LrDpleInternal::hasPlugin(name);
  }

  void LrDpleSolver::loadPlugin(const std::string& name) {
    LrDpleInternal::loadPlugin(name);
  }

  std::string LrDpleSolver::doc(const std::string& name) {
    return LrDpleInternal::getPlugin(name).doc;
  }

  LrDpleSolver::LrDpleSolver(const std::string& name,
                         const LrDpleStructure & st,
                         const std::vector< std::vector<int> > &Hs) {
    assignNode(LrDpleInternal::getPlugin(name).creator(st, Hs));
  }

} // namespace casadi

