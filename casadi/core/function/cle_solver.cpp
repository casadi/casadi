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


#include "cle_solver.hpp"
#include "cle_internal.hpp"

using namespace std;
namespace casadi {

  CleSolver::CleSolver() {

  }

  CleSolver  CleSolver::clone() const {
    CleSolver ret;
    if (!isNull()) ret.assignNode((*this)->clone());
    return ret;
  }

  void CleSolver::printStats(ostream &stream) const {
    (*this)->printStats(stream);
  }

  CleInternal* CleSolver::operator->() {
    return static_cast<CleInternal*>(Function::operator->());
  }

  const CleInternal* CleSolver::operator->() const {
    return static_cast<const CleInternal*>(Function::operator->());
  }

  bool CleSolver::testCast(const SharedObjectNode* ptr) {
    return dynamic_cast<const CleInternal*>(ptr)!=0;
  }

  bool CleSolver::hasPlugin(const std::string& name) {
    return CleInternal::hasPlugin(name);
  }

  void CleSolver::loadPlugin(const std::string& name) {
    CleInternal::loadPlugin(name);
  }

  std::string CleSolver::doc(const std::string& name) {
    return CleInternal::getPlugin(name).doc;
  }

  CleSolver::CleSolver(const std::string& name,
                         const CleStructure& st) {
    assignNode(CleInternal::instantiatePlugin(name, st));
  }

} // namespace casadi

