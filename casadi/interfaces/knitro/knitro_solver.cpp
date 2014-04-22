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

#include "knitro_internal.hpp"

using namespace std;

namespace casadi {

  KnitroSolver::KnitroSolver() {
  }

  KnitroSolver::KnitroSolver(const Function& nlp) {
    assignNode(new KnitroInternal(nlp));
  }

  KnitroInternal* KnitroSolver::operator->() {
    return static_cast<KnitroInternal*>(NLPSolver::operator->());
  }

  const KnitroInternal* KnitroSolver::operator->() const {
    return static_cast<const KnitroInternal*>(NLPSolver::operator->());
  }

  bool KnitroSolver::checkNode() const {
    return dynamic_cast<const KnitroInternal*>(get());
  }

  void KnitroSolver::setIntParam(const std::string& name, int val) {
    (*this)->int_param_[name] = val;
  }

  void KnitroSolver::setDoubleParam(const std::string& name, double val) {
    (*this)->double_param_[name] = val;
  }

  void KnitroSolver::setStringParam(const std::string& name, const std::string& val) {
    (*this)->string_param_[name] = val;
  }

} // namespace casadi
