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


#include "homotopy_nlp_solver.hpp"
#include "homotopy_nlp_internal.hpp"
#include "../sx/sx_tools.hpp"
#include "../mx/mx_tools.hpp"

using namespace std;
namespace casadi {

  HomotopyNlpSolver::HomotopyNlpSolver() {
  }

  HomotopyNLPInternal* HomotopyNlpSolver::operator->() {
    return static_cast<HomotopyNLPInternal*>(Function::operator->());
  }

  const HomotopyNLPInternal* HomotopyNlpSolver::operator->() const {
    return static_cast<const HomotopyNLPInternal*>(Function::operator->());
  }

  bool HomotopyNlpSolver::testCast(const SharedObjectNode* ptr) {
    return dynamic_cast<const HomotopyNLPInternal*>(ptr)!=0;
  }

  bool HomotopyNlpSolver::hasPlugin(const std::string& name) {
    return HomotopyNLPInternal::hasPlugin(name);
  }

  void HomotopyNlpSolver::loadPlugin(const std::string& name) {
    HomotopyNLPInternal::loadPlugin(name);
  }

  std::string HomotopyNlpSolver::doc(const std::string& name) {
    return HomotopyNLPInternal::getPlugin(name).doc;
  }

} // namespace casadi

