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

#include "direct_collocation_internal.hpp"

namespace casadi {

DirectCollocation::DirectCollocation() {
}

DirectCollocation::DirectCollocation(const Function& ffcn, const Function& mfcn,
                                     const Function& cfcn, const Function& rfcn) {
  assignNode(new DirectCollocationInternal(ffcn, mfcn, cfcn, rfcn));
}

const DirectCollocationInternal* DirectCollocation::operator->() const {
  return static_cast<const DirectCollocationInternal*>(Function::operator->());
}

DirectCollocationInternal* DirectCollocation::operator->() {
  return static_cast<DirectCollocationInternal*>(Function::operator->());
}

void DirectCollocation::getGuess(std::vector<double>& V_init) const {
  (*this)->getGuess(V_init);
}

void DirectCollocation::getVariableBounds(std::vector<double>& V_min,
                                          std::vector<double>& V_max) const {
  (*this)->getVariableBounds(V_min, V_max);
}

void DirectCollocation::getConstraintBounds(std::vector<double>& G_min,
                                            std::vector<double>& G_max) const {
  (*this)->getConstraintBounds(G_min, G_max);
}

void DirectCollocation::setOptimalSolution(const std::vector<double> &V_opt) {
  (*this)->setOptimalSolution(V_opt);
}

  NlpSolver DirectCollocation::getNlpSolver() const
  { return isNull() ? NlpSolver(): (*this)->nlp_solver_; }

void DirectCollocation::reportConstraints(std::ostream &stream) {
  (*this)->reportConstraints();
}

} // namespace casadi

