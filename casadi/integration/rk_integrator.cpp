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

#include "rk_integrator_internal.hpp"

using namespace std;

namespace casadi {

  RKIntegrator::RKIntegrator() {
  }

  RKIntegrator::RKIntegrator(const Function& f, const Function& g) {
    assignNode(new RKIntegratorInternal(f, g));
  }

  RKIntegratorInternal* RKIntegrator::operator->() {
    return static_cast<RKIntegratorInternal*>(Integrator::operator->());
  }

  const RKIntegratorInternal* RKIntegrator::operator->() const {
    return static_cast<const RKIntegratorInternal*>(Integrator::operator->());
  }

  bool RKIntegrator::checkNode() const {
    return dynamic_cast<const RKIntegratorInternal*>(get())!=0;
  }

} // namespace casadi
