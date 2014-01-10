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
namespace CasADi{

  RKIntegratorInternal::RKIntegratorInternal(const FX& f, const FX& g) : RKBaseInternal(f,g){
    addOption("interpolation_order",           OT_INTEGER,  4,  "Order of the interpolating polynomials");
  }

  void RKIntegratorInternal::deepCopyMembers(std::map<SharedObjectNode*,SharedObject>& already_copied){
    RKBaseInternal::deepCopyMembers(already_copied);
  }

  RKIntegratorInternal::~RKIntegratorInternal(){
  }

  void RKIntegratorInternal::init(){
    // Call the base class init
    RKBaseInternal::init();

    casadi_error("Not implemented");
  }
  
} // namespace CasADi
