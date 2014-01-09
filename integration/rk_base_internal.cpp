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

#include "rk_base_internal.hpp"
#include "symbolic/stl_vector_tools.hpp"
#include "symbolic/matrix/sparsity_tools.hpp"
#include "symbolic/matrix/matrix_tools.hpp"
#include "symbolic/sx/sx_tools.hpp"
#include "symbolic/fx/sx_function.hpp"
#include "symbolic/mx/mx_tools.hpp"

using namespace std;
namespace CasADi{

  RKBaseInternal::RKBaseInternal(const FX& f, const FX& g) : IntegratorInternal(f,g){
    addOption("number_of_finite_elements",     OT_INTEGER,  20, "Number of finite elements");
  }

  void RKBaseInternal::deepCopyMembers(std::map<SharedObjectNode*,SharedObject>& already_copied){    
    IntegratorInternal::deepCopyMembers(already_copied);
  }

  RKBaseInternal::~RKBaseInternal(){
  }

  void RKBaseInternal::init(){
    // Call the base class init
    IntegratorInternal::init();
  
    // Number of finite elements
    int nk = getOption("number_of_finite_elements");
  
    casadi_error("Not implemented");
  }
  
  void RKBaseInternal::reset(){
    // Call the base class method
    IntegratorInternal::reset();
  
    casadi_error("Not implemented");
  }

  void RKBaseInternal::resetB(){
    casadi_error("Not implemented");
  }

  void RKBaseInternal::integrate(double t_out){
    casadi_error("Not implemented");
  }

  void RKBaseInternal::integrateB(double t_out){
    casadi_error("Not implemented");
  }

} // namespace CasADi
