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

#include "generic_integrator_internal.hpp"
#include "../stl_vector_tools.hpp"
#include "../matrix/matrix_tools.hpp"
#include "../mx/mx_tools.hpp"
#include "../sx/sx_tools.hpp"
#include "mx_function.hpp"
#include "sx_function.hpp"

//INPUTSCHEME(IntegratorInput)
//OUTPUTSCHEME(IntegratorOutput)

using namespace std;
namespace CasADi{

  GenericIntegratorInternal::GenericIntegratorInternal(const FX& f, const FX& g){
    //inputScheme_ = SCHEME_IntegratorInput;
    //outputScheme_ = SCHEME_IntegratorOutput;
  }

  GenericIntegratorInternal::~GenericIntegratorInternal(){ 
  }

  void GenericIntegratorInternal::init(){
    
    // Call the base class method
    FXInternal::init();
  }

  void GenericIntegratorInternal::deepCopyMembers(std::map<SharedObjectNode*,SharedObject>& already_copied){
    FXInternal::deepCopyMembers(already_copied);
  }


} // namespace CasADi


