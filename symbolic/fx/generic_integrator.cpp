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

#include "generic_integrator.hpp"
#include "generic_integrator_internal.hpp"
#include <cassert>

using namespace std;
namespace CasADi{

  GenericIntegrator::GenericIntegrator(){
  }
  
  GenericIntegratorInternal* GenericIntegrator::operator->(){
    return static_cast<GenericIntegratorInternal*>(FX::operator->());
  }

  const GenericIntegratorInternal* GenericIntegrator::operator->() const{
    return static_cast<const GenericIntegratorInternal*>(FX::operator->()); 
  }
      
  bool GenericIntegrator::checkNode() const{
    return dynamic_cast<const GenericIntegratorInternal*>(get())!=0;
  }
 
} // namespace CasADi

