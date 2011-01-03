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

#include "acado_internal.hpp"

#include <acado_optimal_control.hpp>

#include <cassert>
#include <limits>

using namespace std;

namespace CasADi{

AcadoInterface::AcadoInterface(){ 
}

AcadoInterface::AcadoInterface(const FX& ffcn, const FX& mfcn, const FX& cfcn, const FX& rfcn){
  assignNode(new AcadoInternal(ffcn, mfcn, cfcn, rfcn));
}

void AcadoInterface::setIntegrators(const vector<Integrator>& integrators){
  (*this)->setIntegrators(integrators);
}

AcadoInternal* AcadoInterface::operator->(){
  return (AcadoInternal*)(FX::operator->());
}

const AcadoInternal* AcadoInterface::operator->() const{
   return (const AcadoInternal*)(FX::operator->()); 
}

bool AcadoInterface::checkNode() const{
  return dynamic_cast<const AcadoInternal*>(get());
}


} // namespace CasADi

