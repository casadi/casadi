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

#include "ocp_internal.hpp"

using namespace std;

namespace CasADi{

OCP::OCP(){
}

OCP::OCP(const std::vector<FX>& L, const std::vector<FX>& F, const std::vector<FX>& H, const std::vector<FX>& G){
  assignNode(new OCPInternal(L,F,H,G));
}
  
const OCPInternal* OCP::operator->() const{
  return (const OCPInternal*)FX::operator->();
}

OCPInternal* OCP::operator->(){
  return (OCPInternal*)FX::operator->();
}

} // namespace CasADi

