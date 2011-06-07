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

#include "c_function_internal.hpp"

#include <iostream>
#include <fstream>
#include <sstream>

namespace CasADi{

using namespace std;

CFunction::CFunction(){
}

CFunction::CFunction(CFunctionWrapper c_fcn){
  assignNode(new CFunctionInternal(c_fcn));
}

CFunctionInternal* CFunction::operator->(){
  return (CFunctionInternal*)(FX::operator->());
}

const CFunctionInternal* CFunction::operator->() const{
   return (const CFunctionInternal*)(FX::operator->()); 
}
  
bool CFunction::checkNode() const{
  return dynamic_cast<const CFunctionInternal*>(get())!=0;
}

void CFunction::setUserData(void* user_data){
  (*this)->setUserData(user_data);
}


} // namespace CasADi

