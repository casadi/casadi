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

#include "transpose.hpp"
#include "../stl_vector_tools.hpp"

using namespace std;

namespace CasADi{


Transpose::Transpose(const MX& x){
  setDependencies(x);
  setSize(x.size2(),x.size1());
}

Transpose* Transpose::clone() const{
  return new Transpose(*this);
}

void Transpose::init(){
  // Base class initializations
  Reordering::init();

  // Mapping of the non-zeros: move to constructor?
  CRSSparsity sp = dep(0)->output().sparsity().transpose(nzind_);
  if(sp.size()!=sp.numel())
    throw CasadiException("Error in Transpose::init: sparse nodes not yet allowed");
  
  // Where to store this information?
}

void Transpose::print(std::ostream &stream) const{
  stream << "trans(" << dep(0) << ")";
}


} // namespace CasADi
