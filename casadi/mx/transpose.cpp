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
  
  // Mapping: NOTE: does not handle sparsity correctly, should get the reordering from the Matrix<> class (which will calculate it using binsorting)
  for(int k=0; k<size1()*size2(); ++k)
    nzind_[k] = (k/size2()) + (k%size2())*size1();
}

void Transpose::print(std::ostream &stream) const{
  stream << "trans(" << dep(0) << ")";
}


} // namespace CasADi
