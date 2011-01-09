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

#include "flatten.hpp"
#include <algorithm>

using namespace std;

namespace CasADi{

Flatten::Flatten(const MX& x){
  setDependencies(x);
  setSize(x.numel(),1);
}

Flatten* Flatten::clone() const{
  return new Flatten(*this);
}

void Flatten::print(std::ostream &stream) const{
  stream << "flatten(" << dep(0) << ")";
}

void Flatten::evaluate(int fsens_order, int asens_order){
  // All non-zero elements remains the same
  if(fsens_order==0){
    copy(input(0).begin(),input(0).end(),output().begin());
  } else {
    copy(fwdSeed(0).begin(),fwdSeed(0).end(),fwdSens().begin());
  }
  
  if(asens_order>0){
    transform(adjSeed().begin(),adjSeed().end(),
              adjSens(0).begin(),
              adjSens(0).begin(),
              plus<double>());
  }
}

} // namespace CasADi
