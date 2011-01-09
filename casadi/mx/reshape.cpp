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

#include "reshape.hpp"
#include <algorithm>

using namespace std;

namespace CasADi{

Reshape::Reshape(const MX& x, int n, int m){
  if (n*m != x.numel()) {
    throw CasadiException("MX::reshape: size must be same before and after reshaping");
  }
  setDependencies(x);
  setSize(n,m);
}

Reshape* Reshape::clone() const{
  return new Reshape(*this);
}

void Reshape::print(std::ostream &stream) const{
  stream << "reshape(" << dep(0) << ",[" << size1() << "," << size2() << "])";
}

void Reshape::evaluate(int fsens_order, int asens_order){
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
