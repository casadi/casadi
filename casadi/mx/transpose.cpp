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
#include <cassert>

using namespace std;

namespace CasADi{


Transpose::Transpose(const MX& x) : Reordering(x){
  setSize(x.size2(),x.size1());
}

Transpose* Transpose::clone() const{
  return new Transpose(*this);
}

int Transpose::k2k(int k) {
  return (k/size2()) + (k%size2())*size1();
}

void Transpose::print(std::ostream &stream) const{
  stream << "trans(" << dep(0) << ")";
}

void Transpose::evaluate(int fsens_order, int asens_order){
  if(fsens_order==0 || asens_order==0);
  
  if(fsens_order==0){
  // Get references to the terms
  const vector<double>& arg = input(0);
  vector<double>& res = output();
  
  // carry out the transpose
  for(int i=0; i<size1(); ++i)
    for(int j=0; j<size2(); ++j)
      res[j + i*size2()] = arg[i + j*size1()];
  } else {

    // Get references to the terms
    const vector<double>& arg = fwdSeed(0);
    vector<double>& res = fwdSens();
  
    // carry out the transpose
    for(int i=0; i<size1(); ++i)
      for(int j=0; j<size2(); ++j)
        res[j + i*size2()] = arg[i + j*size1()];
  }
  
  if(asens_order>0){
    // Get references to the terms
    vector<double>& arg = adjSens(0);
    const vector<double>& res = adjSeed();
  
    // carry out the transpose
    for(int i=0; i<size1(); ++i)
      for(int j=0; j<size2(); ++j)
        arg[i + j*size1()] += res[j + i*size2()];
    
  }
}


} // namespace CasADi
