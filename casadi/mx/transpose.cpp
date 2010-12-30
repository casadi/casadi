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
  sz.nrow = x.size2();
  sz.ncol = x.size1();
}

Transpose* Transpose::clone() const{
  return new Transpose(*this);
}

int Transpose::k2k(int k) {
  return (k/sz.ncol) + (k%sz.ncol)*sz.nrow;
}

void Transpose::print(std::ostream &stream) const{
  stream << "trans(" << dep(0) << ")";
}

void Transpose::evaluate(int fsens_order, int asens_order){
  if(fsens_order==0 || asens_order==0);
  
  if(fsens_order==0){
  // Get references to the terms
  const vector<double>& arg = dep(0)->val(0);
  vector<double>& res = val(0);
  
  // carry out the transpose
  for(int i=0; i<sz.nrow; ++i)
    for(int j=0; j<sz.ncol; ++j)
      res[j + i*sz.ncol] = arg[i + j*sz.nrow];
  } else {

    // Get references to the terms
    const vector<double>& arg = dep(0)->val(1);
    vector<double>& res = val(1);
  
    // carry out the transpose
    for(int i=0; i<sz.nrow; ++i)
      for(int j=0; j<sz.ncol; ++j)
        res[j + i*sz.ncol] = arg[i + j*sz.nrow];
  }
  
  if(asens_order>0){
    // Get references to the terms
    vector<double>& arg = dep(0)->val(1);
    const vector<double>& res = val(1);
  
    // carry out the transpose
    for(int i=0; i<sz.nrow; ++i)
      for(int j=0; j<sz.ncol; ++j)
        arg[i + j*sz.nrow] += res[j + i*sz.ncol];
    
  }
}


} // namespace CasADi
