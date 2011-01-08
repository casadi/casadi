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
  nrow_ = x.size2();
  ncol_ = x.size1();
}

Transpose* Transpose::clone() const{
  return new Transpose(*this);
}

int Transpose::k2k(int k) {
  return (k/ncol_) + (k%ncol_)*nrow_;
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
  for(int i=0; i<nrow_; ++i)
    for(int j=0; j<ncol_; ++j)
      res[j + i*ncol_] = arg[i + j*nrow_];
  } else {

    // Get references to the terms
    const vector<double>& arg = dep(0)->val(1);
    vector<double>& res = val(1);
  
    // carry out the transpose
    for(int i=0; i<nrow_; ++i)
      for(int j=0; j<ncol_; ++j)
        res[j + i*ncol_] = arg[i + j*nrow_];
  }
  
  if(asens_order>0){
    // Get references to the terms
    vector<double>& arg = dep(0)->val(1);
    const vector<double>& res = val(1);
  
    // carry out the transpose
    for(int i=0; i<nrow_; ++i)
      for(int j=0; j<ncol_; ++j)
        arg[i + j*nrow_] += res[j + i*ncol_];
    
  }
}


} // namespace CasADi
