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

#include "matrix_element.hpp"
#include <cassert>
using namespace std;

namespace CasADi{

MatrixElement::MatrixElement(const MX& x, int i_, int j_) : MXNode(x), i(i_), j(j_){
  assert(i>=0 && i<x.size1());
  assert(j>=0 && j<x.size2());
  nrow_ = 1;
  ncol_ = 1;
}

MatrixElement* MatrixElement::clone() const{
  return new MatrixElement(*this);
}

void MatrixElement::print(std::ostream &stream) const{
  if(dep(0).size2()==1)
    stream << dep(0) << "[" << i << "]";
  else
    stream << dep(0) << "(" << i << "," << j << ")";
}

void MatrixElement::evaluate(int fsens_order, int asens_order){
  assert(fsens_order==0 || asens_order==0);
  
  if(fsens_order==0){
  // Get references to the terms
  const vector<double>& arg = dep(0)->val(0); // first term
  vector<double>& res = val(0);
  
  // carry out the assignment
  res[0] = arg[j+i*dep(0).size2()];
  } else {
    // Get references to the terms
    const vector<double>& arg = dep(0)->val(1); // first term
    vector<double>& res = val(1);
  
    // carry out the assignment
    res[0] = arg[j+i*dep(0).size2()];
  }
  
  if(asens_order>0){
    // Get references to the terms
    vector<double>& arg = dep(0)->val(1); // first term
    const vector<double>& res = val(1);
  
    // carry out the addition
    // where's the plus sign?
    arg[j+i*dep(0).size2()] = res[0];
    
  }
}

} // namespace CasADi
