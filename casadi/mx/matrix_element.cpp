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
using namespace std;

namespace CasADi{

MatrixElement::MatrixElement(const MX& x, const std::vector<int>& ii, const std::vector<int>& jj) : ii_(ii), jj_(jj){
  setDependencies(x);
  setSize(ii.size(),jj.size());
}

MatrixElement* MatrixElement::clone() const{
  return new MatrixElement(*this);
}

void MatrixElement::print(std::ostream &stream) const{
  if(dep(0).size2()==1)
    if(ii_.size()==1)
      stream << dep(0) << "[" << ii_[0] << "]";
    else
      stream << dep(0) << "[" << ii_ << "]";
  else
    if(ii_.size()==1 && jj_.size()==1)
      stream << dep(0) << "(" << ii_[0] << "," << jj_[0] << ")";
    else
      stream << dep(0) << "(" << ii_ << "," << jj_ << ")";
}

void MatrixElement::evaluate(int fsens_order, int asens_order){
  if(ii_.size() != 1 || jj_.size() != 1)
    throw CasadiException("MatrixElement::evaluate: only scalar submatrices implemented");
  
  int i = ii_[0], j = jj_[0];
  
  // Get references to the terms
  const vector<double>& arg = input(0); // first term
  vector<double>& res = output();
  
  // carry out the assignment
  res[0] = arg[j+i*dep(0).size2()];
  
  if(fsens_order>0){
    // Get references to the terms
    const vector<double>& fseed = fwdSeed(0); // first term
    vector<double>& fsens = fwdSens(0);
  
    // carry out the assignment
    fsens[0] = fseed[j+i*dep(0).size2()];
  }
  
  if(asens_order>0){
    // Get references to the terms
    const vector<double>& aseed = adjSeed();
    vector<double>& asens = adjSens(0); // first term
  
    // carry out the addition
    asens[j+i*dep(0).size2()] += aseed[0];
  }
}

} // namespace CasADi
