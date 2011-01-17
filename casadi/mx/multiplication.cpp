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

#include "multiplication.hpp"
#include <vector>
#include <cassert>

using namespace std;

namespace CasADi{

Multiplication::Multiplication(const MX& x, const MX& y){
  setDependencies(x,y);
  if(x.size2() != y.size1()) throw CasadiException("Multiplication::dimension mismatch");
  setSize(x.size1(),y.size2());
  ni = dep(0).size1();
  nk = dep(0).size2();
  nj = dep(1).size2();
}

Multiplication* Multiplication::clone() const{
  return new Multiplication(*this);
}

void Multiplication::print(std::ostream &stream) const{
  stream << "prod(" << dep(0) << "," << dep(1) << ")";
}

void Multiplication::matrix_matrix_mult(const vector<double>& t1, const vector<double>& t2, vector<double>& t3){
  // NOTE: Remove as it does not exploit sparsity - use a method in Matrix<> instead!
  // t3 = t1*t2
  for(int i=0; i<ni; ++i)
    for(int j=0; j<nj; ++j)
      for(int k=0; k<nk; ++k)
        t3[i + ni*j] += t1[i + ni*k]*t2[k + nk*j];
}

void Multiplication::matrix_matrix_mult1(vector<double>& t1, const vector<double>& t2, const vector<double>& t3){
  // NOTE: Remove as it does not exploit sparsity - use a method in Matrix<> instead!
  // t1 = t3*trans(t2)
  for(int i=0; i<ni; ++i)
    for(int k=0; k<nk; ++k)
      for(int j=0; j<nj; ++j)
        t1[i + ni*k] += t3[i + ni*j]*t2[k + nk*j];
}

void Multiplication::matrix_matrix_mult2(const vector<double>& t1, vector<double>& t2, const vector<double>& t3){
  // NOTE: Remove as it does not exploit sparsity - use a method in Matrix<> instead!
  // t2 = trans(t1)*t3
  for(int k=0; k<nk; ++k)
    for(int j=0; j<nj; ++j)
      for(int i=0; i<ni; ++i)
        t2[k + nk*j] += t1[i + ni*k]*t3[i + ni*j];
}


void Multiplication::evaluate(int fsens_order, int asens_order){
  assert(fsens_order==0 || asens_order==0);
  
  if(fsens_order==0){
  // Erase result
  for(vector<double>::iterator it=output().begin(); it!=output().end(); ++it) *it = 0; 

  // Matrix multiplication
  matrix_matrix_mult(input(0),input(1),output());
} else {

  //   // Erase result
//   for(vector<double>::iterator it=res.begin(); it!=res.end(); ++it) *it = 0; 

  // Matrix multiplication, first argument
  matrix_matrix_mult(fwdSeed(0),input(1),fwdSens());

  // Matrix multiplication, second argument
  matrix_matrix_mult(input(0),fwdSeed(1),fwdSens());
}

if(asens_order>0){
  // Matrix multiplication, first argument
  matrix_matrix_mult1(adjSens(0),input(1),adjSeed());

  // Matrix multiplication, second argument
  matrix_matrix_mult2(input(0),adjSens(1),adjSeed());
}
}


} // namespace CasADi

