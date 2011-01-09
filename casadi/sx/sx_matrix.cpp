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

#include "sx_matrix.hpp"
#include "sx_tools.hpp"

using namespace std;

namespace CasADi{


SXMatrix::SXMatrix(const std::vector<SX>& x) : Matrix<SX>(x){
}

SXMatrix::SXMatrix(const std::vector<double>& x) : Matrix<SX>(x){
}

SXMatrix::SXMatrix(const std::vector<SX>& x,  int n, int m) : Matrix<SX>(x,n,m){
}

SXMatrix::SXMatrix(const std::vector<double>& x,  int n, int m) : Matrix<SX>(x,n,m){
}

SXMatrix::SXMatrix() : Matrix<SX>(){
}

SXMatrix::SXMatrix(double val) : Matrix<SX>(val){
}

SXMatrix::SXMatrix(int n, int m) : Matrix<SX>(n,m){
}

SXMatrix::SXMatrix(int n, int m, const SX& val) : Matrix<SX>(n,m,val){
}

SXMatrix::SXMatrix(const string& name, int n, int m) : Matrix<SX>(create_symbolic(name,n*m),n,m){
}

SXMatrix::SXMatrix(const SX& scalar) : Matrix<SX>(scalar){
}

SXMatrix& operator+=(SXMatrix &ex, const SXMatrix &expr){
  return ex = ex + expr;
}

SXMatrix operator+(const SXMatrix &x, const SXMatrix &y){
  return binary(ADD_NODE,x,y);
}

SXMatrix& operator-=(SXMatrix &ex, const SXMatrix &expr){
 return ex = ex - expr;
}

SXMatrix operator-(SXMatrix &ex){
  return unary(NEG_NODE,ex);
}

SXMatrix operator-(const SXMatrix &x, const SXMatrix &y){ 
  return binary(SUB_NODE,x,y);
}

SXMatrix& operator*=(SXMatrix &ex, const SXMatrix &expr){
 return ex = ex * expr;
}

SXMatrix operator*(const SXMatrix &x, const SXMatrix &y){
  return binary(MUL_NODE,x,y);
}


SXMatrix& operator/=(SXMatrix &ex, const SXMatrix &expr){
  return ex = ex / expr;
}

SXMatrix operator/(const SXMatrix &x, const SXMatrix &y){
  return binary(DIV_NODE,x,y);
}

ostream& operator<<(ostream &stream, const SXMatrix &mat){
  mat.print(stream);
  return stream;
}


} //namespace CasADi

