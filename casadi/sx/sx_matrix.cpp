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

SXMatrix::SXMatrix(const SXMatrix::Element &el) : Matrix<SX>(el.mat.getElement(el.i,el.j)){
}

// destructor
SXMatrix::~SXMatrix(){
}

SX& SXMatrix::operator[](int k){
  return at(k); //  return Element(*this,k%size1(), k/size1());
}

const SX& SXMatrix::operator[](int k) const{
  return at(k); //  return getElement(k%size1(), k/size1());
}

const std::vector<SX> SXMatrix::operator[](const std::vector<int>& ind) const{
  std::vector<SX> ret;
  for(std::vector<int>::const_iterator it=ind.begin(); it!=ind.end(); ++it)
    ret.push_back(at(*it));
  return ret;
}

const SX SXMatrix::operator()(int i, int j) const{
  return getElement(i,j);
}


SXMatrix::Element SXMatrix::operator()(int i, int j){
  return Element(*this,i,j);
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

SXMatrix::Element::Element(SXMatrix& mat_, int i_, int j_) : mat(mat_), i(i_), j(j_){ 
}

SXNode* const SXMatrix::Element::get() const{
  return mat.getElement(i,j).get();
}

const SXNode* SXMatrix::Element::operator->() const{
  return mat.getElement(i,j).get();  
}

SXNode* SXMatrix::Element::operator->(){
  return mat.getElement(i,j).get();  
}

SX& SXMatrix::Element::operator=(const SXMatrix::Element &y){
  return mat.getElementRef(i,j) = y.mat.getElement(y.i,y.j);
}

SX& SXMatrix::Element::operator=(const SXMatrix &y){
  if(!y.scalar()) throw CasadiException("operator=: argument not scalar");
  return mat.getElementRef(i,j) = y(0);
}


ostream& operator<<(ostream &stream, const SXMatrix &mat){
  mat.print(stream);
  return stream;
}


} //namespace CasADi

namespace std{
using namespace CasADi;


SXMatrix sin(const SXMatrix& x){
  return unary(SIN_NODE,x);
}

SXMatrix cos(const SXMatrix& x){
  return unary(COS_NODE,x);
}

SXMatrix tan(const SXMatrix& x){
  return unary(TAN_NODE,x);
}

SXMatrix atan(const SXMatrix& x){
  return unary(ATAN_NODE,x);
}

SXMatrix asin(const SXMatrix& x){
  return unary(ASIN_NODE,x);
}

SXMatrix acos(const SXMatrix& x){
  return unary(ACOS_NODE,x);
}

SXMatrix exp(const SXMatrix& x){
  return unary(EXP_NODE,x);
}

SXMatrix log(const SXMatrix& x){
  return unary(LOG_NODE,x);
}

SXMatrix pow(const SXMatrix& x, const SXMatrix& n){
  return binary(POW_NODE,x,n);
}

SXMatrix sqrt(const SXMatrix& x){
  return unary(SQRT_NODE,x);
}

SXMatrix fmin(const SXMatrix& x, const SXMatrix& y){
  return binary(FMIN_NODE,x,y);
}

SXMatrix fmax(const SXMatrix& x, const SXMatrix& y){
  return binary(FMAX_NODE,x,y);
}

SXMatrix floor(const SXMatrix& x){
  return unary(FLOOR_NODE,x);
}

SXMatrix ceil(const SXMatrix& x){
  return unary(CEIL_NODE,x);
}

SXMatrix erf(const SXMatrix& x){
  return unary(ERF_NODE,x);
}


} // namespace std

