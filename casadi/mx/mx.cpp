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

#include "mx.hpp"
#include "mx_node.hpp"
#include "scalar_matrix_op.hpp"
#include "matrix_scalar_op.hpp"
#include "matrix_matrix_op.hpp"
#include "unary_op.hpp"
#include "matrix_element.hpp"
#include "../fx/sx_function.hpp"
#include "evaluation.hpp"
#include "symbolic_mx_node.hpp"
#include "mx_constant.hpp"


namespace CasADi{

MX::~MX(){
}

MX::MX(){
}

MX::MX(double x){
  assignNode(new MXConstant(x));
}

MX::MX(const Matrix<double> &x){
  assignNode(new MXConstant(x));
}

MX::MX(const std::vector<double> &x){
  assignNode(new MXConstant(x));
}

MX::MX(const std::string& name, int n, int m){
  assignNode(new SymbolicMatrix(name,n,m));
}

const MX MX::__getitem__(const vector<int>& I) const{
  if(I.size()!=2) throw CasADi::CasadiException("__getitem__: not 2D"); 
  return (*this)(I[0],I[1]);
}

const MX MX::__getitem__(int k) const{
  // change this!
  return (*this)(k/size2(),k%size2());
}

MX& MX::__setitem__(int k, const MX& el){ 
  (*this)[k] = el;
  return *this;
}

MX& MX::__setitem__(const std::vector<int> &I, const MX&  el){ 
  if(I.size()!=2) throw CasADi::CasadiException("__setitem__: not 2D"); 
  (*this)(I[0],I[1]) = el;
  return *this;
}




const MX MX::getSub(const std::vector<int>& ii, const std::vector<int>& jj) const{
  casadi_assert(ii.size()==1 && jj.size()==1);
  
  MX ret;
  ret.assignNode(new MatrixElement(*this,ii[0],jj[0]));
  return ret;
}

void MX::setSub(const std::vector<int>& ii, const std::vector<int>& jj, const MX& el){
  throw CasadiException("MX::setSub: not implemented");
}

const MX MX::operator()(int i, int j) const{
  MX ret;
  ret.assignNode(new MatrixElement(*this,i, j));
  return ret;  
}

const MX MX::operator[](int k) const{
  return __getitem__(k);
}

SubMatrix<MX> MX::operator()(int i, int j){
  return SubMatrix<MX>(*this,vector<int>(1,i), vector<int>(1,j));
}
 
SubMatrix<MX> MX::operator[](int k){
  // change this!!!
  return SubMatrix<MX>(*this,vector<int>(1,k/size2()), vector<int>(1,k%size2()));
}

int MX::size() const{
  return sparsity().size();
}

int MX::size1() const{
  return sparsity().size1();
}

int MX::numel() const{
  return sparsity().numel();
}

int MX::size2() const{
  return sparsity().size2();
}

MX operator+(const MX &x, const MX &y){
  return MX::binary(ADD_NODE,x,y);
}

MX operator-(const MX &x, const MX &y){
  return MX::binary(SUB_NODE,x,y);
}

MX MX::binary(int op, const MX &x, const MX &y){
  if(x.numel()==1)
    return scalar_matrix(op,x,y);
  else if(y.numel()==1)  
    return matrix_scalar(op,x,y);
  else
    return matrix_matrix(op,x,y);
}

MX MX::unary(int op, const MX &x){
  MX ret;
  ret.assignNode(new UnaryOp(OPERATION(op),x));
  return ret;
}

MX MX::scalar_matrix(int op, const MX &x, const MX &y){
  MX ret;
  ret.assignNode(new ScalarMatrixOp(OPERATION(op),x,y));  
  return ret;
}

MX MX::matrix_scalar(int op, const MX &x, const MX &y){
  MX ret;
  ret.assignNode(new MatrixScalarOp(OPERATION(op),x,y));  
  return ret;
}

MX MX::matrix_matrix(int op, const MX &x, const MX &y){
  MX ret;
  ret.assignNode(new MatrixMatrixOp(OPERATION(op),x,y)); 
  return ret;
}

MX operator*(const MX &x, const MX &y){
  return MX::binary(MUL_NODE,x,y);
}

MX operator/(const MX &x, const MX &y){
  return MX::binary(DIV_NODE,x,y);
}

MXNode* MX::operator->(){
  return (MXNode*)SharedObject::operator->();
}

const MXNode* MX::operator->() const{
  return (const MXNode*)SharedObject::operator->();
}

MX& MX::operator+=(const MX &y){
  return *this = *this + y;
}

MX& MX::operator-=(const MX &y){
  return *this = *this - y;
}

MX& MX::operator*=(const MX &y){
  return *this = *this * y;
}

MX& MX::operator/=(const MX &y){
  return *this = *this / y;
}

bool MX::empty() const{
  return numel()==0;
}

MX MX::zeros(int nrow, int ncol){
  return MX(Matrix<double>(nrow,ncol,0));
}

MX MX::ones(int nrow, int ncol){
  return MX(Matrix<double>(nrow,ncol,1));
}

MX MX::operator-() const{
  return unary(NEG_NODE,*this);
}

const CRSSparsity& MX::sparsity() const{
  return (*this)->sparsity();
}

} // namespace CasADi


using namespace CasADi;
namespace std{

MX exp(const MX& x){
  return MX::unary(EXP_NODE,x);
}

MX log(const MX& x){
  return MX::unary(LOG_NODE,x);
}

MX sqrt(const MX& x){
  return MX::unary(SQRT_NODE,x);
}

MX sin(const MX& x){
  return MX::unary(SIN_NODE,x);
}

MX cos(const MX& x){;  
  return MX::unary(COS_NODE,x);
}

MX tan(const MX& x){
  return MX::unary(TAN_NODE,x);
}

MX atan(const MX& x){
  return MX::unary(ATAN_NODE,x);
}

MX asin(const MX& x){
  return MX::unary(ASIN_NODE,x);
}

MX acos(const MX& x){
  return MX::unary(ACOS_NODE,x);
}

MX pow(const MX& x, const MX& n){
  return MX::binary(POW_NODE,x,n);
}

MX erf(const MX& x){
  return MX::unary(ERF_NODE,x);
}

MX floor(const MX& x){
  return MX::unary(FLOOR_NODE,x);
}

MX ceil(const MX& x){
  return MX::unary(CEIL_NODE,x);
}

MX fmin(const MX& x, const MX& y){
  return MX::binary(FMIN_NODE,x,y);
}

MX fmax(const MX& x, const MX& y){
  return MX::binary(FMAX_NODE,x,y);
}



} // namespace std
