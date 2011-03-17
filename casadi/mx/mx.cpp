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
#include "mapping.hpp"
#include "../fx/sx_function.hpp"
#include "evaluation.hpp"
#include "symbolic_mx_node.hpp"
#include "mx_constant.hpp"
#include "mx_zero.hpp"
#include "mx_tools.hpp"
#include "../stl_vector_tools.hpp"
#include "../matrix/matrix_tools.hpp"

using namespace std;
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

MX::MX(int nrow, int ncol){
  assignNode(new MXZero(nrow,ncol));
}

const MX MX::getitem(const vector<int>& I) const{
  if(I.size()!=2) throw CasADi::CasadiException("getitem: not 2D"); 
  return (*this)(I[0],I[1]);
}

const MX MX::getitem(int k) const{
  // change this!
  return (*this)[k];
}

MX& MX::setitem(int k, const MX& el){ 
  (*this)[k] = el;
  return *this;
}

MX& MX::setitem(const std::vector<int> &I, const MX&  el){ 
  if(I.size()!=2) throw CasADi::CasadiException("setitem: not 2D"); 
  (*this)(I[0],I[1]) = el;
  return *this;
}

const MX MX::getitem(const std::vector< std::vector<int> > &II) const{
  casadi_assert_message(II.size()==2,"Index vector must be two-dimensional");
  return (*this)(II[0],II[1]);
}

MX& MX::setitem(const std::vector< std::vector<int> > &II, const MX&  m){
  casadi_assert_message(II.size()==2,"Index vector must be two-dimensional");
  setSub(II[0],II[1],m);
}

MX MX::getSub(const std::vector<int>& ii, const std::vector<int>& jj) const{
  // Nonzero mapping from submatrix to full
  vector<int> mapping;
  
  // Get the sparsity pattern
  CRSSparsity sp = sparsity().getSub(ii,jj,mapping);
 
  // Create return MX
  MX ret;
  ret.assignNode(new Mapping(sp));
  ret->addDependency(*this,mapping);
  return ret;
}

void MX::setSub(const std::vector<int>& ii, const std::vector<int>& jj, const MX& el){
  casadi_assert_message(ii.size()==el.size1(),"right hand size must match dimension of left hand side in assignment");
  casadi_assert_message(jj.size()==el.size2(),"right hand size must match dimension of left hand side in assignment");
  if(dense() && el.dense()){
    // Dense mode
    for(int i=0; i<ii.size(); ++i) {
      for(int j=0; j<jj.size(); ++j) {
        (*this)[ii[i]*size2() + jj[j]]=el[i*el.size2()+j];
      }
    }
  } else {
    // Sparse mode

    // Remove submatrix to be replaced
    erase(ii,jj);

    // Extend el to the same dimension as this
    MX el_ext = el;
    el_ext.enlarge(size1(),size2(),ii,jj);

    // Unite the sparsity patterns
    *this = unite(*this,el_ext);
  }
}

MX MX::getNZ(const std::vector<int>& kk) const{
  CRSSparsity sp(kk.size(),1,true);
  MX ret;
  ret.assignNode(new Mapping(sp));
  ret->addDependency(*this,kk);
  return ret;
}

void MX::setNZ(const std::vector<int>& kk, const MX& el){
  MX ret;
  ret.assignNode(new Mapping(sparsity()));
  ret->addDependency(*this,range(size()));
  ret->addDependency(el,range(kk.size()),kk);
  *this = ret;
}

const MX MX::operator()(int i, int j) const{
  const CRSSparsity& sp = CRSSparsity::scalarSparsity;
  int ind = sparsity().getNZ(i,j);
  casadi_assert(ind>=0);
  
  MX ret;
  ret.assignNode(new Mapping(sp));
  ret->addDependency(*this,vector<int>(1,ind));
  return ret;
}

const MX MX::operator[](int k) const{
  return getNZ(vector<int>(1,k));
}

const MX MX::operator()(const std::vector<int>& ii, const std::vector<int>& jj) const{
  return getSub(ii,jj);
}

SubMatrix<MX > MX::operator()(const std::vector<int>& ii, const std::vector<int>& jj){
  return SubMatrix<MX>(*this,ii, jj);
}

SubMatrix<MX> MX::operator()(int i, int j){
  return operator()(vector<int>(1,i),vector<int>(1,j));
}

const MX MX::operator()(const std::vector<int>& ii, int j) const{
  return operator()(ii,vector<int>(1,j));
}

const MX MX::operator()(int i, const std::vector<int>& jj) const{
  return operator()(vector<int>(1,i),jj);
}

SubMatrix<MX > MX::operator()(const std::vector<int>& ii, int j){
  return operator()(ii, vector<int>(1,j));
}

SubMatrix<MX > MX::operator()(int i, const std::vector<int>& jj){
  return operator()(vector<int>(1,i), jj);
}

NonZeros<MX> MX::operator[](int k){
  return NonZeros<MX>(*this,vector<int>(1,k));
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
  return MX::binary(ADD,x,y);
}

MX operator-(const MX &x, const MX &y){
  return MX::binary(SUB,x,y);
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
  ret.assignNode(new UnaryOp(Operation(op),x));
  return ret;
}

MX MX::scalar_matrix(int op, const MX &x, const MX &y){
  MX ret;
  ret.assignNode(new ScalarMatrixOp(Operation(op),x,y));  
  return ret;
}

MX MX::matrix_scalar(int op, const MX &x, const MX &y){
  MX ret;
  ret.assignNode(new MatrixScalarOp(Operation(op),x,y));  
  return ret;
}

MX MX::matrix_matrix(int op, const MX &x, const MX &y){
  MX ret;
  ret.assignNode(new MatrixMatrixOp(Operation(op),x,y)); 
  return ret;
}

MX operator*(const MX &x, const MX &y){
  return MX::binary(MUL,x,y);
}

MX operator/(const MX &x, const MX &y){
  return MX::binary(DIV,x,y);
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

bool MX::dense() const{
  return numel()==size();
}

MX MX::zeros(int nrow, int ncol){
  return MX(nrow,ncol);
}

MX MX::ones(int nrow, int ncol){
  return MX(Matrix<double>(nrow,ncol,1));
}

MX MX::eye(int n){
  Matrix<double> I(n,n);
  for(int i=0; i<n; ++i)
    I(i,i) = 1;
  return MX(I);
}

MX MX::operator-() const{
  return unary(NEG,*this);
}

const CRSSparsity& MX::sparsity() const{
  return (*this)->sparsity();
}

void MX::erase(const std::vector<int>& ii, const std::vector<int>& jj){
  // Get sparsity of the new matrix
  CRSSparsity sp = sparsity();
  
  // Erase from sparsity pattern
  std::vector<int> mapping = sp.erase(ii,jj);
  
  // Create new matrix
  if(mapping.size()!=size()){
    MX ret;
    ret.assignNode(new Mapping(sp));
    ret->addDependency(*this,mapping);
    *this = ret;
  }
}

void MX::enlarge(int nrow, int ncol, const std::vector<int>& ii, const std::vector<int>& jj){
  CRSSparsity sp = sparsity();
  sp.enlarge(nrow,ncol,ii,jj);
  
  MX ret;
  ret.assignNode(new Mapping(sp));
  ret->addDependency(*this,range(size()));
  
  *this = ret;
}

MX::MX(int nrow, int ncol, const MX& val){
  // Make sure that val is scalar
  casadi_assert(val.numel()==1);
  casadi_assert(val.size()==1);
  
  CRSSparsity sp(nrow,ncol,true);
  assignNode(new Mapping(sp));
  (*this)->addDependency(val,vector<int>(sp.size(),0));
}

} // namespace CasADi


using namespace CasADi;
namespace std{

MX exp(const MX& x){
  return MX::unary(EXP,x);
}

MX log(const MX& x){
  return MX::unary(LOG,x);
}

MX sqrt(const MX& x){
  return MX::unary(SQRT,x);
}

MX sin(const MX& x){
  return MX::unary(SIN,x);
}

MX cos(const MX& x){;  
  return MX::unary(COS,x);
}

MX tan(const MX& x){
  return MX::unary(TAN,x);
}

MX atan(const MX& x){
  return MX::unary(ATAN,x);
}

MX asin(const MX& x){
  return MX::unary(ASIN,x);
}

MX acos(const MX& x){
  return MX::unary(ACOS,x);
}

MX pow(const MX& x, const MX& n){
  return MX::binary(POW,x,n);
}

MX erf(const MX& x){
  return MX::unary(ERF,x);
}

MX floor(const MX& x){
  return MX::unary(FLOOR,x);
}

MX ceil(const MX& x){
  return MX::unary(CEIL,x);
}

MX fmin(const MX& x, const MX& y){
  return MX::binary(FMIN,x,y);
}

MX fmax(const MX& x, const MX& y){
  return MX::binary(FMAX,x,y);
}



} // namespace std
