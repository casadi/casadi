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

MX::MX(const vector<double> &x){
  assignNode(new MXConstant(x));
}

MX::MX(const string& name, int n, int m){
  assignNode(new SymbolicMatrix(name,n,m));
}

MX::MX(const string& name, const CRSSparsity & sp){
  assignNode(new SymbolicMatrix(name,sp));
}

MX::MX(int nrow, int ncol){
  assignNode(new Mapping(CRSSparsity(nrow,ncol)));
}

MX MX::create(MXNode* node){
  MX ret;
  ret.assignNode(node);
  return ret;
}

const MX MX::getSub(int i, const std::vector<int>& j) const{
  return getSub(vector<int>(1,i),j);
}

const MX MX::getSub(const std::vector<int>& i, int j) const{
  return getSub(i,vector<int>(1,j));
}

const MX MX::getSub(const vector<int>& ii, const vector<int>& jj) const{
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

const MX MX::getSub(int i, int j) const{
  int ind = sparsity().getNZ(i,j);

  MX ret;
  if (ind>=0) {
    const CRSSparsity& sp = CRSSparsity::scalarSparsity;
    ret.assignNode(new Mapping(sp));
    ret->addDependency(*this,vector<int>(1,ind));
  } else {
    const CRSSparsity& sp = CRSSparsity::scalarSparsitySparse;
    ret.assignNode(new Mapping(sp));
    ret->addDependency(*this,vector<int>(0));
  }
  return ret;
}

void MX::setSub(int i, int j, const MX& el){
  setSub(vector<int>(1,i),vector<int>(1,j),el);
}

void MX::setSub(int i, const std::vector<int>& j, const MX& el){
  setSub(vector<int>(1,i),j,el);
}

void MX::setSub(const std::vector<int>& i, int j, const MX& el){
  setSub(i,vector<int>(1,j),el);
}

void MX::setSub(const vector<int>& ii, const vector<int>& jj, const MX& el){
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

MX MX::getNZ(int k) const{
  if (k<0) k+=size();
  if (k>=size()) {
    stringstream ss;
    ss << "MX::getNZ: requested at(" <<  k << "), but that is out of bounds:  " << dimString() << ".";
    throw CasadiException(ss.str());
  }
  return getNZ(vector<int>(1,k));
}

MX MX::getNZ(const vector<int>& k) const{
  CRSSparsity sp(k.size(),1,true);
  MX ret;
  ret.assignNode(new Mapping(sp));
  ret->addDependency(*this,k);
  return ret;
}

void MX::setNZ(int k, const MX& el){
  if (k<0) k+=size();
  if (k>=size()) {
    stringstream ss;
    ss << "MX::setNZ: requested at(" <<  k << "), but that is out of bounds:  " << dimString() << ".";
    throw CasadiException(ss.str());
  }
  setNZ(vector<int>(1,k),el);
}

void MX::setNZ(const vector<int>& k, const MX& el){
  if (k.size()!=el.size() && el.size()!=1) {
    std::stringstream ss;
    ss << "MX::setNZ: length of non-zero indices (" << k.size() << ") " << std::endl;
    ss << "must match size of rhs (" << el.size() << ")." << std::endl;
    throw CasadiException(ss.str());
  }
  MX ret;
  ret.assignNode(new Mapping(sparsity()));
  ret->addDependency(*this,range(size()));
  if (el.size()==1) {
    ret->addDependency(el,std::vector<int>(k.size(),0),k);
  } else {
    ret->addDependency(el,range(k.size()),k);
  }
  *this = ret;
}

const MX MX::at(int k) const {
  return getNZ(k); 
}

/// Access a non-zero element
NonZeros<MX,int> MX::at(int k) {
  return NonZeros<MX,int>(*this,k);
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
  return create(new UnaryOp(Operation(op),x));
}

MX MX::scalar_matrix(int op, const MX &x, const MX &y){
  return create(new ScalarMatrixOp(Operation(op),x,y)); 
}

MX MX::matrix_scalar(int op, const MX &x, const MX &y){
  return create(new MatrixScalarOp(Operation(op),x,y));
}

MX MX::matrix_matrix(int op, const MX &x, const MX &y){
  return create(new MatrixMatrixOp(Operation(op),x,y)); 
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

MX MX::zeros(const CRSSparsity& sparsity){
  return DMatrix(sparsity,0);
}

MX MX::ones(int nrow, int ncol){
  return MX(Matrix<double>(nrow,ncol,1));
}

MX MX::eye(int n){
  Matrix<double> I(CRSSparsity::createDiagonal(n),1);
  return MX(I);
}


MX MX::operator-() const{
  return unary(NEG,*this);
}

const CRSSparsity& MX::sparsity() const{
  return (*this)->sparsity();
}

void MX::erase(const vector<int>& ii, const vector<int>& jj){
  // Get sparsity of the new matrix
  CRSSparsity sp = sparsity();
  
  // Erase from sparsity pattern
  vector<int> mapping = sp.erase(ii,jj);
  
  // Create new matrix
  if(mapping.size()!=size()){
    MX ret;
    ret.assignNode(new Mapping(sp));
    ret->addDependency(*this,mapping);
    *this = ret;
  }
}

void MX::enlarge(int nrow, int ncol, const vector<int>& ii, const vector<int>& jj){
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

Matrix<double> MX::eval(const vector<Matrix<double> >& x){
  casadi_assert_message((*this)->ndep()==x.size(),"incorrect number of arguments");
  return (*this)->eval(x);
}

Matrix<SX> MX::eval(const vector<Matrix<SX> >& x){
  casadi_assert_message((*this)->ndep()==x.size(),"incorrect number of arguments");
  return (*this)->eval(x);
}

MX MX::eval(const vector<MX>& x){
  casadi_assert_message((*this)->ndep()==x.size(),"incorrect number of arguments");
  return (*this)->eval(x);
}

std::string MX::dimString() const {
  std::stringstream ss;
  ss << "(" << size1() << "x" << size2() << "=" << numel() << "|" << size() << ")";
  return ss.str();
}

MX MX::jac(int iind){
  return (*this)->jac(iind);
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
  if(n->isConstant()){
    return MX::binary(CONSTPOW,x,n);
  } else {
    return MX::binary(POW,x,n);
  }
}

MX constpow(const MX& x, const MX& n){
  return MX::binary(CONSTPOW,x,n);
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
