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
#include "mx_tools.hpp"
#include "unary_mx.hpp"
#include "binary_mx.hpp"
#include "mapping.hpp"
#include "../fx/sx_function.hpp"
#include "evaluation_mx.hpp"
#include "symbolic_mx.hpp"
#include "constant_mx.hpp"
#include "mx_tools.hpp"
#include "../stl_vector_tools.hpp"
#include "../matrix/matrix_tools.hpp"
#include "multiplication.hpp"
#include "densification.hpp"
#include "norm.hpp"
#include "../casadi_math.hpp"

using namespace std;
namespace CasADi{

MX::~MX(){
}

MX::MX(){
}

MX::MX(double x){
  assignNode(new ConstantMX(x));
}

MX::MX(const Matrix<double> &x){
  assignNode(new ConstantMX(x));
}

MX::MX(const vector<double> &x){
  assignNode(new ConstantMX(x));
}

MX::MX(const string& name, int n, int m){
  assignNode(new SymbolicMX(name,n,m));
}

MX::MX(const std::string& name,const std::pair<int,int> &nm) {
  assignNode(new SymbolicMX(name,nm.first,nm.second));
}

MX::MX(const string& name, const CRSSparsity & sp){
  assignNode(new SymbolicMX(name,sp));
}

MX::MX(int nrow, int ncol){
  assignNode(new Mapping(CRSSparsity(nrow,ncol)));
}

MX::MX(const CRSSparsity& sp, const MX& val){
  // Make sure that val is dense and scalar
  casadi_assert(val.scalar());
  
  // Dense matrix if val dense
  if(val.dense()){
    // Create mapping
    assignNode(new Mapping(sp));
    (*this)->assign(val,vector<int>(sp.size(),0));
    simplifyMapping(*this);
  } else {
    // Empty matrix
    *this = sparse(sp.size1(),sp.size2());
  }
}

MX MX::create(MXNode* node){
  MX ret;
  ret.assignNode(node);
  return ret;
}

bool MX::__nonzero__() const {
  if (isNull()) {casadi_error("Cannot determine truth value of null MX.");}
  return (operator->())->__nonzero__();
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
  MX ret = MX::create(new Mapping(sp));
  ret->assign(*this,mapping);
  simplifyMapping(ret);
  return ret;
}

const MX MX::getSub(const Matrix<int>& k, int dummy) const {
  MX ret = MX::create(new Mapping(k.sparsity()));
  int s = size();
  const std::vector<int> & d = k.data();
  for (int i=0;i< k.size();i++) {
    if (d[i]>=s) casadi_error("MX::getSub: a non-zero element at position " << d[i] << " was requested, but MX is only " << dimString());
  }
  ret->assign(*this,k.data());
  simplifyMapping(ret);
  return ret;
}

const MX MX::getSub(int i, int j) const{
  int ind = sparsity().getNZ(i,j);

  MX ret;
  if (ind>=0) {
    const CRSSparsity& sp = CRSSparsity::scalarSparsity;
    ret.assignNode(new Mapping(sp));
    ret->assign(*this,vector<int>(1,ind));
  } else {
    const CRSSparsity& sp = CRSSparsity::scalarSparsitySparse;
    ret.assignNode(new Mapping(sp));
    ret->assign(*this,vector<int>(0));
  }
  simplifyMapping(ret);
  return ret;
}

const MX MX::getSub(const std::vector<int>& ii, const Matrix<int>& k) const{
  std::vector< int > cols = range(size2());
  std::vector< MX > temp;

  for (int i=0;i<ii.size();++i) {
    MX m(k.sparsity(),MX(0));
    for (int j=0;j<m.size();++j) {
      m[j] = getSub(ii.at(i),k.at(j));
    }
    temp.push_back(m);
  }
  MX ret = vertcat(temp);
  simplifyMapping(ret);
  return ret;
}

const MX MX::getSub(const Matrix<int>& k, const std::vector<int>& jj) const{
  std::vector< int > rows = range(size1());
  std::vector< MX > temp;

  for (int j=0;j<jj.size();++j) {
    MX m(k.sparsity(),MX(0));
    for (int i=0;i<m.size();++i) {
      m[i] = getSub(k.at(i),jj.at(j));
    }
    temp.push_back(m);
  }
  MX ret = horzcat(temp);
  simplifyMapping(ret);
  return ret;
}

const MX MX::getSub(const Matrix<int>& i, const Matrix<int>& j) const {
   casadi_assert_message(i.sparsity()==j.sparsity(),"getSub(Imatrix i, Imatrix j): sparsities must match. Got " << i.dimString() << " and " << j.dimString() << ".");

   MX ret(i.sparsity(),MX(0));
   for (int k=0;k<i.size();++k) {
    ret[k] = getSub(i.at(k),j.at(k));
   }
   simplifyMapping(ret);
   return ret;
}

const MX MX::getSub(const CRSSparsity& sp, int dummy) const {
  casadi_assert_message(size1()==sp.size1() && size2()==sp.size2(),"getSub(CRSSparsity sp): shape mismatch. This matrix has shape " << size1() << " x " << size2() << ", but supplied sparsity index has shape " << sp.size1() << " x " << sp.size2() << "." );
  vector<unsigned char> mappingc; // Mapping that will be filled by patternunion
  
  sparsity().patternUnion(sp,mappingc,true, false, true);
  vector<int> inz;
  vector<int> onz;
  
  int k = 0;     // Flat index into non-zeros of this matrix
  int j = 0;     // Flat index into non-zeros of the resultant matrix
  for (int i=0;i<mappingc.size();++i) {
    if (mappingc[i] & 1) { // If the original matrix has a non-zero entry in the union
      if (!(mappingc[i] & 4)) {
         // If this non-zero entry appears in the intersection, add it to the mapping
         inz.push_back(k);
         onz.push_back(j);
      }
      k++;                 // Increment the original matrix' non-zero index counter
    }
    if (mappingc[i] & 2) j++;
  }

  MX ret = MX::create(new Mapping(sp));
  ret->assign(*this,inz,onz);

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

void MX::setSub(const Matrix<int>& k, const MX& el){
  // Allow el to be a 1x1
  if (el.dense() && el.scalar()) {
    if (k.numel()>1) {
      setSub(k,MX(k.sparsity(),el));
      return;
    }
  }
  
  casadi_assert_message(k.sparsity()==el.sparsity(),"Sparsity mismatch." << "lhs is " << k.dimString() << ", while rhs is " << el.dimString());
  
  casadi_error("MX::setSub not implemented yet");
}


void MX::setSub(const vector<int>& ii, const vector<int>& jj, const MX& el){
  // Allow el to be a 1x1
  if (el.dense() && el.scalar()) {
    if (ii.size()>1 || jj.size()>1) {
      setSub(ii,jj,MX(ii.size(),jj.size(),el));
      return;
    }
  }
  
  casadi_assert_message(ii.size()==el.size1(),"Dimension mismatch." << "lhs is " << ii.size() << " x " << jj.size() << ", while rhs is " << el.dimString());
  casadi_assert_message(jj.size()==el.size2(),"Dimension mismatch." << "lhs is " << ii.size() << " x " << jj.size() << ", while rhs is " << el.dimString());

  if(dense() && el.dense()){
    // Dense mode
    int ld = size2(), ld_el = el.size2(); // leading dimensions
    for(int i=0; i<ii.size(); ++i) {
      for(int j=0; j<jj.size(); ++j) {
        (*this)[ii[i]*ld + jj[j]]=el[i*ld_el+j];
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

void MX::setSub(const Matrix<int>& i, const std::vector<int>& jj, const MX& el) {
  // If el is scalar
  if(el.scalar() && (jj.size() > 1 || i.size() > 1)){
    setSub(i,jj,repmat(MX(i.sparsity(),el),1,jj.size()));
    return;
  }

  if (!inBounds(jj,size2())) {
    casadi_error("setSub[i,jj] out of bounds. Your jj contains " << *std::min_element(jj.begin(),jj.end()) << " up to " << *std::max_element(jj.begin(),jj.end()) << ", which is outside of the matrix shape " << dimString() << ".");
  }
  
  //CRSSparsity result_sparsity = repmat(i,1,jj.size()).sparsity();
  CRSSparsity result_sparsity = horzcat(std::vector< Matrix<int> >(jj.size(),i)).sparsity();
  
  casadi_assert_message(result_sparsity == el.sparsity(),"setSub(Imatrix" << i.dimString() << ",Ivector(length=" << jj.size() << "),Matrix<T>)::Dimension mismatch. The sparsity of repmat(Imatrix,1," << jj.size() << ") = " << result_sparsity.dimString()  << " must match the sparsity of MX = "  << el.dimString() << ".");

  std::vector<int> slice_i = range(i.size1());
  
  for(int k=0; k<jj.size(); ++k) {
     MX el_k = el(slice_i,range(k*i.size2(),(k+1)*i.size2()));
     for (int j=0;j<i.size();++j) {
       (*this)(i.at(j),jj[k])=el_k[j];
     }
  }
  
}



void MX::setSub(const std::vector<int>& ii, const Matrix<int>& j, const MX& el) {
  // If el is scalar
  if(el.scalar() && (ii.size() > 1 || j.size() > 1)){
    setSub(ii,j,repmat(MX(j.sparsity(),el),ii.size(),1));
    return;
  }

  if (!inBounds(ii,size1())) {
    casadi_error("setSub[ii,j] out of bounds. Your ii contains " << *std::min_element(ii.begin(),ii.end()) << " up to " << *std::max_element(ii.begin(),ii.end()) << ", which is outside of the matrix shape " << dimString() << ".");
  }
  
  //CRSSparsity result_sparsity = repmat(j,ii.size(),1).sparsity();
  CRSSparsity result_sparsity = vertcat(std::vector< Matrix<int> >(ii.size(),j)).sparsity();
  
  casadi_assert_message(result_sparsity == el.sparsity(),"setSub(Ivector(length=" << ii.size() << "),Imatrix" << j.dimString() << ",MX)::Dimension mismatch. The sparsity of repmat(Imatrix," << ii.size() << ",1) = " << result_sparsity.dimString() << " must match the sparsity of Matrix<T> = " << el.dimString() << ".");
  
  std::vector<int> slice_j = range(j.size2());
  
  for(int k=0; k<ii.size(); ++k) {
     MX el_k = el(range(k*j.size1(),(k+1)*j.size1()),slice_j);
     for (int i=0;i<j.size();++i) {
       (*this)(ii[k],j.at(i))=el_k[i];
     }
  }
  
}


void MX::setSub(const Matrix<int>& i, const Matrix<int>& j, const MX& el) {
   casadi_assert_message(i.sparsity()==j.sparsity(),"setSub(Imatrix i, Imatrix j, Imatrix el): sparsities must match. Got " << i.dimString() << " for i and " << j.dimString() << " for j.");

  // If el is scalar
  if(el.scalar() && i.numel() > 1){
    setSub(i,j,MX(i.sparsity(),el));
    return;
  }
  
  casadi_assert_message(el.sparsity()==i.sparsity(),"setSub(Imatrix i, Imatrix j, MX el): sparsities must match. Got " << el.dimString() << " for el and " << j.dimString() << " for i and j.");
  
  for(int k=0; k<i.size(); ++k) {
     (*this)(i.at(k),j.at(k)) = el[k]; 
  }
}

void MX::setSub(const CRSSparsity& sp, int dummy, const MX& el) {
   casadi_assert_message(size1()==sp.size1() && size2()==sp.size2(),"getSub(CRSSparsity sp): shape mismatch. This matrix has shape " << size1() << " x " << size2() << ", but supplied sparsity index has shape " << sp.size1() << " x " << sp.size2() << "." );
   casadi_error("Not implemented yet");
}

MX MX::getNZ(int k) const{
  if (k<0) k+=size();
  casadi_assert_message(k<size(),"MX::getNZ: requested at(" <<  k << "), but that is out of bounds:  " << dimString() << ".");
  return getNZ(vector<int>(1,k));
}

MX MX::getNZ(const vector<int>& k) const{
  CRSSparsity sp(k.size(),1,true);
  MX ret;
  
  for (int i=0;i<k.size();i++) {
    casadi_assert_message(k[i] < size(),"Mapping::assign: index vector reaches " << k[i] << ", while dependant is only of size " << size());
  }
  
  ret.assignNode(new Mapping(sp));
  ret->assign(*this,k);
  simplifyMapping(ret);
  return ret;
}

MX MX::getNZ(const Matrix<int>& k) const{
  CRSSparsity sp(k.size(),1,true);
  MX ret;

  ret.assignNode(new Mapping(sp));
  ret->assign(*this,k.data());
  
  return ret;
}

void MX::setNZ(int k, const MX& el){
  if (k<0) k+=size();
  casadi_assert_message(k<size(),"MX::setNZ: requested at(" <<  k << "), but that is out of bounds:  " << dimString() << ".");
  setNZ(vector<int>(1,k),el);
}

void MX::setNZ(const vector<int>& k, const MX& el){
  casadi_assert_message(k.size()==el.size() || el.size()==1,
    "MX::setNZ: length of non-zero indices (" << k.size() << ") " <<
    "must match size of rhs (" << el.size() << ")."
  );
  
  // Call recursively if points both objects point to the same node
  if(this==&el){
    MX el2 = el;
    setNZ(k,el2);
    return;
  }
  
  // Assert correctness
  for(int i=0; i<k.size(); ++i){
    casadi_assert_message(k[i] < size(), "Mapping::assign: index vector reaches " << k[i] << ", while dependant is only of size " << size());
  }
  
  // Quick return if no assignments to be made
  if(k.empty()) return;

  // Make mapping if not already so
  if(!isMapping()){
    MX x;
    x.assignNode(new Mapping(sparsity()));
    x->assign(*this,range(size()));
    *this = x;
  }
  
  // Make sure that the assignment does not affect other nodes
  makeUnique(false);
  
  // Input nonzeros (see Mapping class)
  std::vector<int> inz(k.size(),0);
  if(!(el.scalar() && el.dense())){
    for(int i=0; i<inz.size(); ++i) inz[i]=i;
  }
  
  // Add to mapping
  (*this)->assign(el,inz,k);
  simplifyMapping(*this);
}

void MX::setNZ(const Matrix<int>& kk, const MX& m){
  if (m.size()==1 && m.numel()==1) {
    setNZ(kk.data(),m);
    return;
  }
  casadi_assert_message(kk.sparsity()==m.sparsity(),"Matrix<T>::setNZ: sparsity of IMatrix index " << kk.dimString() << " " << std::endl << "must match sparsity of rhs " << m.dimString() << ".");
  setNZ(kk.data(),m);
}

const MX MX::at(int k) const {
  return getNZ(k); 
}

/// Access a non-zero element
NonZeros<MX,int> MX::at(int k) {
  return NonZeros<MX,int>(*this,k);
}

MX MX::binary(int op, const MX &x, const MX &y){
  // Make sure that dimensions match
  casadi_assert_message((x.scalar() || y.scalar() || (x.size1()==y.size1() && x.size2()==y.size2())),"Dimension mismatch." << "lhs is " << x.dimString() << ", while rhs is " << y.dimString());
  
  // Quick return if zero
  if((operation_checker<F0XChecker>(op) && isZero(x)) || 
    (operation_checker<FX0Checker>(op) && isZero(y))){
    return sparse(std::max(x.size1(),y.size1()),std::max(x.size2(),y.size2()));
  }
  
  // Create binary node
  if(x.scalar())
    return scalar_matrix(op,x,y);
  else if(y.scalar())  
    return matrix_scalar(op,x,y);
  else
    return matrix_matrix(op,x,y);
}

MX MX::unary(int op, const MX &x){
  // Quick return if zero
  if(operation_checker<F0XChecker>(op) && isZero(x)){
    return sparse(x.size1(),x.size2());
  } else {
    return create(new UnaryMX(Operation(op),x));
  }
}

MX MX::scalar_matrix(int op, const MX &x, const MX &y){
  // Check if the scalar is sparse (i.e. zero)
  if(x.size()==0){
    return scalar_matrix(op,0,y);
  } else {
    // Check if it is ok to loop over nonzeros only
    if(y.dense() || operation_checker<FX0Checker>(op)){
      // Loop over nonzeros
      return create(new ScalarNonzerosOp(Operation(op),x,y));
    } else {
      // Put a densification node in between
      return scalar_matrix(op,x,densify(y));
    }
  }
}

MX MX::matrix_scalar(int op, const MX &x, const MX &y){
  // Check if the scalar is sparse (i.e. zero)
  if(y.size()==0){
    return matrix_scalar(op,x,0);
  } else {
    // Check if it is ok to loop over nonzeros only
    if(x.dense() || operation_checker<F0XChecker>(op)){
      // Loop over nonzeros
      return create(new NonzerosScalarOp(Operation(op),x,y));
    } else {
      // Put a densification node in between
      return matrix_scalar(op,densify(x),y);
    }
  }
}

MX MX::matrix_matrix(int op, const MX &x, const MX &y){
  // Check if we can carry out the operation only on the nonzeros
  if((x.dense() && y.dense()) ||
     (operation_checker<F00Checker>(op) && x.sparsity()==y.sparsity())){
    // Loop over nonzeros only
    return create(new NonzerosNonzerosOp(Operation(op),x,y)); 
  } else {
    // Sparse matrix-matrix operation necessary
    return create(new SparseSparseOp(Operation(op),x,y)); 
  }
}

MXNode* MX::operator->(){
  return (MXNode*)SharedObject::operator->();
}

const MXNode* MX::operator->() const{
  return (const MXNode*)SharedObject::operator->();
}

MX MX::repmat(const MX& x, const std::pair<int, int> &nm){
  return repmat(x,nm.first,nm.second);
}

MX MX::repmat(const MX& x, int nrow, int ncol){
  if(x.scalar()){
    return MX(nrow,ncol,x);
  } else {
    casadi_assert_message(0,"not implemented");
    return MX();
  }
}

MX MX::sparse(int nrow, int ncol){
  return MX(nrow,ncol);
}

MX MX::sparse(const std::pair<int, int> &nm){
  return sparse(nm.first,nm.second);
}

MX MX::zeros(int nrow, int ncol){
  return MX(Matrix<double>::zeros(nrow,ncol));
}

MX MX::zeros(const std::pair<int, int> &nm){
  return zeros(nm.first,nm.second);
}

MX MX::ones(int nrow, int ncol){
  return MX(Matrix<double>::ones(nrow,ncol));
}

MX MX::ones(const std::pair<int, int> &nm){
  return ones(nm.first,nm.second);
}

MX MX::inf(int nrow, int ncol){
  return MX(Matrix<double>::inf(nrow,ncol));
}

MX MX::inf(const std::pair<int, int> &nm){
  return inf(nm.first,nm.second);
}

MX MX::nan(int nrow, int ncol){
  return MX(Matrix<double>::nan(nrow,ncol));
}

MX MX::nan(const std::pair<int, int> &nm){
  return nan(nm.first,nm.second);
}

MX MX::eye(int n){
  Matrix<double> I(CRSSparsity::createDiagonal(n),1);
  return MX(I);
}


MX MX::operator-() const{
  if((*this)->isOperation(OP_NEG)){
    return (*this)->dep(0);
  } else {
    return unary(OP_NEG,*this);
  }
}

MX::MX(const MX& x) : CachedObject(x){
}

const CRSSparsity& MX::sparsity() const{
  return (*this)->sparsity();
}

CRSSparsity& MX::sparsityRef(){
  // Since we can potentially change the behavior of the MX node, we must make a deep copy if there are other references
  makeUnique();
  
  // Return the reference, again, deep copy if multiple references
  (*this)->sparsity_.makeUnique();
  return (*this)->sparsity_;
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
    ret->assign(*this,mapping);
    simplifyMapping(ret);
    *this = ret;
  }
}

void MX::enlarge(int nrow, int ncol, const vector<int>& ii, const vector<int>& jj){
  CRSSparsity sp = sparsity();
  sp.enlarge(nrow,ncol,ii,jj);
  
  MX ret;
  ret.assignNode(new Mapping(sp));
  ret->assign(*this,range(size()));
  simplifyMapping(ret);

  *this = ret;
}

MX::MX(int nrow, int ncol, const MX& val){
  // Make sure that val is scalar
  casadi_assert(val.scalar());
  casadi_assert(val.dense());
  
  CRSSparsity sp(nrow,ncol,true);
  assignNode(new Mapping(sp));
  (*this)->assign(val,vector<int>(sp.size(),0));
}

std::string MX::dimString() const {
  std::stringstream ss;
  ss << "(" << size1() << "x" << size2() << "=" << numel() << "|" << size() << ")";
  return ss.str();
}

MX MX::mul(const MX& y) const{
  const MX& x = *this;

  // Check if we can simplify the product
  if(isIdentity(x)){
    return y;
  } else if(isIdentity(y)){
    return x;
  } else if(isZero(x) || isZero(y)){
    // See if one of the arguments can be used as result
    if(x.size()==0 && y.size1()==y.size2())
      return x;
    else if(y.size()==0 && x.size1()==x.size2())
      return y;
    else
      return MX::sparse(x.size1(),y.size2());
  } else if(x.scalar() || y.scalar()){
    return x*y;
  } else if(x.sparsity().diagonal() && y.size2()==1){
    return diag(x)*y;
  } else if(y.sparsity().diagonal() && x.size1()==1){
    return x*trans(diag(y));
  } else {
    return MX::create(new Multiplication(x,trans(y)));
  }
}

MX MX::inner_prod(const MX& y) const{
  const MX& x = *this;
  casadi_assert_message(x.size2()==1,"inner_prod: first factor not a vector");
  casadi_assert_message(y.size2()==1, "inner_prod: second factor not a vector");
  casadi_assert_message(x.size1()==y.size1(),"inner_prod: dimension mismatch");
  MX sum = 0;
  for(int i=0; i<x.size1(); ++i)
    sum += x(i)*y(i);
  return sum;
}

MX MX::outer_prod(const MX& y) const{
  return mul(trans(y));
}

MX MX::__pow__(const MX& n) const{
  if(n->isConstant()){
    return MX::binary(OP_CONSTPOW,*this,n);
  } else {
    return MX::binary(OP_POW,*this,n);
  }
}

MX MX::constpow(const MX& b) const{    
  return binary(OP_CONSTPOW,*this,b);
}

MX MX::fmin(const MX& b) const{       
  return binary(OP_FMIN,*this,b);
}

MX MX::fmax(const MX& b) const{       
  return binary(OP_FMAX,*this,b);
}

MX MX::arctan2(const MX& b) const{       
  return binary(OP_ATAN2,*this,b);
}

MX MX::printme(const MX& b) const{ 
  return binary(OP_PRINTME,*this,b);
}

MX MX::exp() const{ 
  return unary(OP_EXP,*this);
}

MX MX::log() const{ 
  return unary(OP_LOG,*this);
}

MX MX::log10() const{ 
  return log()*(1/std::log(10.));
}

MX MX::sqrt() const{ 
  return unary(OP_SQRT,*this);
}

MX MX::sin() const{ 
  return unary(OP_SIN,*this);
}

MX MX::cos() const{ 
  return unary(OP_COS,*this);
}

MX MX::tan() const{ 
  return unary(OP_TAN,*this);
}

MX MX::arcsin() const{ 
  return unary(OP_ASIN,*this);
}

MX MX::arccos() const{ 
  return unary(OP_ACOS,*this);
}

MX MX::arctan() const{ 
  return unary(OP_ATAN,*this);
}

MX MX::sinh() const{ 
  return unary(OP_SINH,*this);
}

MX MX::cosh() const{ 
  return unary(OP_COSH,*this);
}

MX MX::tanh() const{ 
  return unary(OP_TANH,*this);
}

MX MX::arcsinh() const{ 
  return unary(OP_ASINH,*this);
}

MX MX::arccosh() const{ 
  return unary(OP_ACOSH,*this);
}

MX MX::arctanh() const{ 
  return unary(OP_ATANH,*this);
}

MX MX::floor() const{ 
  return unary(OP_FLOOR,*this);
}

MX MX::ceil() const{ 
  return unary(OP_CEIL,*this);
}

MX MX::fabs() const{ 
  return unary(OP_FABS,*this);
}

MX MX::sign() const{ 
  return unary(OP_SIGN,*this);
}

MX MX::erfinv() const{ 
  return unary(OP_ERFINV,*this);
}

MX MX::erf() const{ 
  return unary(OP_ERF,*this);
}

MX MX::logic_not() const{ 
  return unary(OP_NOT,*this);
}

MX MX::__add__(const MX& y) const{
  const MX& x = *this;
  bool samedim = x.size1()==y.size1() && x.size2()==y.size2();
  if((samedim || x.scalar()) && isZero(x)){
    return y;
  } else if((samedim || y.scalar()) && isZero(y)){
    return x;
  } else if(y->isOperation(OP_NEG)){
    return x - y->dep(0);
  } else if(x->isOperation(OP_NEG)){
    return y - x->dep(0);
  } else if(x->isOperation(OP_SUB) && y.get()==x->dep(1).get()){
    return x->dep(0);
  } else if(y->isOperation(OP_SUB) && x.get()==y->dep(1).get()){
    return y->dep(0);
  } else {
    return MX::binary(OP_ADD,x,y);
  }
}

MX MX::__sub__(const MX& y) const{
  const MX& x = *this;
  bool samedim = x.size1()==y.size1() && x.size2()==y.size2();
  if((samedim || x.scalar()) && isZero(x)){
    return -y;
  } else if((samedim || y.scalar()) && isZero(y)){
    return x;
  } else if(y->isOperation(OP_NEG)){
    return x+y->dep(0);
  } else if(y.get()==x.get()){
    return MX::sparse(x.size1(),x.size2());
  } else {
    return MX::binary(OP_SUB,x,y);
  }
}

MX MX::__mul__(const MX& y) const{
  const MX& x = *this;
  bool samedim = x.size1()==y.size1() && x.size2()==y.size2();
  if((samedim || x.scalar()) && isOne(x)){
    return y;
  } else if((samedim || x.scalar()) && isMinusOne(x)){
    return -y;
  } else if((samedim || y.scalar()) && isOne(y)){
    return x;
  } else if((samedim || y.scalar()) && isMinusOne(y)){
    return -x;
  } else {
    return MX::binary(OP_MUL,x,y);
  }
}

MX MX::__div__(const MX& y) const{
  const MX& x = *this;
  bool samedim = x.size1()==y.size1() && x.size2()==y.size2();
  if((samedim || y.scalar()) && isOne(y)){
    return x;
  } else {
    return MX::binary(OP_DIV,x,y);
  }
}

MX MX::__lt__(const MX& y) const{
  return MX::binary(OP_LT,*this,y);
}
  
MX MX::__le__(const MX& y) const{
  return MX::binary(OP_LE,*this,y);
}
  
MX MX::__eq__(const MX& y) const{
  return MX::binary(OP_EQ,*this,y);
}
  
MX MX::__ne__(const MX& y) const{
  return MX::binary(OP_NE,*this,y);
}
  
MX MX::logic_and(const MX& y) const{
  return MX::binary(OP_AND,*this,y);
}
  
MX MX::logic_or(const MX& y) const{
  return MX::binary(OP_OR,*this,y);
}
  
MX MX::if_else_zero(const MX& y) const{
  return MX::binary(OP_IF_ELSE_ZERO,*this,y);
}
  
MX MX::__constpow__(const MX& b) const { return (*this).constpow(b);}
MX MX::__mrdivide__(const MX& b) const { if (b.scalar()) return *this/b; throw CasadiException("mrdivide: Not implemented");}
MX MX::__mpower__(const MX& b) const   { return pow(*this,b); throw CasadiException("mpower: Not implemented");}

void MX::append(const MX& y){
  *this = vertcat(*this,y);
}

long MX::max_num_calls_in_print_ = 10000;

void MX::setMaxNumCallsInPrint(long num){
  max_num_calls_in_print_ = num;
}

long MX::getMaxNumCallsInPrint(){
  return max_num_calls_in_print_;
}

MX MX::getDep(int ch) const { return !isNull() ? (*this)->dep(ch) : MX(); }

int MX::getNdeps() const { return !isNull() ? (*this)->ndep() : 0; }
  
std::string MX::getName() const { return !isNull() ? (*this)->getName() : "null"; }

bool 	MX::isSymbolic () const { return !isNull() ? (*this)->isSymbolic() : false; }
bool 	MX::isConstant () const { return !isNull() ? (*this)->isConstant() : false; }
bool 	MX::isMapping () const { return !isNull() ? (*this)->isMapping() : false; }
bool 	MX::isEvaluation () const { return !isNull() ? (*this)->isEvaluation() : false; }
bool 	MX::isEvaluationOutput () const { return !isNull() ? (*this)->isOutputNode() : false; }

int 	MX::getEvaluationOutput () const { return !isNull() ? (*this)->getFunctionOutput() : -1; }

    
bool 	MX::isOperation (int op) const { return !isNull() ? (*this)->isOperation(op) : false; }
bool 	MX::isMultiplication () const { return !isNull() ? (*this)->isMultiplication() : false; }

bool 	MX::isNorm () const { return !isNull() ? dynamic_cast<const Norm*>(get())!=0 : false; }

bool 	MX::isDensification () const { return !isNull() ? dynamic_cast<const Densification*>(get())!=0 : false; }


FX MX::getFunction () {  return (*this)->getFunction(); }
 	
const Matrix<double> & MX::getConstant() const {
  casadi_assert_message(isConstant(),"MX::getConstant: must be constant");
  return dynamic_cast<const ConstantMX*>(get())->x_;
}
 	  
bool MX::isBinary() const { return !isNull() ? dynamic_cast<const BinaryMX*>(get()) != 0 : false;  }

bool MX::isUnary() const { return !isNull() ? dynamic_cast<const UnaryMX*>(get()) != 0 : false;  }
 	
int MX::getOp() const {
  casadi_assert_message(isBinary() || isUnary(),"MX::getOp: must be binary or unary operation");
  if (isBinary()) {
    return dynamic_cast<const BinaryMX*>(get())->op_;
  } else {
    return dynamic_cast<const UnaryMX*>(get())->op_;
  }
}
 	
bool MX::isCommutative() const {
  if (isUnary()) return true;
  casadi_assert_message(isBinary() || isUnary(),"MX::isCommutative: must be binary or unary operation");
  return operation_checker<CommChecker>(getOp());
}

long MX::__hash__() const {
   if (isNull()) return 0;
   // FIXME: This is bad coding, the pointer might on certain architectures be larger than the long,
   // giving an error of type "error: cast from 'const SharedObjectInternal*' to 'int' loses precision"
   // The solution is to convert to intptr_t, but this is not yet C++ standard
   return long(get());
}

Matrix<int> MX::mapping(int iind) const {
  const Mapping * m = dynamic_cast<const Mapping*>(get());
  casadi_assert_message(m!=0, "mapping: argument MX should point to a Mapping node");
  return m->mapping(iind);
}

std::vector<int> MX::getDepInd() const {
  const Mapping * m = dynamic_cast<const Mapping*>(get());
  casadi_assert_message(m!=0, "mapping: argument MX should point to a Mapping node");
  return m->getDepInd();
}

int MX::getTemp() const{
  return (*this)->temp;
}
    
void MX::setTemp(int t){
  (*this)->temp = t;
}

 	  
} // namespace CasADi
