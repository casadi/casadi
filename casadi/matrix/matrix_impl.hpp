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

#ifndef MATRIX_IMPL_HPP
#define MATRIX_IMPL_HPP

// The declaration of the class is in a separate file
#include "matrix.hpp"

namespace CasADi{
// Implementations

template<class T>
const T& Matrix<T>::elem(int i, int j) const{
  int ind = sparsity().getNZ(i,j);
  if(ind==-1)
    return casadi_limits<T>::zero;
  else
    return at(ind);
}

template<class T>
T& Matrix<T>::elem(int i, int j){
  int oldsize = sparsity().size();
  int ind = sparsityRef().getNZ(i,j);
  if(oldsize != sparsity().size())
    data().insert(begin()+ind,0);
  return at(ind);
}

template<class T>
const Matrix<T> Matrix<T>::getSub(int i, int j) const{
  return elem(i,j);
}

template<class T>
const Matrix<T> Matrix<T>::getSub(const std::vector<int>& ii, const std::vector<int>& jj) const{
  // Nonzero mapping from submatrix to full
  std::vector<int> mapping;
  
  // Get the sparsity pattern
  CRSSparsity sp = sparsity().getSub(ii,jj,mapping);

  // Create return object
  Matrix<T> ret(sp);
  
  // Copy nonzeros
  for(int k=0; k<mapping.size(); ++k)
    ret.data()[k] = data()[mapping[k]];
  
  // Return (RVO)
  return ret;
}

template<class T>
const Matrix<T> Matrix<T>::getSub(const std::vector<int>& ii, const Matrix<int>& k) const{
  std::vector< int > cols = range(size2());
  std::vector< Matrix<T> > temp;

  for (int i=0;i<ii.size();++i) {
    Matrix<T> m = k;
    for (int j=0;j<m.size();++j) {
      m.data()[j] = elem(ii.at(i),k.at(j));
    }
    temp.push_back(m);
  }
  
  return vertcat(temp);
}

template<class T>
const Matrix<T> Matrix<T>::getSub(const Matrix<int>& k, const std::vector<int>& jj) const{
  std::vector< int > rows = range(size1());
  std::vector< Matrix<T> > temp;

  for (int j=0;j<jj.size();++j) {
    Matrix<T> m = k;
    for (int i=0;i<m.size();++i) {
      m.data()[i] = elem(k.at(i),jj.at(j));
    }
    temp.push_back(m);
  }
  
  return horzcat(temp);
}

template<class T>
const Matrix<T> Matrix<T>::getSub(const Matrix<int>& i, const Matrix<int>& j) const {
   casadi_assert_message(i.sparsity()==j.sparsity(),"getSub(Imatrix i, Imatrix j): sparsities must match. Got " << i.dimString() << " and " << j.dimString() << ".");

   Matrix<T> ret(i.sparsity());
   for (int k=0;k<i.size();++k) {
    ret.data()[k] = elem(i.at(k),j.at(k));
   }

   return ret;
}

template<class T>
const Matrix<T> Matrix<T>::getSub(const CRSSparsity& sp, int dummy) const {
  casadi_assert_message(size1()==sp.size1() && size2()==sp.size2(),"getSub(CRSSparsity sp): shape mismatch. This matrix has shape " << size1() << " x " << size2() << ", but supplied sparsity index has shape " << sp.size1() << " x " << sp.size2() << "." );
  Matrix<T> ret(sp);

  std::vector<unsigned char> mapping; // Mapping that will be filled by patternunion

  sparsity().patternUnion(sp,mapping,true, false, true);

  int k = 0;     // Flat index into non-zeros of this matrix
  int j = 0;     // Flat index into non-zeros of the resultant matrix;
  for (int i=0;i<mapping.size();++i) {
    if (mapping[i] & 1) { // If the original matrix has a non-zero entry in the union
      if (!(mapping[i] & 4)) ret[j] = data()[k]; // If this non-zero entry appears in the intersection, add it to the mapping
      k++;                 // Increment the original matrix' non-zero index counter
    }
    if (mapping[i] & 2) j++;
  }

  return ret;
}

template<class T>
void Matrix<T>::setSub(int i, int j, const Matrix<T>& el){
  if(el.dense()){
    elem(i,j) = el.toScalar();
  } else {
    setSub(std::vector<int>(1,i),std::vector<int>(1,j),el);
  }
}

template<class T>
void Matrix<T>::setSub(const std::vector<int>& ii, const std::vector<int>& jj, const Matrix<T>& el){
  casadi_assert_message(el.numel()==1 || (ii.size() == el.size1() && jj.size() == el.size2()),"Dimension mismatch." << std::endl << "lhs is " << ii.size() << " x " << jj.size() << ", while rhs is " << el.dimString());
  
  // If m is scalar
  if(el.numel() != ii.size() * jj.size()){
    setSub(ii,jj,Matrix<T>(ii.size(),jj.size(),el.toScalar()));
    return;
  }

  if(dense() && el.dense()){
    // Dense mode
    for(int i=0; i<ii.size(); ++i) {
      for(int j=0; j<jj.size(); ++j) {
        data()[ii[i]*size2() + jj[j]]=el.data()[i*el.size2()+j];
      }
    }
  } else {
    // Sparse mode

    // Remove submatrix to be replaced
    erase(ii,jj);

    // Extend el to the same dimension as this
    Matrix<T> el_ext = el;
    el_ext.enlarge(size1(),size2(),ii,jj);

    // Unite the sparsity patterns
    *this = unite(*this,el_ext);
  }
}

template<class T>
void Matrix<T>::setSub(const Matrix<int>& i, const std::vector<int>& jj, const Matrix<T>& el) {
  // If el is scalar
  if(el.scalar() && (jj.size() > 1 || i.size() > 1)){
    setSub(i,jj,repmat(Matrix<T>(i.sparsity(),el.toScalar()),1,jj.size()));
    return;
  }
  
  CRSSparsity result_sparsity = repmat(i,1,jj.size()).sparsity();
  
  
  casadi_assert_message(result_sparsity == el.sparsity(),"setSub(Imatrix" << i.dimString() << ",Ivector(length=" << jj.size() << "),Matrix<T>)::Dimension mismatch. The sparsity of repmat(Imatrix,1," << jj.size() << ") = " << result_sparsity.dimString()  << " must match the sparsity of Matrix<T> = "  << el.dimString() << ".");
  
  
  std::vector<int> slice_i = range(i.size1());
  
  for(int k=0; k<jj.size(); ++k) {
     Matrix<T> el_k = el(slice_i,range(k*i.size2(),(k+1)*i.size2()));
     for (int j=0;j<i.size();++j) {
       elem(i.at(j),jj[k])=el_k.at(j);
     }
  }
  
}

template<class T>
void Matrix<T>::setSub(const std::vector<int>& ii, const Matrix<int>& j, const Matrix<T>& el) {
  
  // If el is scalar
  if(el.scalar() && (ii.size() > 1 || j.size() > 1)){
    setSub(ii,j,repmat(Matrix<T>(j.sparsity(),el.toScalar()),ii.size(),1));
    return;
  }
  
  CRSSparsity result_sparsity = repmat(j,ii.size(),1).sparsity();
  
  
  casadi_assert_message(result_sparsity == el.sparsity(),"setSub(Ivector(length=" << ii.size() << "),Imatrix" << j.dimString() << ",Matrix<T>)::Dimension mismatch. The sparsity of repmat(Imatrix," << ii.size() << ",1) = " << result_sparsity.dimString() << " must match the sparsity of Matrix<T> = " << el.dimString() << ".");
  
  std::vector<int> slice_j = range(j.size2());
  
  for(int k=0; k<ii.size(); ++k) {
     Matrix<T> el_k = el(range(k*j.size1(),(k+1)*j.size1()),slice_j);
     for (int i=0;i<j.size();++i) {
       elem(ii[k],j.at(i))=el_k.at(i);
     }
  }
  
}


template<class T>
void Matrix<T>::setSub(const Matrix<int>& i, const Matrix<int>& j, const Matrix<T>& el) {
   casadi_assert_message(i.sparsity()==j.sparsity(),"setSub(Imatrix i, Imatrix j, Imatrix el): sparsities must match. Got " << i.dimString() << " for i and " << j.dimString() << " for j.");

  // If el is scalar
  if(el.scalar() && i.numel() > 1){
    setSub(i,j,Matrix<T>(i.sparsity(),el.toScalar()));
    return;
  }
  
  casadi_assert_message(el.sparsity()==i.sparsity(),"setSub(Imatrix i, Imatrix j, Imatrix el): sparsities must match. Got " << el.dimString() << " for el and " << j.dimString() << " for i and j.");
  
  for(int k=0; k<i.size(); ++k) {
     elem(i.at(k),j.at(k)) = el.at(k); 
  }
}

template<class T>
void Matrix<T>::setSub(const CRSSparsity& sp, int dummy, const Matrix<T>& el) {
   casadi_assert_message(size1()==sp.size1() && size2()==sp.size2(),"getSub(CRSSparsity sp): shape mismatch. This matrix has shape " << size1() << " x " << size2() << ", but supplied sparsity index has shape " << sp.size1() << " x " << sp.size2() << "." );
   casadi_error("Not implemented yet");
}

template<class T>
const Matrix<T> Matrix<T>::getNZ(int k) const{
  if (k<0) k+=size();
  return at(k);
}

template<class T>
const Matrix<T> Matrix<T>::getNZ(const std::vector<int>& k) const{
  Matrix<T> ret(k.size(),1,0);
  for(int el=0; el<k.size(); ++el)
    ret.data()[el] = data().at(k[el]);
  
  return ret;
}

template<class T>
const Matrix<T> Matrix<T>::getNZ(const Matrix<int>& k) const{
  Matrix<T> ret(k.sparsity(),0);
  for(int el=0; el<k.size(); ++el)
    ret.data()[el] = data().at(k.at(el));
  
  return ret;
}

template<class T>
void Matrix<T>::setNZ(int k, const Matrix<T>& m){
  if (k<0) k+=size();
  at(k) = m.toScalar();
}

template<class T>
void Matrix<T>::setNZ(const std::vector<int>& kk, const Matrix<T>& m){
  if (m.size()==1 && m.numel()==1) {
    // Allow scalar assignment:
    // m[:2]=3
    for(int k=0; k<kk.size(); ++k)
      data().at(kk[k]) = m.data()[0];
    return;
  }
  
  casadi_assert_message(kk.size()==m.size(),"Matrix<T>::setNZ: length of non-zero indices (" << kk.size() << ") " << std::endl << "must match size of rhs (" << m.size() << ").");
  
  for(int k=0; k<kk.size(); ++k)
    data().at(kk[k]) = m.data()[k];
}

template<class T>
void Matrix<T>::setNZ(const Matrix<int>& kk, const Matrix<T>& m){
  if (m.size()==1 && m.numel()==1) {
    // Allow scalar assignment:
    // m[:2]=3
    for(int k=0; k<kk.size(); ++k)
      data().at(kk.at(k)) = m.data()[0];
    return;
  }
  casadi_assert_message(kk.sparsity()==m.sparsity(),"Matrix<T>::setNZ: sparsity of IMatrix index " << kk.dimString() << " " << std::endl << "must match sparsity of rhs " << m.dimString() << ".");
  
  for(int k=0; k<kk.size(); ++k)
    data().at(kk.at(k)) = m.data()[k];
}

template<class T>
int Matrix<T>::size() const{
  return data().size();
}

template<class T>
int Matrix<T>::sizeU() const{
  return sparsity_.sizeU();
}

template<class T>
int Matrix<T>::sizeL() const{
  return sparsity_.sizeL();
}

template<class T>
int Matrix<T>::size1() const{
  return sparsity_.size1();
}

template<class T>
int Matrix<T>::size2() const{
  return sparsity_.size2();
}

template<class T>
int Matrix<T>::numel() const{
  return size1()*size2();
}

template<class T>
std::pair<int,int> Matrix<T>::shape() const{
  return sparsity_.shape();
}

        
template<class T>
void Matrix<T>::makeDense(int n, int m, const T& val){
  // Quick return if already dense
  if(n*m == size())
    return;
  
  if(size1()!=n || size2()!=m){
    // Also resize
    sparsity_ = CRSSparsity(n,m,true);
    std::fill(data().begin(),data().end(),val);
    data().resize(n*m, val);
  } else {
    // Create a new data vector
    data().resize(n*m,val);
    
    // Loop over the rows in reverse order
    for(int i=n-1; i>=0; --i){
      // Loop over nonzero elements in reverse order
      for(int el=rowind(i+1)-1; el>=rowind(i); --el){
        // Column
        int j = col(el);
        
        // Swap the old position with the new position
        if(el!=j+i*m){
          data()[j+i*m] = data()[el];
          data()[el] = val;
        }
      }
    }
      
    // Save the new sparsity pattern
    sparsity_ = CRSSparsity(n,m,true);
  }
}

template<class T>
bool Matrix<T>::empty() const{
  return numel()==0;
}

template<class T>
bool Matrix<T>::scalar() const{
  return numel()==1;
}

template<class T>
bool Matrix<T>::vector() const{
  return size2()==1;
}

template<class T>
bool Matrix<T>::dense() const{
  return size()==numel();
}

template<class T>
Matrix<T>::Matrix(){
  sparsity_ = CRSSparsity(0,0,false);
}

template<class T>
Matrix<T>::Matrix(const Matrix<T>& m){
  data_ = m.data_;
  sparsity_ = m.sparsity_;
}

template<class T>
Matrix<T>::Matrix(const std::vector<T>& x) : data_(x){
  sparsity_ = CRSSparsity(x.size(),1,true);
}

template<class T>
Matrix<T>::Matrix(const std::vector<T>& x, int n, int m) : data_(x){
  casadi_assert_message(x.size() == n*m, "Dimension mismatch." << std::endl << "You supplied a vector of length " << x.size() << ", but " << n << " x " << m << " = " << n*m);
  sparsity_ = CRSSparsity(n,m,true);
}

template<class T>
Matrix<T>& Matrix<T>::operator=(const Matrix<T>& m){
  sparsity_ = m.sparsity_;
  data_ = m.data_;
  return *this;
}

template<class T>
Matrix<T>::Matrix(int n, int m){
  sparsity_ = CRSSparsity(n,m,false);
}

template<class T>
Matrix<T>::Matrix(int n, int m, const T& val){
  sparsity_ = CRSSparsity(n,m,true);
  data_.resize(n*m, val);
}

template<class T>
void Matrix<T>::makeEmpty(int n, int m){
  sparsity_ = CRSSparsity(n,m,false);
  data().clear();
}

template<class T>
std::string Matrix<T>::className(){ return std::string("Matrix<") + typeName<T>() + std::string(">"); }

template<class T>
void Matrix<T>::printScalar(std::ostream &stream) const {
  casadi_assert_message(numel()==1, "Not a scalar");
  stream << toScalar();
}
  
template<class T>
void Matrix<T>::printVector(std::ostream &stream) const {
  casadi_assert_message(vector(),"Not a vector");
  stream << "[";
  
  // Loop over rows
  for(int i=0; i<size1(); ++i){
    // Add delimitor
    if(i!=0) stream << ",";
    
    // Check if nonzero
    if(rowind(i)==rowind(i+1)){
      stream << "00";
    } else {
      stream << data()[rowind(i)];
    }
  }
  stream << "]";  
}

template<class T>
void Matrix<T>::printMatrix(std::ostream &stream) const{
  for(int i=0; i<size1(); ++i){
    if(i==0)
      stream << "[[";
    else
      stream << " [";
    int j=0;
    int maxj=size2()-1;
    
    for(int el=rowind(i); el<rowind(i+1); ++el){
      // Print leading zeros
      for(;j<col(el); ++j)
        stream << "00" << (j==maxj? " " : ",  ");
      
      // Print element
      stream << data()[el] << (j==maxj? " " : ",  ");
      j++;
    }
    
    // Print trailing zeros
    for(;j<size2(); ++j)
      stream << "00" << (j==maxj? " " : ",  ");
    
    // New row
    if(i==size1()-1)
      stream << "]]" << std::endl;
    else
      stream << "]" << std::endl;
  }
}

template<class T>
void Matrix<T>::printDense(std::ostream &stream) const{
  stream << className() << "(rows = " << size1() << ", cols = " << size2() << "):" << std::endl;
  printMatrix(stream);
}


template<class T>
void Matrix<T>::printSparse(std::ostream &stream) const {
  stream << className() << "(rows = " << size1() << ", cols = " << size2() << ", nnz = " << size() << ")";
  if (size()>0) stream << ":";
  for(int i=0; i<size1(); ++i)
    for(int el=rowind(i); el<rowind(i+1); ++el){
      int j=col(el);
      stream << "[" << i << "," << j << "] -> " << data()[el] << std::endl;
    }
}

template<class T>
void Matrix<T>::print(std::ostream &stream) const{
  if(dense()){
    if (size()==1){
      printScalar(stream);
    } else if (size2()==1){
      printVector(stream);
    } else {
      printDense(stream);
    }
  } else {
    printSparse(stream);
  }
}

template<class T>
void Matrix<T>::repr(std::ostream &stream) const{
  stream << className() << "(";
  if (empty()){
    stream << "[]";
  } else if (numel()==1 && size()==1){
    printScalar(stream);
  } else if (size2()==1){
    printVector(stream);
  } else {
    stream << std::endl;
    printMatrix(stream);
  }
  stream << ")";
}

template<class T>
const std::vector<int>& Matrix<T>::col() const{
  return sparsity().col();
}

template<class T>
const std::vector<int>& Matrix<T>::rowind() const{
  return sparsity_.rowind();
}

template<class T>
int Matrix<T>::col(int el) const{
  return sparsity_.col(el);
}

template<class T>
int Matrix<T>::rowind(int row) const{
  return sparsity_.rowind(row);
}

template<class T>
void Matrix<T>::reserve(int nnz){
  reserve(nnz,size1());
}

template<class T>
void Matrix<T>::reserve(int nnz, int nrow){
  data().reserve(nnz);
  sparsity_.reserve(nnz,nrow);
}

template<class T>
void Matrix<T>::resize(int n_, int m_){
  sparsity_.resize(n_,m_);
}

template<class T>
void Matrix<T>::clear(){
  sparsity_ = CRSSparsity(0,0,false);
  data().clear();
}

template<class T>
Matrix<T>::Matrix(double val){
  sparsity_ = CRSSparsity(1,1,true);
  data().resize(1,val);
}

template<class T>
Matrix<T>::Matrix(int n, int m, const std::vector<int>& col, const std::vector<int>& rowind, const std::vector<T>& d) : data_(d){
  sparsity_ = CRSSparsity(n,m,col,rowind);
  if(data_.size() != sparsity_.size())
    data_.resize(sparsity_.size());
  sanityCheck(true);
}

template<class T>
Matrix<T>::Matrix(const std::vector< std::vector<T> >& d){
  int n=d.size();
  int m=-1;
  for (int i=0;i<n;i++) {
    if (m==-1) {
      m = d[i].size();
    } else {
      casadi_assert_message(m==d[i].size(), 
        "Matrix<T>::Matrix(const std::vector< std::vector<T> >& d): shape mismatch" << std::endl <<
        "Attempting to construct a matrix from a nested list." << std::endl <<
        "I got convinced that the desired size is ("<< n << " x " << m << " ), but now I encounter a vector of size (" << d[i].size() <<  " )" << std::endl
      );
    }
  }
  sparsity_ = CRSSparsity(n,m,true);
  
  data().resize(n*m);

  for (int i=0;i<n;i++) {
    copy(d[i].begin(),d[i].end(),begin()+i*m);
  }
  
}

template<class T>
Matrix<T>::Matrix(const CRSSparsity& sparsity, const T& val){
  sparsity_ = sparsity;
  data().resize(sparsity_.size(),val);
}

template<class T>
Matrix<T>::Matrix(const CRSSparsity& sparsity, const std::vector<T>& d) : data_(d) {
  casadi_assert_message(sparsity.size()==d.size(),"Size mismatch." << std::endl << "You supplied a sparsity of " << sparsity.dimString() << ", but the supplied vector is of length " << d.size());
  sparsity_ = sparsity;
}

template<class T>
void Matrix<T>::setZero(){
  setAll(0);
}

template<class T>
void Matrix<T>::setAll(const T& val){
  std::fill(begin(),end(),val);
}

template<class T>
Matrix<T> Matrix<T>::unary_old(T (*fcn)(const T&)) const{
  Matrix<T> temp;
  temp.unary_old(fcn,*this);
  return temp;
}

// template<class T>
// Matrix<T> Matrix<T>::unary_old(int op, const Matrix<T> &x){
//   // Return value
//   Matrix<T> ret(x.sparsity());
//   
//   // Do the operation on all non-zero elements
//   for(int el=0; el<size(); ++el)
//     casadi_math<T>::fun[op](x.data()[el],0,ret.data()[el]);
// 
//   // Check the value of the structural zero-entries, if there are any
//   if(!x.dense()){
//     // Get the value for the structural zeros
//     T fcn_0;
//     casadi_math<T>::fun[op](0,0,fcn_0);
//     if(!casadi_limits<T>::isZero(fcn_0)){
//       ret.makeDense(ret.size1(),ret.size2(),fcn_0);
//     }
//   }
//     
//   return ret;
// }

template<class T>
void Matrix<T>::unary_old(T (*fcn)(const T&), const Matrix<T>& x){
  // First check the value of the zero-entries
  T fcn_0 = fcn(0);
  if(casadi_limits<T>::isZero(fcn_0)){
    // Copy the matrix, including the sparsity pattern
    *this = x;
    
    // Do the operation on all non-zero elements
    for(int el=0; el<size(); ++el)
      data()[el] = fcn(x.data()[el]);
    
  } else {
    makeDense(x.size1(),x.size2(),fcn_0);
    for(int i=0; i<size1(); ++i){ // loop over rows
      for(int el=x.rowind(i); el<x.rowind(i+1); ++el){ // loop over non-zero elements
        int j = x.col(el);
        data()[j+i*size2()] = fcn(x.data()[el]);
      }
    }
  }
}

template<class T>
Matrix<T> Matrix<T>::binary_old(T (*fcn)(const T&, const T&), const Matrix<T>& y) const{
  Matrix<T> temp;
  temp.binary_old(fcn,*this,y);
  return temp;
}

template<class T>
void Matrix<T>::binary_old(T (*fcn)(const T&, const T&), const Matrix<T> &x, const Matrix<T> &y){
  if(x.scalar())
    if(y.scalar())
      *this = fcn(x.toScalar(),y.toScalar());
    else
      scalar_matrix_old(fcn,x.toScalar(),y);
  else if(y.scalar())
    matrix_scalar_old(fcn,x,y.toScalar());
  else
    matrix_matrix_old(fcn,x,y);
}

template<class T>
void Matrix<T>::scalar_matrix_old(T (*fcn)(const T&, const T&), const T& x, const Matrix<T>& y){
  T fcn_x_0 = fcn(x,0);
  if(casadi_limits<T>::isZero(fcn_x_0)){
    // Start with an empty matrix: all elements are added to the end!
    makeEmpty(y.size1(),y.size2());
    
    // Reserve space for the nonzeros
    reserve(y.size());
  } else {
    // Start with an dense matrix: element access in constant time
    makeDense(y.size1(),y.size2(),fcn_x_0);
  }
  
  for(int i=0; i<size1(); ++i){ // loop over rows
    for(int el=y.rowind(i); el<y.rowind(i+1); ++el){
      int j = y.col(el);
      elem(i,j) = fcn(x,y.at(el));
    }
  }
}

template<class T>
void Matrix<T>::matrix_scalar_old(T (*fcn)(const T&, const T&), const Matrix<T>& x, const T& y){
  T fcn_0_y = fcn(0,y);
  if(casadi_limits<T>::isZero(fcn_0_y)){
    // Start with an empty matrix: all elements are added to the end!
    makeEmpty(x.size1(),x.size2());
    
    // Reserve space for the nonzeros
    reserve(x.size());
    
  } else {
    // Start with an dense matrix: element access in constant time
    makeDense(x.size1(),x.size2(),fcn_0_y);
  }
  
  for(int i=0; i<size1(); ++i){ // loop over rows
    for(int el=x.rowind(i); el<x.rowind(i+1); ++el){
      int j = x.col(el);
      elem(i,j) = fcn(x.at(el),y);
    }
  }
}

template<class T>
void Matrix<T>::matrix_matrix_old(T (*fcn)(const T&, const T&), const Matrix<T>& x, const Matrix<T>& y){
  casadi_assert_message(x.size1() == y.size1() && x.size2() == y.size2(),
    "matrix_matrix: dimension mismatch." << std::endl << "Left argument has shape " << x.dimString() << ", right has shape " << y.dimString()
  ); 
  
  // Make a deep copy if *this and x or y is the same object
  if(this == &x)
    *this = x;
  else if(this == &y)
    *this = y;

  T fcn_0_0 = fcn(0,0);
  if(casadi_limits<T>::isZero(fcn_0_0)){
    // Start with an empty matrix: all elements are added to the end!
    makeEmpty(0,x.size2());
  
    // Reserve space for the nonzeros
    reserve(std::max(x.size(),y.size()),x.size1());
  
    for(int i=0; i<x.size1(); ++i){ // loop over rows
      // Resize matrix
      resize(i+1,x.size2());
      
      int el1 = x.rowind(i);
      int el2 = y.rowind(i);
      int k1 = x.rowind(i+1);
      int k2 = y.rowind(i+1);
      while(el1 < k1 || el2 < k2){
        int j1 = (el1 < k1) ? x.col(el1) : size2() ;
        int j2 = (el2 < k2) ? y.col(el2) : size2() ;
        
        if(j1==j2)
          elem(i,j1) = fcn(x.data()[el1++],y.data()[el2++]); 
        else if(j1>j2)
          elem(i,j2) = fcn(0,y.data()[el2++]);
        else
          elem(i,j1) = fcn(x.data()[el1++],0);
        }
    }
  
  } else {
    
    // Start with an dense matrix: element access in constant time
    makeDense(x.size1(),x.size2(),fcn_0_0);

    for(int i=0; i<size1(); ++i){ // loop over rows
    int el1 = x.rowind(i);
    int el2 = y.rowind(i);
    int k1 = x.rowind(i+1);
    int k2 = y.rowind(i+1);
    while(el1 < k1 || el2 < k2){
      int j1 = (el1 < k1) ? x.col(el1) : size2() ;
      int j2 = (el2 < k2) ? y.col(el2) : size2() ;
      
      if(j1==j2)
        elem(i,j1) = fcn(x.data()[el1++],y.data()[el2++]); 
      else if(j1>j2)
        elem(i,j2) = fcn(0,y.data()[el2++]);
      else
        elem(i,j1) = fcn(x.data()[el1++],0);
      }
    }
  }
}

template<class T>
Matrix<T> Matrix<T>::operator-() const{
  Matrix<T> temp;
  temp.unary_old(casadi_operators<T>::neg,*this);
  return temp;
}

template<class T>
Matrix<T> Matrix<T>::operator+() const{
  return *this;
}

template<class T>
Matrix<T> Matrix<T>::__add__(const Matrix<T> &y) const{
  Matrix<T> r;
  r.binary_old(casadi_operators<T>::add,*this,y);
  return r;
}

template<class T>
Matrix<T> Matrix<T>::__sub__(const Matrix<T> &y) const{
  Matrix<T> r;
  r.binary_old(casadi_operators<T>::sub,*this,y);
  return r;
}

template<class T>
Matrix<T> Matrix<T>::__mul__(const Matrix<T> &y) const{
  Matrix<T> r;
  r.binary_old(casadi_operators<T>::mul,*this,y);
  return r;
}

template<class T>
Matrix<T> Matrix<T>::__div__(const Matrix<T> &y) const{
  Matrix<T> r;
  r.binary_old(casadi_operators<T>::div,*this,y);
  return r;
}

// template<class T>
// Matrix<T>& Matrix<T>::operator+=(const double &y){
//   Matrix<T> x = *this;
//   binary_old(casadi_operators<T>::add,x,y);
//   return *this;
// }
// 
// template<class T>
// Matrix<T>& Matrix<T>::operator-=(const double &y){
//   Matrix<T> x = *this;
//   binary_old(casadi_operators<T>::sub,x,y);
//   return *this;
// }
// 
// template<class T>
// Matrix<T>& Matrix<T>::operator*=(const double &y){
//   Matrix<T> x = *this;
//   binary_old(casadi_operators<T>::mul,x,y);
//   return *this;
// }
// 
// template<class T>
// Matrix<T>& Matrix<T>::operator/=(const double &y){
//   Matrix<T> x = *this;
//   binary_old(casadi_operators<T>::div,x,y);
//   return *this;
// }
template<class T>
Matrix<T> Matrix<T>::__mrdivide__(const Matrix<T>& b) const { if (b.numel()==1) return *this/b; throw CasadiException("mrdivide: Not implemented");}

template<class T>
Matrix<T> Matrix<T>::__mldivide__(const Matrix<T>& b) const { if (b.numel()==1) return *this/b; throw CasadiException("mldivide: Not implemented");}

template<class T>
Matrix<T> Matrix<T>::__mpower__(const Matrix<T>& b) const { if (b.numel()==1) return (*this).__pow__(b); throw CasadiException("mpower: Not implemented");}

template<class T>
CRSSparsity& Matrix<T>::sparsityRef(){
  sparsity_.makeUnique();
  return sparsity_;
}

template<class T>
void Matrix<T>::getBand(int kl, int ku, int ldres, T *res) const{
  // delete the content of the matrix
  for(int j=0; j<size2(); ++j) // loop over columns
    for(int s=0; s<kl+ku+1; ++s) // loop over the subdiagonals
      res[s + ldres*j] = 0;
  
  // loop over rows
  for(int i=0; i<size1(); ++i){ 
    
    // loop over the non-zero elements
    for(int el=rowind(i); el<rowind(i+1); ++el){ 
      int j=col(el);  // column
      
      // Check if we have not yet inside the band
      if(j<i-kl) continue;

      // Check if we are already outside the band
      if(j>i+ku) break;

      // Get the subdiagonal
      int s = i - j + ku;

      // Store the element
      res[s + ldres*j] = data()[el];
    }
  }
}

template<class T>
void Matrix<T>::set(T val, Sparsity sp){
  setArray(&val,1,DENSE);
}
    
template<class T>
void Matrix<T>::get(T& val, Sparsity sp) const{
  getArray(&val,1,DENSE);
}

template<class T>
void Matrix<T>::set(const std::vector<T>& val, Sparsity sp){
  setArray(val.empty() ? 0 : &val.front(),val.size(),sp);
}

template<class T>
void Matrix<T>::get(std::vector<T>& val, Sparsity sp) const{
  getArray(val.empty() ? 0 : &val.front(),val.size(),sp);
}

template<class T>
void Matrix<T>::set(const Matrix<T>& val, Sparsity sp){
  if (!((size1()==0 || size2()==0) && (val.size1()==0 || val.size2()==0))) {
    casadi_assert_message(size1()==val.size1() && size2()==val.size2(),"Matrix<T>::set(Matrix<T>): shape mismatch. lhs is matrix of shape " << dimString() << ", while rhs is shape " << val.dimString() << ".");
  }
  //casadi_assert_message(val.sparsity()==sparsity(),"Matrix<T>::set(Imatrix): sparsity mismatch. Lhs is matrix of shape " << dimString() << ", while rhs is shape " << val.dimString() << ".");
  set(val.data(),sp);
}

template<class T>
void Matrix<T>::get(Matrix<T>& val, Sparsity sp) const{
  if (!((size1()==0 || size2()==0) && (val.size1()==0 || val.size2()==0))) {
    casadi_assert_message(size1()==val.size1() && size2()==val.size2(),"Matrix<T>::get(Matrix<T>): shape mismatch. lhs is matrix of shape " << dimString() << ", while rhs is shape " << val.dimString() << ".");
  }
  //casadi_assert_message(val.sparsity()==sparsity(),"Matrix<T>::get(Imatrix): sparsity mismatch. Lhs is matrix of shape " << dimString() << ", while rhs is shape " << val.dimString() << ".");
  get(val.data(),sp);
}

template<class T>
void Matrix<T>::set(const T* val, Sparsity sp){
  int len = sp==SPARSE ? size() : sp==DENSE ? numel() : sp==SPARSESYM ? sizeL() : -1;
  setArray(val,len,sp);
}

template<class T>
void Matrix<T>::get(T* val, Sparsity sp) const{
  int len = sp==SPARSE ? size() : sp==DENSE ? numel() : sp==SPARSESYM ? sizeL() : -1;
  getArray(val,len,sp);
}

template<class T>
void Matrix<T>::getArray(T* val, int len, Sparsity sp) const{
  const std::vector<T> &v = data();
  if(sp==SPARSE || (sp==DENSE && numel()==v.size())){
    casadi_assert_message(len==size(),
          "Matrix<T>::getArray: Dimension mismatch." << std::endl <<
          "I am a matrix of shape " << size1() << " x " << size2() << " = " << numel() << " ready to fill something up with " << size() << " non-zero elements," << std::endl <<
          "but you supplied something with " << len << " non-zero elements."); 
    copy(v.begin(),v.end(),val);
  } else if(sp==DENSE){
    casadi_assert_message(len==numel(),
        "Matrix<T>::getArray: Dimension mismatch." << std::endl <<
        "I am a dense matrix of shape " << size1() << " x " << size2() << " ready to fill something up with " << numel() << " elements," << std::endl <<
        "but you supplied something with " << len << " elements."); 
    int k=0; // index of the result
    for(int i=0; i<size1(); ++i) // loop over rows
      for(int el=rowind(i); el<rowind(i+1); ++el){ // loop over the non-zero elements
        int j=col(el);  // column
        for(; k<i*size2()+j; ++k)
          val[k] = 0; // add zeros before the non-zero element

      // add the non-zero element
      val[k] = v[el];
      k++;
    }
    // add sparse zeros at the end of the matrix
    for(; k<numel(); ++k)
     val[k] = 0;
  } else if(sp==SPARSESYM){
    // copy to the result vector
    int nz = 0;
    for(int row=0; row<size1(); ++row){
      // Loop over the elements in the row
      for(int el=rowind(row); el<rowind(row+1); ++el){ // loop over the non-zero elements
        if(col(el) > row) break; // break inner loop (only lower triangular part is used)
        val[nz++] = v[el];
      }
    }
  } else {
    casadi_error("Matrix<T>::getArray: not SPARSE or DENSE");
  }
}

/**
Set stride to zero for unstrided acces
*/
template<class T>
void Matrix<T>::getStridedArray(T* val, int len, int stride1, int stride2, Sparsity sp) const{
  if (stride1==0 || stride2==0 || (stride2==1 && stride1==size2())) 
    return getArray(val, len, sp);
    
  const std::vector<T> &v = data();
  if(sp==SPARSE){
    throw CasadiException("Matrix<T>::getArray: strided SPARSE not implemented");
  } else if(sp==DENSE && numel()==v.size()) {
    for (int i=0;i<size1();i++) {
      for (int j=0;j<size2();j++) {
        val[i*stride1+j*stride2] = v[i*size2()+j];
      }
    }
  } else if(sp==DENSE){
    throw CasadiException("Matrix<T>::getArray: strided sparse DMatrix->dense not implemented");
  } else if(sp==SPARSESYM){
    throw CasadiException("Matrix<T>::getArray: strided SPARSESYM not implemented");
  } else {
    throw CasadiException("Matrix<T>::getArray: not SPARSE or DENSE");
  }

}

template<class T>
void Matrix<T>::setArray(const T* val, int len, Sparsity sp){
  std::vector<T> &v = data();
  if(sp==SPARSE || (sp==DENSE && numel()==size())){
    casadi_assert_message(len==size(),
          "Matrix<T>::setArray: Dimension mismatch." << std::endl <<
          "I am a matrix of shape " << size1() << " x " << size2() << " = " << numel() << " ready to be filled up with " << size() << " non-zero elements," << std::endl <<
          "but you supplied something with " << len << " non-zero elements."); 
    copy(val,val+len,v.begin());
  } else if(sp==DENSE){
    casadi_assert_message(len==numel(),
          "Matrix<T>::setArray: Dimension mismatch." << std::endl <<
          "I am a dense matrix of shape " << size1() << " x " << size2() << " ready to be filled up with " << numel() << " elements," << std::endl <<
          "but you supplied something with " << len << " elements."); 
    for(int i=0; i<size1(); ++i) // loop over rows
      for(int el=rowind(i); el<rowind(i+1); ++el){ // loop over the non-zero elements
        // column
        int j=col(el);
        
        // Set the element
        v[el] = val[i*size2()+j];
    }
  } else {
    throw CasadiException("Matrix<T>::setArray: not SPARSE or DENSE");
  }
}

template<class T>
void Matrix<T>::getArray(T* val) const{
  getArray(val,size(),SPARSE);
}

template<class T>
void Matrix<T>::setArray(const T* val){
  setArray(val,size(),SPARSE);
}






template<class T>
int Matrix<T>::size(Sparsity sp) const{
  if(sp==SPARSE){
    return size();
  } else if(sp==SPARSESYM){
    return sizeL();
  } else if(sp==DENSE){
    return numel();
  } else if(sp==DENSESYM){
    return (numel()+size1())/2;
  } else {
      throw CasadiException("Matrix<T>::size(Sparsity): unknown sparsity");
  }
}

template<class T>
Matrix<T> Matrix<T>::__pow__(const Matrix<T>& y) const{
  return binary_old(CasADi::casadi_operators<T>::pow,y);
}

template<class T>
Matrix<T> Matrix<T>::__constpow__(const Matrix<T>& y) const{
  return binary_old(CasADi::casadi_operators<T>::constpow,y);
}

template<class T>
Matrix<T> Matrix<T>::sin() const{
  return unary_old(CasADi::casadi_operators<T>::sin);
}

template<class T>
Matrix<T> Matrix<T>::cos() const{
  return unary_old(CasADi::casadi_operators<T>::cos);
}

template<class T>
Matrix<T> Matrix<T>::tan() const{
  return unary_old(CasADi::casadi_operators<T>::tan);
}

template<class T>
Matrix<T> Matrix<T>::erf() const{
  return unary_old(CasADi::casadi_operators<T>::erf);
}

template<class T>
Matrix<T> Matrix<T>::arcsin() const{
  return unary_old(CasADi::casadi_operators<T>::asin);
}

template<class T>
Matrix<T> Matrix<T>::arccos() const{
  return unary_old(CasADi::casadi_operators<T>::acos);
}

template<class T>
Matrix<T> Matrix<T>::arctan() const{
  return unary_old(CasADi::casadi_operators<T>::atan);
}

template<class T>
Matrix<T> Matrix<T>::sinh() const{
  return unary_old(CasADi::casadi_operators<T>::sinh);
}

template<class T>
Matrix<T> Matrix<T>::cosh() const{
  return unary_old(CasADi::casadi_operators<T>::cosh);
}

template<class T>
Matrix<T> Matrix<T>::tanh() const{
  return unary_old(CasADi::casadi_operators<T>::tanh);
}

template<class T>
Matrix<T> Matrix<T>::exp() const{
  return unary_old(CasADi::casadi_operators<T>::exp);
}

template<class T>
Matrix<T> Matrix<T>::log() const{
  return unary_old(CasADi::casadi_operators<T>::log);
}

template<class T>
Matrix<T> Matrix<T>::log10() const{
  return log()*(1/std::log(10.));
}

template<class T>
Matrix<T> Matrix<T>::sqrt() const{
  return unary_old(CasADi::casadi_operators<T>::sqrt);
}

template<class T>
Matrix<T> Matrix<T>::floor() const{
  return unary_old(CasADi::casadi_operators<T>::floor);
}

template<class T>
Matrix<T> Matrix<T>::ceil() const{
  return unary_old(CasADi::casadi_operators<T>::ceil);
}

template<class T>
Matrix<T> Matrix<T>::fabs() const{
  return unary_old(CasADi::casadi_operators<T>::fabs);
}

template<class T>
Matrix<T> Matrix<T>::sign() const{
  return unary_old(CasADi::casadi_operators<T>::sign);
}

template<class T>
Matrix<T> Matrix<T>::erfinv() const{
  return unary_old(CasADi::casadi_operators<T>::erfinv);
}

template<class T>
Matrix<T> Matrix<T>::fmin(const Matrix<T>& y) const{
  return binary_old(CasADi::casadi_operators<T>::fmin, y);
}

template<class T>
Matrix<T> Matrix<T>::fmax(const Matrix<T>& y) const{
  return binary_old(CasADi::casadi_operators<T>::fmax, y);
}

template<class T>
Matrix<T> Matrix<T>::printme(const Matrix<T>& y) const{
  return binary_old(CasADi::casadi_operators<T>::printme, y);
}

template<class T>
std::vector<T>& Matrix<T>::data(){
  return data_;  
}
    
template<class T>
const std::vector<T>& Matrix<T>::data() const{
  return data_;  
}

template<class T>
void Matrix<T>::erase(const std::vector<int>& ii, const std::vector<int>& jj){
  // Erase from sparsity pattern
  std::vector<int> mapping = sparsityRef().erase(ii,jj);
  
  // Update non-zero entries
  for(int k=0; k<mapping.size(); ++k)
    data()[k] = data()[mapping[k]];
    
  // Truncate nonzero vector
  data().resize(mapping.size());
}

template<class T>
void Matrix<T>::enlarge(int nrow, int ncol, const std::vector<int>& ii, const std::vector<int>& jj){
  sparsityRef().enlarge(nrow,ncol,ii,jj);
}

template<class T>
std::string Matrix<T>::dimString() const {
  std::stringstream ss;
  ss << "(" << size1() << " x " << size2();
  if (size1()!=1 && size2()!=1)  ss << " = " << numel();
  if (numel()==size()) {
    ss << " dense";
  } else {
    ss << " | " << size();
  }
  ss << ")";
  return ss.str();
}

template<class T>
void Matrix<T>::sanityCheck(bool complete) const {
  sparsity_.sanityCheck(complete);
  
  if (data_.size()!=sparsity_.col().size()) {
      std::stringstream s;
      s << "Matrix:Compressed Row Storage is not sane. The following must hold:" << std::endl;
      s << "  data.size() = ncol, but got   col.size()  = " << data_.size() << "   and   ncol = "  << sparsity_.col().size() << std::endl;
      s << "  Note that the signature is as follows: DMatrix (nrow, ncol, col, rowind, data)." << std::endl;
      casadi_error(s.str());
  }
}

template<class T>
Matrix<T> Matrix<T>::mul(const Matrix<T> &y) const{
  // First factor
  const Matrix<T>& x = *this;
  
  // Return object (assure RVO)
  Matrix<T> ret;
  
  if (x.numel()==1 || y.numel()==1){
    // Elementwise multiplication when either x or y is a scalar
    ret = x*y;
    
  } else {
    // Matrix multiplication
    
    casadi_assert_message(x.size2() == y.size1(),
      "Matrix<T>::mul: dimension mismatch. Attemping product of (" << x.size1() << " x " << x.size2() << ") " << std::endl <<
      "with (" << y.size1() << " x " << y.size2() << ") matrix."
    );

    // Form the transpose of y
    Matrix<T> y_trans = y.trans();
  
    // Create the sparsity pattern for the matrix-matrix product
    CRSSparsity spres = x.sparsity().patternProduct(y_trans.sparsity());
  
    // Create the return object
    ret = Matrix<T>(spres, 0);
  
    // Carry out the matrix product
    mul_no_alloc(x,y_trans,ret);
  }
  
  return ret;
}

template<class T>
void Matrix<T>::mul_no_alloc1(Matrix<T> &x, const Matrix<T> &y_trans, const Matrix<T>& z){

  // Direct access to the arrays
  const std::vector<T> &z_data = z.data();
  const std::vector<int> &z_col = z.col();
  const std::vector<int> &z_rowind = z.rowind();
  std::vector<T> &x_data = x.data();
  const std::vector<int> &x_col = x.col();
  const std::vector<T> &y_trans_data = y_trans.data();
  const std::vector<int> &y_row = y_trans.col();
  const std::vector<int> &x_rowind = x.rowind();
  const std::vector<int> &y_colind = y_trans.rowind();

  // loop over the row of the resulting matrix)
  for(int i=0; i<z.size1(); ++i){
    for(int el=z_rowind[i]; el<z_rowind[i+1]; ++el){ // loop over the non-zeros of the resulting matrix
      int j = z_col[el];
      int el1 = x_rowind[i];
      int el2 = y_colind[j];
      while(el1 < x_rowind[i+1] && el2 < y_colind[j+1]){ // loop over non-zero elements
        int j1 = x_col[el1];
        int i2 = y_row[el2];      
        if(j1==i2){
          x_data[el1++] += z_data[el]*y_trans_data[el2++];
        } else if(j1<i2) {
          el1++;
        } else {
          el2++;
        }
      }
    }
  }
}


template<class T>
void Matrix<T>::mul_no_alloc2(const Matrix<T> &x, Matrix<T> &y_trans, const Matrix<T>& z){

  // Direct access to the arrays
  const std::vector<T> &z_data = z.data();
  const std::vector<int> &z_col = z.col();
  const std::vector<int> &z_rowind = z.rowind();
  const std::vector<T> &x_data = x.data();
  const std::vector<int> &x_col = x.col();
  std::vector<T> &y_trans_data = y_trans.data();
  const std::vector<int> &y_row = y_trans.col();
  const std::vector<int> &x_rowind = x.rowind();
  const std::vector<int> &y_colind = y_trans.rowind();

  // loop over the row of the resulting matrix)
  for(int i=0; i<z.size1(); ++i){
    for(int el=z_rowind[i]; el<z_rowind[i+1]; ++el){ // loop over the non-zeros of the resulting matrix
      int j = z_col[el];
      int el1 = x_rowind[i];
      int el2 = y_colind[j];
      while(el1 < x_rowind[i+1] && el2 < y_colind[j+1]){ // loop over non-zero elements
        int j1 = x_col[el1];
        int i2 = y_row[el2];      
        if(j1==i2){
          y_trans_data[el2++] += z_data[el] * x_data[el1++];
        } else if(j1<i2) {
          el1++;
        } else {
          el2++;
        }
      }
    }
  }
}

template<class T>
void Matrix<T>::mul_no_alloc(const Matrix<T> &x, const Matrix<T> &y_trans, Matrix<T>& z){
  // Direct access to the arrays
  const std::vector<int> &z_col = z.col();
  const std::vector<int> &z_rowind = z.rowind();
  const std::vector<int> &x_col = x.col();
  const std::vector<int> &y_row = y_trans.col();
  const std::vector<int> &x_rowind = x.rowind();
  const std::vector<int> &y_colind = y_trans.rowind();

  const std::vector<T> &x_data = x.data();
  const std::vector<T> &y_trans_data = y_trans.data();
  std::vector<T> &z_data = z.data();

  // loop over the rows of the resulting matrix)
  for(int i=0; i<z_rowind.size()-1; ++i){
    for(int el=z_rowind[i]; el<z_rowind[i+1]; ++el){ // loop over the non-zeros of the resulting matrix
      int j = z_col[el];
      int el1 = x_rowind[i];
      int el2 = y_colind[j];
      while(el1 < x_rowind[i+1] && el2 < y_colind[j+1]){ // loop over non-zero elements
        int j1 = x_col[el1];
        int i2 = y_row[el2];      
        if(j1==i2){
          z_data[el] += x_data[el1++] * y_trans_data[el2++];
        } else if(j1<i2) {
          el1++;
        } else {
          el2++;
        }
      }
    }
  }
}

template<class T>
void Matrix<T>::mul_sparsity(Matrix<T> &x, Matrix<T> &y_trans, Matrix<T>& z, bool fwd){
  // Direct access to the arrays
  const std::vector<int> &z_col = z.col();
  const std::vector<int> &z_rowind = z.rowind();
  const std::vector<int> &x_col = x.col();
  const std::vector<int> &y_row = y_trans.col();
  const std::vector<int> &x_rowind = x.rowind();
  const std::vector<int> &y_colind = y_trans.rowind();

  // Convert data array to arrays of integers
  bvec_t *x_data = get_bvec_t(x.data());
  bvec_t *y_trans_data = get_bvec_t(y_trans.data());
  bvec_t *z_data = get_bvec_t(z.data());
  
  // loop over the rows of the resulting matrix)
  for(int i=0; i<z_rowind.size()-1; ++i){
    for(int el=z_rowind[i]; el<z_rowind[i+1]; ++el){ // loop over the non-zeros of the resulting matrix
      int j = z_col[el];
      int el1 = x_rowind[i];
      int el2 = y_colind[j];
      if(fwd) z_data[el] = 0;
      while(el1 < x_rowind[i+1] && el2 < y_colind[j+1]){ // loop over non-zero elements
        int j1 = x_col[el1];
        int i2 = y_row[el2];      
        if(j1==i2){
          // | and not & since we are propagating dependencies
          if(fwd){
            z_data[el] |= x_data[el1] | y_trans_data[el2];
          } else {
            x_data[el1] |= z_data[el];
            y_trans_data[el2] |= z_data[el];
          }
          el1++;
          el2++;
        } else if(j1<i2) {
          el1++;
        } else {
          el2++;
        }
      }
    }
  }
}

template<class T>
Matrix<T> Matrix<T>::trans() const{
  // quick return if empty or scalar
  if(empty() || scalar()) return *this;

  // Create the new sparsity pattern and the mapping
  std::vector<int> mapping;
  CRSSparsity s = sparsity().transpose(mapping);

  // create the return matrix
  Matrix<T> ret(s);
  
  // Copy the content
  for(int i=0; i<mapping.size(); ++i)
    ret[i] = data()[mapping[i]];
  
  return ret;
}

// template<class T>
// Matrix<T>::operator const T() const{
//   return toScalar();
// }

template<class T>
const T Matrix<T>::toScalar() const{
  // Make sure that the matrix is 1-by-1
  casadi_assert_message(scalar(),"Can only convert 1-by-1 matrices to scalars");

  // return zero or the nonzero element
  if(size()==1)
    return data()[0];
  else
    return casadi_limits<T>::zero;
}

template<class T>
Matrix<T> Matrix<T>::binary(int op, const Matrix<T> &x, const Matrix<T> &y){
  if(x.numel()==1)
    return scalar_matrix(op,x,y);
  else if(y.numel()==1)  
    return matrix_scalar(op,x,y);
  else
    return matrix_matrix(op,x,y);
}

template<class T>
void Matrix<T>::binary_no_alloc(void (*fcn)(unsigned char op, const T&, const T&, T&), unsigned char op, const Matrix<T> &x, const Matrix<T> &y, Matrix<T>& r, const std::vector<unsigned char>& mapping){
  std::vector<T>& rd = r.data();
  const std::vector<T> &xd = x.data();
  const std::vector<T> &yd = y.data();

  // Argument values
  T zero = 0;

  // Nonzero counters
  int el0=0, el1=0, el=0;
  
  // Loop over nonzero elements
  for(int i=0; i<mapping.size(); ++i){
    // Check which elements are nonzero
    unsigned char m = mapping[i];
    bool nz0(m & 1);
    bool nz1(m & 2);
    bool skip_nz(m & 4);
    
    // Evaluate
    if(!skip_nz){
      fcn(op, nz0 ? xd[el0] : zero, nz1 ? yd[el1] : zero, rd[el++]);
    }
    
    // Go to next nonzero
    el0 += nz0;
    el1 += nz1;
  }
}

template<class T>
Matrix<T> Matrix<T>::unary(int op, const Matrix<T> &x){
  casadi_assert_message(0,"not implemented");
  Matrix<T> ret;
  return ret;
}

template<class T>
Matrix<T> Matrix<T>::scalar_matrix(int op, const Matrix<T> &x, const Matrix<T> &y){
  casadi_assert_message(0,"not implemented");
  Matrix<T> ret;
  return ret;
}

template<class T>
Matrix<T> Matrix<T>::matrix_scalar(int op, const Matrix<T> &x, const Matrix<T> &y){
  casadi_assert_message(0,"not implemented");
  Matrix<T> ret;
  return ret;
}

template<class T>
Matrix<T> Matrix<T>::matrix_matrix(int op, const Matrix<T> &x, const Matrix<T> &y){
  casadi_assert_message(0,"not implemented");
  Matrix<T> ret;
  return ret;
}

template<class T>
Matrix<T> Matrix<T>::sparse(const std::pair<int,int> &nm){
  return sparse(nm.first,nm.second);
}

template<class T>
Matrix<T> Matrix<T>::sparse(int n, int m){
  return Matrix<T>(n,m);
}

template<class T>
Matrix<T> Matrix<T>::zeros(const std::pair<int,int> &nm){
  return zeros(nm.first,nm.second);
}

template<class T>
Matrix<T> Matrix<T>::zeros(int n, int m){
  return Matrix<T>(n,m,0);
}

template<class T>
Matrix<T> Matrix<T>::ones(const std::pair<int,int> &nm){
  return ones(nm.first,nm.second);
}

template<class T>
Matrix<T> Matrix<T>::ones(int n, int m){
  return Matrix<T>(n,m,1);
}

template<class T>
Matrix<T> Matrix<T>::repmat(const Matrix<T>& x, const std::pair<int,int>& nm){
  return repmat(x,nm.first,nm.second);
}

template<class T>
Matrix<T> Matrix<T>::repmat(const Matrix<T>& x, int nrow, int ncol){
  if(x.scalar()){
    if(x.dense()){
      return Matrix<T>(nrow,ncol,x.toScalar());
    } else {
      return sparse(nrow,ncol);
    }
  } else {
    return horzcat(std::vector< Matrix<T> >(ncol,vertcat(std::vector< Matrix<T> >(nrow,x))));
  }
}

template<class T>
Matrix<T> Matrix<T>::eye(int n){
  return Matrix<T>(CRSSparsity::createDiagonal(n),1);
}

template<class T>
Matrix<T> Matrix<T>::inf(const std::pair<int,int>& nm){
  return inf(nm.first, nm.second);
}

template<class T>
Matrix<T> Matrix<T>::inf(int n, int m){
  casadi_assert_message(std::numeric_limits<T>::has_infinity,"Datatype cannot represent infinity");
  return Matrix<T>(n,m,std::numeric_limits<T>::infinity());
}

template<class T>
Matrix<T> Matrix<T>::nan(const std::pair<int,int>& nm){
  return nan(nm.first, nm.second);
}

template<class T>
Matrix<T> Matrix<T>::nan(int n, int m){
  casadi_assert_message(std::numeric_limits<T>::has_quiet_NaN,"Datatype cannot represent not-a-number");
  return Matrix<T>(n,m,std::numeric_limits<T>::quiet_NaN());
}

template<class T>
void Matrix<T>::append(const Matrix<T>& y){
  // Quick return if we are adding an empty expression
  if(y.empty()) return;

  // Likewise if expr is empty
  if(empty()){
    *this=y;
    return;
  }
  
  // Append the sparsity pattern
  sparsityRef().append(y.sparsity());
  
  // Add the non-zeros
  data().insert(end(),y.begin(),y.end());
}

} // namespace CasADi

#endif // MATRIX_IMPL_HPP

