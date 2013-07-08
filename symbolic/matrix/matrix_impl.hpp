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
#include "matrix_tools.hpp"
#include "sparsity_tools.hpp"

#ifdef WITH_EIGEN3
#include <Eigen/Dense>
#endif

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
int Matrix<T>::stream_precision_ = 6;
template<class T>
int Matrix<T>::stream_width_ = 0;
template<class T>
bool Matrix<T>::stream_scientific_ = false;

template<class T>
T& Matrix<T>::elem(int i, int j){
  int oldsize = sparsity().size();
  int ind = sparsityRef().getNZ(i,j);
  if(oldsize != sparsity().size())
    data().insert(begin()+ind,T(0));
  return at(ind);
}

template<class T>
bool Matrix<T>::__nonzero__() const {
  if (isNull()) {casadi_error("Cannot determine truth value of null Matrix.");}
  if (numel()!=1) {casadi_error("Only scalar Matrix could have a truth value, but you provided a shape" << dimString());}
  return CasADi::__nonzero__(at(0));
}

template<class T>
const Matrix<T> Matrix<T>::sub(int i, int j) const{
  return elem(i,j);
}

template<class T>
const Matrix<T> Matrix<T>::sub(const std::vector<int>& ii, const std::vector<int>& jj) const{
  // Nonzero mapping from submatrix to full
  std::vector<int> mapping;
  
  // Get the sparsity pattern - does bounds checking
  CRSSparsity sp = sparsity().sub(ii,jj,mapping);

  // Create return object
  Matrix<T> ret(sp);
  
  // Copy nonzeros
  for(int k=0; k<mapping.size(); ++k)
    ret.data()[k] = data()[mapping[k]];
  
  // Return (RVO)
  return ret;
}

template<class T>
const Matrix<T> Matrix<T>::sub(const std::vector<int>& ii, const Matrix<int>& k) const{
  std::vector< int > cols = range(size2());
  std::vector< Matrix<T> > temp;
  
  if (!inBounds(ii,size1())) {
    casadi_error("Slicing [ii,k] out of bounds. Your ii contains " << *std::min_element(ii.begin(),ii.end()) << " up to " << *std::max_element(ii.begin(),ii.end()) << ", which is outside of the matrix shape " << dimString() << ".");
  }

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
const Matrix<T> Matrix<T>::sub(const Matrix<int>& k, const std::vector<int>& jj) const{
  std::vector< int > rows = range(size1());
  std::vector< Matrix<T> > temp;

  if (!inBounds(jj,size2())) {
    casadi_error("Slicing [ii,k] out of bounds. Your jj contains " << *std::min_element(jj.begin(),jj.end()) << " up to " << *std::max_element(jj.begin(),jj.end()) << ", which is outside of the matrix shape " << dimString() << ".");
  }
  
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
const Matrix<T> Matrix<T>::sub(const Matrix<int>& i, const Matrix<int>& j) const {
   casadi_assert_message(i.sparsity()==j.sparsity(),"sub(Imatrix i, Imatrix j): sparsities must match. Got " << i.dimString() << " and " << j.dimString() << ".");

   Matrix<T> ret(i.sparsity());
   for (int k=0;k<i.size();++k) {
    ret.data()[k] = elem(i.at(k),j.at(k));
   }

   return ret;
}

template<class T>
const Matrix<T> Matrix<T>::sub(const CRSSparsity& sp, int dummy) const {
  casadi_assert_message(size1()==sp.size1() && size2()==sp.size2(),"sub(CRSSparsity sp): shape mismatch. This matrix has shape " << size1() << " x " << size2() << ", but supplied sparsity index has shape " << sp.size1() << " x " << sp.size2() << "." );
  Matrix<T> ret(sp);

  std::vector<unsigned char> mapping; // Mapping that will be filled by patternunion
  sparsity().patternCombine(sp, false, true, mapping);

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
void Matrix<T>::setSub(const Matrix<T>& m, int i, int j){
  if(m.dense()){
    elem(i,j) = m.toScalar();
  } else {
    setSub(m,std::vector<int>(1,i),std::vector<int>(1,j));
  }
}

template<class T>
void Matrix<T>::setSub(const Matrix<T>& m, const std::vector<int>& ii, const std::vector<int>& jj){
  casadi_assert_message(m.numel()==1 || (ii.size() == m.size1() && jj.size() == m.size2()),"Dimension mismatch." << std::endl << "lhs is " << ii.size() << " x " << jj.size() << ", while rhs is " << m.dimString());

  if (!inBounds(ii,size1())) {
    casadi_error("setSub [.,ii,jj] out of bounds. Your ii contains " << *std::min_element(ii.begin(),ii.end()) << " up to " << *std::max_element(ii.begin(),ii.end()) << ", which is outside of the matrix shape " << dimString() << ".");
  }
  if (!inBounds(jj,size2())) {
    casadi_error("setSub[.,ii,jj] out of bounds. Your jj contains " << *std::min_element(jj.begin(),jj.end()) << " up to " << *std::max_element(jj.begin(),jj.end()) << ", which is outside of the matrix shape " << dimString() << ".");
  }
  
  // If m is scalar
  if(m.numel() != ii.size() * jj.size()){
    setSub(Matrix<T>(ii.size(),jj.size(),m.toScalar()),ii,jj);
    return;
  }

  if(dense() && m.dense()){
    // Dense mode
    for(int i=0; i<ii.size(); ++i) {
      for(int j=0; j<jj.size(); ++j) {
        data()[ii[i]*size2() + jj[j]]=m.data()[i*m.size2()+j];
      }
    }
  } else {
    // Sparse mode

    // Remove submatrix to be replaced
    erase(ii,jj);

    // Extend el to the same dimension as this
    Matrix<T> el_ext = m;
    el_ext.enlarge(size1(),size2(),ii,jj);

    // Unite the sparsity patterns
    *this = unite(*this,el_ext);
  }
}

template<class T>
void Matrix<T>::setSub(const Matrix<T>& m, const Matrix<int>& i, const std::vector<int>& jj) {
  // If el is scalar
  if(m.scalar() && (jj.size() > 1 || i.size() > 1)){
    setSub(repmat(Matrix<T>(i.sparsity(),m.toScalar()),1,jj.size()),i,jj);
    return;
  }

  if (!inBounds(jj,size2())) {
    casadi_error("setSub[.,i,jj] out of bounds. Your jj contains " << *std::min_element(jj.begin(),jj.end()) << " up to " << *std::max_element(jj.begin(),jj.end()) << ", which is outside of the matrix shape " << dimString() << ".");
  }
  
  CRSSparsity result_sparsity = repmat(i,1,jj.size()).sparsity();
  
  
  casadi_assert_message(result_sparsity == m.sparsity(),"setSub(Imatrix" << i.dimString() << ",Ivector(length=" << jj.size() << "),Matrix<T>)::Dimension mismatch. The sparsity of repmat(Imatrix,1," << jj.size() << ") = " << result_sparsity.dimString()  << " must match the sparsity of Matrix<T> = "  << m.dimString() << ".");
  
  
  std::vector<int> slice_i = range(i.size1());
  
  for(int k=0; k<jj.size(); ++k) {
     Matrix<T> el_k = m(slice_i,range(k*i.size2(),(k+1)*i.size2()));
     for (int j=0;j<i.size();++j) {
       elem(i.at(j),jj[k])=el_k.at(j);
     }
  }
  
}

template<class T>
void Matrix<T>::setSub(const Matrix<T>& m, const std::vector<int>& ii, const Matrix<int>& j) {
  
  // If el is scalar
  if(m.scalar() && (ii.size() > 1 || j.size() > 1)){
    setSub(repmat(Matrix<T>(j.sparsity(),m.toScalar()),ii.size(),1),ii,j);
    return;
  }

  if (!inBounds(ii,size1())) {
    casadi_error("setSub[.,ii,j] out of bounds. Your ii contains " << *std::min_element(ii.begin(),ii.end()) << " up to " << *std::max_element(ii.begin(),ii.end()) << ", which is outside of the matrix shape " << dimString() << ".");
  }
  
  CRSSparsity result_sparsity = repmat(j,ii.size(),1).sparsity();
  
  
  casadi_assert_message(result_sparsity == m.sparsity(),"setSub(Ivector(length=" << ii.size() << "),Imatrix" << j.dimString() << ",Matrix<T>)::Dimension mismatch. The sparsity of repmat(Imatrix," << ii.size() << ",1) = " << result_sparsity.dimString() << " must match the sparsity of Matrix<T> = " << m.dimString() << ".");
  
  std::vector<int> slice_j = range(j.size2());
  
  for(int k=0; k<ii.size(); ++k) {
     Matrix<T> el_k = m(range(k*j.size1(),(k+1)*j.size1()),slice_j);
     for (int i=0;i<j.size();++i) {
       elem(ii[k],j.at(i))=el_k.at(i);
     }
  }
  
}


template<class T>
void Matrix<T>::setSub(const Matrix<T>& m, const Matrix<int>& i, const Matrix<int>& j) {
   casadi_assert_message(i.sparsity()==j.sparsity(),"setSub(., Imatrix i, Imatrix j): sparsities must match. Got " << i.dimString() << " for i and " << j.dimString() << " for j.");

  // If m is scalar
  if(m.scalar() && i.numel() > 1){
    setSub(Matrix<T>(i.sparsity(),m.toScalar()),i,j);
    return;
  }
  
  casadi_assert_message(m.sparsity()==i.sparsity(),"setSub(Matrix m, Imatrix i, Imatrix j): sparsities must match. Got " << m.dimString() << " for m and " << j.dimString() << " for i and j.");
  
  for(int k=0; k<i.size(); ++k) {
     elem(i.at(k),j.at(k)) = m.at(k); 
  }
}

template<class T>
void Matrix<T>::setSub(const Matrix<T>& m, const CRSSparsity& sp, int dummy) {
  casadi_assert_message(size1()==sp.size1() && size2()==sp.size2(),"sub(CRSSparsity sp): shape mismatch. This matrix has shape " << size1() << " x " << size2() << ", but supplied sparsity index has shape " << sp.size1() << " x " << sp.size2() << "." );
  // TODO: optimize this for speed
  Matrix<T> elm;
  if (m.scalar()) {
    elm = Matrix<T>(sp,m.at(0));
  } else {
    elm = m.sub(sp);
  }

  for(int i=0; i<sp.rowind().size()-1; ++i){
    for(int k=sp.rowind()[i]; k<sp.rowind()[i+1]; ++k){
      int j=sp.col()[k];
      elem(i,j)=elm.data()[k];
    }
  }
}

template<class T>
const Matrix<T> Matrix<T>::getNZ(const std::vector<int>& k) const{
  try{
    Matrix<T> ret(k.size(),1,0);
    for(int el=0; el<k.size(); ++el)
      ret.data()[el] = data().at(k[el]);
  
    return ret;
  } catch(std::out_of_range& ex){
    std::stringstream ss;
    ss << "Out of range error in Matrix<>::getNZ: " << k << " not all in range [0," << size() << ")";
    throw CasadiException(ss.str());
  }
}

template<class T>
const Matrix<T> Matrix<T>::getNZ(const Matrix<int>& k) const{
  try{
    Matrix<T> ret(k.sparsity(),0);
    for(int el=0; el<k.size(); ++el)
      ret.data()[el] = data().at(k.at(el));
  
    return ret;
  } catch(std::out_of_range& ex){
    std::stringstream ss;
    ss << "Out of range error in Matrix<>::getNZ: " << k << " not all in range [0," << size() << ")";
    throw CasadiException(ss.str());
  }
}

template<class T>
void Matrix<T>::setNZ(int k, const Matrix<T>& m){
  if (k<0) k+=size();
  at(k) = m.toScalar();
}

template<class T>
void Matrix<T>::setNZ(const std::vector<int>& kk, const Matrix<T>& m){
  if (m.scalar()){
    // Assign all elements with the same scalar
    for(int k=0; k<kk.size(); ++k){
      setNZ(kk[k],m);
    }
  } else {
    // Assignment elementwise
    casadi_assert_message(kk.size()==m.size(),"Matrix<T>::setNZ: length of non-zero indices (" << kk.size() << ") " << std::endl << "must match size of rhs (" << m.size() << ").");
    for(int k=0; k<kk.size(); ++k){
      setNZ(kk[k],m[k]);
    }
  }
}

template<class T>
void Matrix<T>::setNZ(const Matrix<int>& kk, const Matrix<T>& m){
  if (m.scalar()){
    // Assign all elements with the same scalar
    for(int k=0; k<kk.size(); ++k){
      setNZ(kk.at(k),m);
    }
  } else if (kk.dense() && !m.dense() && kk.size1()==m.size1() && kk.size2()==m.size2()) {
    const std::vector<int>& col = m.sparsity().col();
    const std::vector<int>& rowind = m.sparsity().rowind();
    for(int i=0; i<rowind.size()-1; ++i){
      for(int k=rowind[i]; k<rowind[i+1]; ++k){
        int j=col[k];
        setNZ(kk.elem(i,j),m[k]);
      }
    }
  } else {
    casadi_assert_message(kk.sparsity()==m.sparsity(),"Matrix<T>::setNZ: sparsity of IMatrix index " << kk.dimString() << " " << std::endl << "must match sparsity of rhs " << m.dimString() << ".");
    for(int k=0; k<kk.size(); ++k){
      setNZ(kk.at(k),m[k]);
    }
  }
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
bool Matrix<T>::vector() const{
  return size2()==1;
}

template<class T>
Matrix<T>::Matrix() : sparsity_(CRSSparsity(0,0,false)){
}

template<class T>
Matrix<T>::Matrix(const Matrix<T>& m) : sparsity_(m.sparsity_), data_(m.data_){
}

template<class T>
Matrix<T>::Matrix(const std::vector<T>& x) : sparsity_(CRSSparsity(x.size(),1,true)), data_(x){
}

template<class T>
Matrix<T>::Matrix(const std::vector<T>& x, int n, int m) : sparsity_(CRSSparsity(n,m,true)), data_(x){
  casadi_assert_message(x.size() == n*m, "Dimension mismatch." << std::endl << "You supplied a vector of length " << x.size() << ", but " << n << " x " << m << " = " << n*m);
}

template<class T>
Matrix<T>& Matrix<T>::operator=(const Matrix<T>& m){
  sparsity_ = m.sparsity_;
  data_ = m.data_;
  return *this;
}

template<class T>
Matrix<T>::Matrix(int n, int m) : sparsity_(CRSSparsity(n,m,false)){
}

template<class T>
Matrix<T>::Matrix(int n, int m, const T& val) : sparsity_(CRSSparsity(n,m,true)), data_(std::vector<T>(n*m, val)){
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
  
  std::streamsize precision = stream.precision();
  std::streamsize width = stream.width();
  std::ios_base::fmtflags flags = stream.flags();
  
  stream.precision(stream_precision_);
  stream.width(stream_width_);
  if (stream_scientific_) {
    stream.setf(std::ios::scientific);
  } else {
    stream.unsetf(std::ios::scientific);
  }
  
  stream << toScalar();
  
  stream.precision(precision);
  stream.width(width);
  stream.flags(flags); 
}
  
template<class T>
void Matrix<T>::printVector(std::ostream &stream) const {
  casadi_assert_message(vector(),"Not a vector");
  
  std::streamsize precision = stream.precision();
  std::streamsize width = stream.width();
  std::ios_base::fmtflags flags = stream.flags();
  
  stream.precision(stream_precision_);
  stream.width(stream_width_);
  if (stream_scientific_) {
    stream.setf(std::ios::scientific);
  } else {
    stream.unsetf(std::ios::scientific);
  }
  
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
    
  stream.precision(precision);
  stream.width(width);
  stream.flags(flags); 
}

template<class T>
void Matrix<T>::printMatrix(std::ostream &stream) const{

  std::streamsize precision = stream.precision();
  std::streamsize width = stream.width();
  std::ios_base::fmtflags flags = stream.flags();
  
  stream.precision(stream_precision_);
  stream.width(stream_width_);
  if (stream_scientific_) {
    stream.setf(std::ios::scientific);
  } else {
    stream.unsetf(std::ios::scientific);
  }
  
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
  
    
  stream.precision(precision);
  stream.width(width);
  stream.flags(flags);
}

template<class T>
void Matrix<T>::printDense(std::ostream &stream) const{
  stream << className() << "(rows = " << size1() << ", cols = " << size2() << "):" << std::endl;
  printMatrix(stream);
}


template<class T>
void Matrix<T>::printSparse(std::ostream &stream) const {
  stream << className() << "(rows = " << size1() << ", cols = " << size2() << ", nnz = " << size() << ")";
  if (size()>0) stream << ":" << std::endl;
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
Matrix<T>::Matrix(double val) : sparsity_(CRSSparsity(1,1,true)), data_(std::vector<T>(1,val)) {
}

template<class T>
Matrix<T>::Matrix(int n, int m, const std::vector<int>& col, const std::vector<int>& rowind, const std::vector<T>& d) : sparsity_(CRSSparsity(n,m,col,rowind)), data_(d){
  if(data_.size() != sparsity_.size())
    data_.resize(sparsity_.size()); // Why not throw an error?
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
Matrix<T>::Matrix(const CRSSparsity& sparsity, const T& val) : sparsity_(sparsity), data_(std::vector<T>(sparsity.size(),val)){
}

template<class T>
Matrix<T>::Matrix(const CRSSparsity& sparsity, const std::vector<T>& d) : sparsity_(sparsity), data_(d) {
  casadi_assert_message(sparsity.size()==d.size(),"Size mismatch." << std::endl << "You supplied a sparsity of " << sparsity.dimString() << ", but the supplied vector is of length " << d.size());
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
Matrix<T> Matrix<T>::unary(int op, const Matrix<T> &x){
  // Return value
  Matrix<T> ret(x.sparsity());
  
  // Nonzeros
  std::vector<T>& ret_data = ret.data();
  const std::vector<T>& x_data = x.data();
  
  // Do the operation on all non-zero elements
  for(int el=0; el<x.size(); ++el){
    casadi_math<T>::fun(op,x_data[el],x_data[el],ret_data[el]);
  }

  // Check the value of the structural zero-entries, if there are any
  if(!x.dense() && !operation_checker<F0XChecker>(op)){
    // Get the value for the structural zeros
    T fcn_0;
    casadi_math<T>::fun(op,0,0,fcn_0);
    if(!casadi_limits<T>::isZero(fcn_0)){ // Remove this if?
      ret.makeDense(ret.size1(),ret.size2(),fcn_0);
    }
  }
    
  return ret;
}

template<class T>
Matrix<T> Matrix<T>::operator-() const{
  return unary(OP_NEG,*this);
}

template<class T>
Matrix<T> Matrix<T>::operator+() const{
  return *this;
}

template<class T>
Matrix<T> Matrix<T>::__add__(const Matrix<T> &y) const{
  return binary(OP_ADD,*this,y);
}

template<class T>
Matrix<T> Matrix<T>::__sub__(const Matrix<T> &y) const{
  return binary(OP_SUB,*this,y);
}

template<class T>
Matrix<T> Matrix<T>::__mul__(const Matrix<T> &y) const{
  return binary(OP_MUL,*this,y);
}

template<class T>
Matrix<T> Matrix<T>::__div__(const Matrix<T> &y) const{
  return binary(OP_DIV,*this,y);
}

template<class T>
Matrix<T> Matrix<T>::__lt__(const Matrix<T> &y) const{
  return binary(OP_LT,*this,y);
}

template<class T>
Matrix<T> Matrix<T>::__le__(const Matrix<T> &y) const{
  return binary(OP_LE,*this,y);
}

template<class T>
Matrix<T> Matrix<T>::__eq__(const Matrix<T> &y) const{
  return binary(OP_EQ,*this,y);
}

template<class T>
Matrix<T> Matrix<T>::__ne__(const Matrix<T> &y) const{
  return binary(OP_NE,*this,y);
}

template<class T>
Matrix<T> Matrix<T>::__mrdivide__(const Matrix<T>& b) const { if (b.numel()==1) return *this/b; throw CasadiException("mrdivide: Not implemented");}

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
        std::fill(data().begin(),data().end(),val);
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
  sparsity().set(getPtr(data()),getPtr(val.data()),val.sparsity());
}

template<class T>
void Matrix<T>::get(Matrix<T>& val, Sparsity sp) const{
  val.set(*this,sp);
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
  if(sp==SPARSE || (sp==DENSE && dense())){
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
    casadi_error("Matrix<T>::getArray: not SPARSE, SPARSESYM  or DENSE");
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
  } else if(sp==SPARSESYM) {
    std::vector<int> mapping;
    sparsity().transpose(mapping,false);
    // copy to the result vector
    int nz = 0;
    for(int row=0; row<size1(); ++row){
      // Loop over the elements in the row
      for(int el=rowind(row); el<rowind(row+1); ++el){ // loop over the non-zero elements
        if(col(el) > row) break; // break inner loop (only lower triangular part is used)
        v[el] = val[nz++];
        v[mapping[el]] = v[el];
      }
    }
  } else {
    throw CasadiException("Matrix<T>::setArray: not SPARSE, SPARSESYM or DENSE");
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
Matrix<T> Matrix<T>::__pow__(const Matrix<T>& y) const{
  return binary(OP_POW,*this,y);
}

template<class T>
Matrix<T> Matrix<T>::__constpow__(const Matrix<T>& y) const{
  return binary(OP_CONSTPOW,*this,y);
}

template<class T>
Matrix<T> Matrix<T>::sin() const{
  return unary(OP_SIN,*this);
}

template<class T>
Matrix<T> Matrix<T>::cos() const{
  return unary(OP_COS,*this);
}

template<class T>
Matrix<T> Matrix<T>::tan() const{
  return unary(OP_TAN,*this);
}

template<class T>
Matrix<T> Matrix<T>::erf() const{
  return unary(OP_ERF,*this);
}

template<class T>
Matrix<T> Matrix<T>::arcsin() const{
  return unary(OP_ASIN,*this);
}

template<class T>
Matrix<T> Matrix<T>::arccos() const{
  return unary(OP_ACOS,*this);
}

template<class T>
Matrix<T> Matrix<T>::arctan() const{
  return unary(OP_ATAN,*this);
}

template<class T>
Matrix<T> Matrix<T>::sinh() const{
  return unary(OP_SINH,*this);
}

template<class T>
Matrix<T> Matrix<T>::cosh() const{
  return unary(OP_COSH,*this);
}

template<class T>
Matrix<T> Matrix<T>::tanh() const{
  return unary(OP_TANH,*this);
}

template<class T>
Matrix<T> Matrix<T>::arcsinh() const{
  return unary(OP_ASINH,*this);
}

template<class T>
Matrix<T> Matrix<T>::arccosh() const{
  return unary(OP_ACOSH,*this);
}

template<class T>
Matrix<T> Matrix<T>::arctanh() const{
  return unary(OP_ATANH,*this);
}

template<class T>
Matrix<T> Matrix<T>::exp() const{
  return unary(OP_EXP,*this);
}

template<class T>
Matrix<T> Matrix<T>::log() const{
  return unary(OP_LOG,*this);
}

template<class T>
Matrix<T> Matrix<T>::log10() const{
  return log()*(1/std::log(10.));
}

template<class T>
Matrix<T> Matrix<T>::sqrt() const{
  return unary(OP_SQRT,*this);
}

template<class T>
Matrix<T> Matrix<T>::floor() const{
  return unary(OP_FLOOR,*this);
}

template<class T>
Matrix<T> Matrix<T>::ceil() const{
  return unary(OP_CEIL,*this);
}

template<class T>
Matrix<T> Matrix<T>::fabs() const{
  return unary(OP_FABS,*this);
}

template<class T>
Matrix<T> Matrix<T>::sign() const{
  return unary(OP_SIGN,*this);
}

template<class T>
Matrix<T> Matrix<T>::erfinv() const{
  return unary(OP_ERFINV,*this);
}

template<class T>
Matrix<T> Matrix<T>::fmin(const Matrix<T>& y) const{
  return binary(OP_FMIN,*this,y);
}

template<class T>
Matrix<T> Matrix<T>::arctan2(const Matrix<T>& y) const{
  return binary(OP_ATAN2,*this,y);
}

template<class T>
Matrix<T> Matrix<T>::fmax(const Matrix<T>& y) const{
  return binary(OP_FMAX,*this,y);
}

template<class T>
Matrix<T> Matrix<T>::printme(const Matrix<T>& y) const{
  return binary(OP_PRINTME,*this,y);
}

template<class T>
Matrix<T> Matrix<T>::logic_not() const{
  return unary(OP_NOT,*this);
}

template<class T>
Matrix<T> Matrix<T>::logic_and(const Matrix<T>& y) const{
  return binary(OP_AND,*this,y);
}

template<class T>
Matrix<T> Matrix<T>::logic_or(const Matrix<T>& y) const{
  return binary(OP_OR,*this,y);
}

template<class T>
Matrix<T> Matrix<T>::if_else_zero(const Matrix<T>& y) const{
  return binary(OP_IF_ELSE_ZERO,*this,y);
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
void Matrix<T>::remove(const std::vector<int>& ii, const std::vector<int>& jj) {
  if (!inBounds(ii,size1())) {
    casadi_error("Remove(ii,jj) out of bounds. Your ii contains " << *std::min_element(ii.begin(),ii.end()) << " up to " << *std::max_element(ii.begin(),ii.end()) << ", which is outside of the matrix shape " << dimString() << ".");
  }
  if (!inBounds(jj,size2())) {
    casadi_error("Remove(ii,jj) out of bounds. Your jj contains " << *std::min_element(jj.begin(),jj.end()) << " up to " << *std::max_element(jj.begin(),jj.end()) << ", which is outside of the matrix shape " << dimString() << ".");
  }
  
  // Remove by performing a complementary slice
  std::vector<int> iic = complement(ii,size1());
  std::vector<int> jjc = complement(jj,size2());
  
  Matrix<T> ret = operator()(iic,jjc);
  
  operator=(ret);
  
}

template<class T>
void Matrix<T>::enlarge(int nrow, int ncol, const std::vector<int>& ii, const std::vector<int>& jj){
  sparsityRef().enlarge(nrow,ncol,ii,jj);
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
Matrix<T> Matrix<T>::mul(const Matrix<T> &y) const {
  return this->mul_smart(y);
}

template<class T>
Matrix<T> Matrix<T>::mul_full(const Matrix<T> &y) const{
  // First factor
  const Matrix<T>& x = *this;
  
  // Return object (assure RVO)
  Matrix<T> ret;

  // Matrix multiplication

  // Form the transpose of y
  Matrix<T> y_trans = y.trans();

  // Create the sparsity pattern for the matrix-matrix product
  CRSSparsity spres = x.sparsity().patternProduct(y_trans.sparsity());

  // Create the return object
  ret = Matrix<T>(spres, 0);

  // Carry out the matrix product
  mul_no_alloc_nt(x,y_trans,ret);
  
  return ret;
}

template<class T>
void Matrix<T>::mul_no_alloc_nn(const Matrix<T> &x, const Matrix<T> &y, Matrix<T>& z){
  // Assert dimensions
  casadi_assert(x.size1()==z.size1());
  casadi_assert(y.size2()==z.size2());
  casadi_assert(x.size2()==y.size1());
  
  // Direct access to the arrays
  const std::vector<int> &x_rowind = x.rowind();
  const std::vector<int> &x_col = x.col();
  const std::vector<T> &x_data = x.data();
  const std::vector<int> &y_rowind = y.rowind();
  const std::vector<int> &y_col = y.col();
  const std::vector<T> &y_data = y.data();
  const std::vector<int> &z_rowind = z.rowind();
  const std::vector<int> &z_col = z.col();
  std::vector<T> &z_data = z.data();

  // loop over the rows of the first argument
  for(int i=0; i<x_rowind.size()-1; ++i){
    for(int el=x_rowind[i]; el<x_rowind[i+1]; ++el){ // loop over the non-zeros of the first argument
      int j = x_col[el];
      int el1 = z_rowind[i];
      int el2 = y_rowind[j];
      while(el1 < z_rowind[i+1] && el2 < y_rowind[j+1]){ // loop over matching non-zero elements
        int j1 = z_col[el1];
        int i2 = y_col[el2];      
        if(j1==i2){
          z_data[el1++] += x_data[el]*y_data[el2++];
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
  void Matrix<T>::mul_no_alloc_nn(const Matrix<T> &x, const std::vector<T> &y, std::vector<T>& z){
    // Assert dimensions
    casadi_assert(x.size1()==z.size());
    casadi_assert(x.size2()==y.size());
    
    // Direct access to the arrays
    const std::vector<int> &x_rowind = x.rowind();
    const std::vector<int> &x_col = x.col();
    const std::vector<T> &x_data = x.data();
    
    // loop over the rows of the matrix
    for(int i=0; i<x_rowind.size()-1; ++i){
      for(int el=x_rowind[i]; el<x_rowind[i+1]; ++el){ // loop over the non-zeros of the matrix
        int j = x_col[el];
        
        // Perform operation
        z[i] += x_data[el] * y[j];
      }
    }
  }
  
  template<class T>
  void Matrix<T>::mul_no_alloc_tn(const Matrix<T>& x_trans, const std::vector<T> &y, std::vector<T> &z){
    // Assert dimensions
    casadi_assert(x_trans.size2()==z.size());
    casadi_assert(x_trans.size1()==y.size());
    
    // Direct access to the arrays
    const std::vector<int> &x_colind = x_trans.rowind();
    const std::vector<int> &x_row = x_trans.col();
    const std::vector<T> &x_trans_data = x_trans.data();
    
    // loop over the columns of the matrix
    for(int i=0; i<x_colind.size()-1; ++i){
      for(int el=x_colind[i]; el<x_colind[i+1]; ++el){ // loop over the non-zeros of the matrix
        int j = x_row[el];
        z[j] += x_trans_data[el] * y[i];
      }
    }
  }
  
template<class T>
void Matrix<T>::mul_no_alloc_tn(const Matrix<T>& x_trans, const Matrix<T> &y, Matrix<T> &z){
  // Assert dimensions
  casadi_assert(x_trans.size2()==z.size1());
  casadi_assert(y.size2()==z.size2());
  casadi_assert(x_trans.size1()==y.size1());
  
  // Direct access to the arrays
  const std::vector<int> &x_colind = x_trans.rowind();
  const std::vector<int> &x_row = x_trans.col();
  const std::vector<T> &x_trans_data = x_trans.data();
  const std::vector<int> &y_rowind = y.rowind();
  const std::vector<int> &y_col = y.col();
  const std::vector<T> &y_data = y.data();
  const std::vector<int> &z_rowind = z.rowind();
  const std::vector<int> &z_col = z.col();
  std::vector<T> &z_data = z.data();

  // loop over the columns of the first argument
  for(int i=0; i<x_colind.size()-1; ++i){
    for(int el=x_colind[i]; el<x_colind[i+1]; ++el){ // loop over the non-zeros of the first argument
      int j = x_row[el];
      int el1 = y_rowind[i];
      int el2 = z_rowind[j];
      while(el1 < y_rowind[i+1] && el2 < z_rowind[j+1]){ // loop over matching non-zero elements
        int j1 = y_col[el1];
        int i2 = z_col[el2];      
        if(j1==i2){
          z_data[el2++] += x_trans_data[el] * y_data[el1++];
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
void Matrix<T>::mul_no_alloc_nt(const Matrix<T> &x, const Matrix<T> &y_trans, Matrix<T>& z){
  // Assert dimensions
  casadi_assert(x.size1()==z.size1());
  casadi_assert(y_trans.size1()==z.size2());
  casadi_assert(x.size2()==y_trans.size2());
  
  // Direct access to the arrays
  const std::vector<int> &x_rowind = x.rowind();
  const std::vector<int> &x_col = x.col();
  const std::vector<T> &x_data = x.data();
  const std::vector<int> &y_colind = y_trans.rowind();
  const std::vector<int> &y_row = y_trans.col();
  const std::vector<T> &y_trans_data = y_trans.data();
  const std::vector<int> &z_rowind = z.rowind();
  const std::vector<int> &z_col = z.col();
  std::vector<T> &z_data = z.data();

  #ifdef WITH_EIGEN3
  // NOTE: this doesn't belong here. It should be put at some higher level. Also, is this really fast than the implementation below?
  if (x.dense() && y_trans.dense() && z.dense()) {
    Eigen::Map< const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic , Eigen::RowMajor > > X(&x_data[0],x.size1(),x.size2());
    Eigen::Map< const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic > > Y(&y_trans_data[0],y_trans.size2(),y_trans.size1());
    Eigen::Map< Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic , Eigen::RowMajor> > Z(&z.data()[0],z.size1(),z.size2());
    Z = X*Y;
    return;
  }
  #endif
  
  // loop over the rows of the resulting matrix
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
template<bool Fwd>
void Matrix<T>::mul_sparsity(Matrix<T> &x, Matrix<T> &y_trans, Matrix<T>& z){
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
      while(el1 < x_rowind[i+1] && el2 < y_colind[j+1]){ // loop over non-zero elements
        int j1 = x_col[el1];
        int i2 = y_row[el2];      
        if(j1==i2){
          // | and not & since we are propagating dependencies
          if(Fwd){
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
T Matrix<T>::quad_form(const Matrix<T>& A, const std::vector<T>& x){
  // Assert dimensions
  casadi_assert(x.size()==A.size1() && x.size()==A.size2());
  
  // Access the internal data of A
  const std::vector<int> &A_rowind = A.rowind();
  const std::vector<int> &A_col = A.col();
  const std::vector<T> &A_data = A.data();
  
  // Return value
  T ret=0;

  // Loop over the rows of A
  for(int i=0; i<x.size(); ++i){
    // Loop over the nonzeros of A
    for(int el=A_rowind[i]; el<A_rowind[i+1]; ++el){
      // Get column
      int j = A_col[el];
      
      // Add contribution
      ret += x[i]*A_data[el]*x[j];
    }
  }
  
  return ret;
}

template<class T>
Matrix<T> Matrix<T>::trans() const{
  // quick return if empty or scalar
  if((size1()==0 && size2()==0) || scalar()) return *this;

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
Matrix<T> Matrix<T>::scalar_matrix(int op, const Matrix<T> &x, const Matrix<T> &y){
  // Return value
  Matrix<T> ret(y.sparsity());
  
  // Nonzeros
  std::vector<T>& ret_data = ret.data();
  const std::vector<T>& x_data = x.data();
  const T& x_val = x_data.empty() ? casadi_limits<T>::zero : x.front();
  const std::vector<T>& y_data = y.data();
  
  // Do the operation on all non-zero elements
  for(int el=0; el<y.size(); ++el){
    casadi_math<T>::fun(op,x_val,y_data[el],ret_data[el]);
  }

  // Check the value of the structural zero-entries, if there are any
  if(!y.dense() && !operation_checker<FX0Checker>(op)){
    // Get the value for the structural zeros
    T fcn_0;
    casadi_math<T>::fun(op,x_val,casadi_limits<T>::zero,fcn_0);
    if(!casadi_limits<T>::isZero(fcn_0)){ // Remove this if?
      ret.makeDense(ret.size1(),ret.size2(),fcn_0);
    }
  }
    
  return ret;
}

template<class T>
Matrix<T> Matrix<T>::matrix_scalar(int op, const Matrix<T> &x, const Matrix<T> &y){
  // Return value
  Matrix<T> ret(x.sparsity());
  
  // Nonzeros
  std::vector<T>& ret_data = ret.data();
  const std::vector<T>& x_data = x.data();
  const std::vector<T>& y_data = y.data();
  const T& y_val = y_data.empty() ? casadi_limits<T>::zero : y.front();
  
  // Do the operation on all non-zero elements
  for(int el=0; el<x.size(); ++el){
    casadi_math<T>::fun(op,x_data[el],y_val,ret_data[el]);
  }

  // Check the value of the structural zero-entries, if there are any
  if(!x.dense() && !operation_checker<F0XChecker>(op)){
    // Get the value for the structural zeros
    T fcn_0;
    casadi_math<T>::fun(op,casadi_limits<T>::zero,y_val,fcn_0);
    if(!casadi_limits<T>::isZero(fcn_0)){ // Remove this if?
      ret.makeDense(ret.size1(),ret.size2(),fcn_0);
    }
  }
    
  return ret;
}

template<class T>
Matrix<T> Matrix<T>::matrix_matrix(int op, const Matrix<T> &x, const Matrix<T> &y){

  if (!(x.size1() == y.size1() && x.size2() == y.size2())) {
    std::stringstream ss;
    casadi_math<T>::print(op,ss,"lhs","rhs");
    casadi_error("matrix_matrix: dimension mismatch in element-wise matrix operation " << ss.str() <<"." << std::endl << "Left argument has shape " << x.dimString() << ", right has shape " << y.dimString() << ". They should be equal."
    ); 
  }

  // Get the sparsity pattern of the result (ignoring structural zeros giving rise to nonzero result)
  const CRSSparsity& x_sp = x.sparsity();
  const CRSSparsity& y_sp = y.sparsity();
  CRSSparsity r_sp = x_sp.patternCombine(y_sp, operation_checker<F0XChecker>(op), operation_checker<FX0Checker>(op));

  // Return value
  Matrix<T> r(r_sp);
  
  // Perform the operations elementwise
  if(x_sp==y_sp){
    // Matching sparsities
    casadi_math<T>::fun(op,getPtr(x.data()),getPtr(y.data()),getPtr(r.data()),r_sp.size());
  } else if(y_sp==r_sp){
    // Project first argument
    Matrix<T> x_mod = x(r_sp);
    casadi_math<T>::fun(op,getPtr(x_mod.data()),getPtr(y.data()),getPtr(r.data()),r_sp.size());
  } else if(x_sp==r_sp){
    // Project second argument
    Matrix<T> y_mod = y(r_sp);
    casadi_math<T>::fun(op,getPtr(x.data()),getPtr(y_mod.data()),getPtr(r.data()),r_sp.size());
  } else {
    // Project both arguments
    Matrix<T> x_mod = x(r_sp);
    Matrix<T> y_mod = y(r_sp);
    casadi_math<T>::fun(op,getPtr(x_mod.data()),getPtr(y_mod.data()),getPtr(r.data()),r_sp.size());
  }

  // Handle structural zeros giving rise to nonzero result, e.g. cos(0) == 1
  if(!r.dense() && !operation_checker<F00Checker>(op)){
    // Get the value for the structural zeros
    T fcn_0;
    casadi_math<T>::fun(op,casadi_limits<T>::zero,casadi_limits<T>::zero,fcn_0);
    r.makeDense(r.size1(),r.size2(),fcn_0);
  }
  
  return r;
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
Matrix<T> Matrix<T>::sparse(const std::vector<int>& row, const std::vector<int>& col, const std::vector<T>& d) {
  return sparse(row,col,d,*std::max_element(row.begin(),row.end()),*std::max_element(col.begin(),col.end()));
}

template<class T>
Matrix<T> Matrix<T>::sparse(const std::vector<int>& row, const std::vector<int>& col, const std::vector<T>& d, const std::pair<int,int>& nm) {
  return sparse(row,col,d,nm.first,nm.second);
}

template<class T>
Matrix<T> Matrix<T>::sparse(const std::vector<int>& row, const std::vector<int>& col, const std::vector<T>& d, int n, int m) {
  casadi_assert_message(row.size()==col.size() && row.size()==d.size(),"Argument error in Matrix<T>::sparse(row,col,d): supplied lists must all be of equal length, but got: " << row.size() << ", " << col.size()  << " and " << d.size());
  std::vector<int> mapping;
  Matrix<T> ret(sp_triplet(n,m,row,col,mapping),0);
  
  for (int k=0;k<mapping.size();++k) ret.data()[k] = d[mapping[k]];
  
  return ret;
}

template<class T>
Matrix<T> Matrix<T>::zeros(const CRSSparsity& sp){
  return Matrix<T>(sp,0);
}

template<class T>
Matrix<T> Matrix<T>::zeros(const std::pair<int,int> &nm){
  return zeros(nm.first,nm.second);
}

template<class T>
Matrix<T> Matrix<T>::zeros(int n, int m){
  return zeros(sp_dense(n,m));
}

template<class T>
Matrix<T> Matrix<T>::ones(const CRSSparsity& sp){
  return Matrix<T>(sp,1);
}

template<class T>
Matrix<T> Matrix<T>::ones(const std::pair<int,int> &nm){
  return ones(nm.first,nm.second);
}

template<class T>
Matrix<T> Matrix<T>::ones(int n, int m){
  return ones(sp_dense(n,m));
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
Matrix<T> Matrix<T>::inf(const CRSSparsity& sp){
  casadi_assert_message(std::numeric_limits<T>::has_infinity,"Datatype cannot represent infinity");
  return Matrix<T>(sp,std::numeric_limits<T>::infinity());
}


template<class T>
Matrix<T> Matrix<T>::inf(const std::pair<int,int>& nm){
  return inf(nm.first, nm.second);
}

template<class T>
Matrix<T> Matrix<T>::inf(int n, int m){
  return inf(sp_dense(n,m));
}

template<class T>
Matrix<T> Matrix<T>::nan(const CRSSparsity& sp){
  casadi_assert_message(std::numeric_limits<T>::has_quiet_NaN,"Datatype cannot represent not-a-number");
  return Matrix<T>(sp,std::numeric_limits<T>::quiet_NaN());
}

template<class T>
Matrix<T> Matrix<T>::nan(const std::pair<int,int>& nm){
  return nan(nm.first, nm.second);
}

template<class T>
Matrix<T> Matrix<T>::nan(int n, int m){
  return nan(sp_dense(n,m));
}

template<class T>
void Matrix<T>::append(const Matrix<T>& y){

  // Quick return if expr is empty
  if(size1()==0 && size2()==0){
    *this=y;
    return;
  }

  // Quick return if empty
  if(y.size1()==0 && y.size2()==0) return;
  
  // Append the sparsity pattern
  sparsityRef().append(y.sparsity());
  
  // Add the non-zeros
  data().insert(end(),y.begin(),y.end());
}

template<class T>
NonZeroIterator<T>::NonZeroIterator(const Matrix<T> & m) 
                     : m_(m) {
  nz.i  = 0;
  nz.j  = 0;
  nz.k  = 0;
}

template<class T>
bool NonZeroIterator<T>::operator==(const NonZeroIterator<T>& rhs) {return (m_ == rhs.m_) && (nz.k==rhs.nz.k);}

template<class T>
NonZero<T>& NonZeroIterator<T>::operator*() {
  
  return nz;
}
    
template<class T>
NonZeroIterator<T>& NonZeroIterator<T>::operator++() {
  nz.k ++;
  
  if (nz.k < m_.size()) {
    nz.j = m_.col()[nz.k];
    nz.el = m_.data()[nz.k];
    while (nz.k>=m_.rowind(nz.i)) {nz.i++; }
  }
  return *this;
}

template<class T>
NonZeroIterator<T> NonZeroIterator<T>::begin() {
  NonZeroIterator<T> it = NonZeroIterator<T>(m_);
  return it;

}
template<class T>
NonZeroIterator<T> NonZeroIterator<T>::end() {
  NonZeroIterator<T> it = NonZeroIterator<T>(m_);
  it.nz.k = m_.size()-1;
  return it;
}


} // namespace CasADi

#endif // MATRIX_IMPL_HPP

