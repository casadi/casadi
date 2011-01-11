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

#ifndef MATRIX_HPP
#define MATRIX_HPP

#include <vector>
#include "../casadi_exception.hpp"
#include "../printable_object.hpp"
#include "../casadi_limits.hpp"
#include "element.hpp"
#include "crs_sparsity.hpp"

namespace CasADi{
  
  /** Sparsity format for getting and setting inputs and outputs */
  enum Sparsity{SPARSE,SPARSESYM,DENSE,DENSESYM};

  /** \brief General sparse matrix class
  \author Joel Andersson 
  \date 2010	
*/
template<class T>
class Matrix : public std::vector<T>, public PrintableObject{

  public:
    /** \brief  constructors */
    /// empty 0-by-0 matrix constructor
    Matrix();
    
    /// Copy constructor (normal)
    Matrix(const Matrix<T>& m);
    
    #ifndef SWIG
    /// Copy constructor (possible swap)
    Matrix(Matrix<T>& m);
    
    /// Assignment (normal)
    Matrix<T>& operator=(const Matrix<T>& m);
    
    /// Assignment (possible swap)
    Matrix<T>& operator=(Matrix<T>& m);
    #endif // SWIG
    
    /// empty n-by-m matrix constructor
    Matrix(int n, int m);
    
    /// dense n-by-m matrix filled with val constructor
    Matrix(int n, int m, const T& val);

    /// sparse n-by-m matrix filled with given sparsity
    Matrix(int n, int m, const std::vector<int>& col, const std::vector<int>& rowind);

    /// sparse matrix with a given sparsity
    explicit Matrix(const CRSSparsity& sparsity);
    
    /// This constructor enables implicit type conversion from a scalar type
    Matrix(const T &val);

#ifndef SWIG    
    /** \brief  Create an expression from an stl vector  */
    template<typename A>
    Matrix(const std::vector<A>& x){
      sparsity_ = CRSSparsity(x.size(),1,true);
      std::vector<T>::resize(x.size());
      copy(x.begin(),x.end(),std::vector<T>::begin());
    }

    /** \brief  Create a non-vector expression from an stl vector */
    template<typename A>
    Matrix(const std::vector<A>& x,  int n, int m){
      if(x.size() != n*m) throw CasadiException("Matrix::Matrix(const std::vector<T>& x,  int n, int m): dimension mismatch");
      sparsity_ = CRSSparsity(n,m,true);
      std::vector<T>::resize(x.size());
      copy(x.begin(),x.end(),std::vector<T>::begin());
    }
#endif // SWIG

    // get the number of non-zeros
    //int size() const;        

    /// get the number of elements
    int numel() const;

    /// get the first dimension
    int size1() const;       
    
    /// get the second dimension
    int size2() const;

    //@{
    /// Check type of matrix
    bool empty() const; // is the matrix empty
    bool scalar() const; // is the matrix scalar
    bool vector() const; // is the matrix a vector
    //@}


#ifndef SWIG    
    /// get an element
    const T getElement(int i=0, int j=0) const;
    
    /// set an element
    void setElement(int i, int j, const T& el);

    /// get a reference to an element
    T& getElementRef(int i=0, int j=0);
  
    /// Access an element 
    Element<Matrix<T>,T> operator()(int i, int j=0){ return Element<Matrix<T>,T>(*this,i,j); }

    /// Const access an element 
    const T operator()(int i, int j=0) const{ return getElement(i,j); }
#endif // SWIG

#if 0
    /// Get a non-zero entry
    const T getitem(int k) const;
    
    /// Get a matrix entry
    const T getitem(int I[2]) const;
    
    /// Get a slice
    const T getslice(int start1, int stop1, int stride1, int start2, int stop2, int stride2);
    
    /// Set a non-zero entry
    void setitem(int k, const T& el);
    
    /// Set a matrix entry
    void setitem(int i, int j, const T&  el);

    /// Set a slice
    void setslice(int start1, int stop1, int stride1, int start2, int stop2, int stride2, const std::vector<T>& el);
#endif    


    /// Python: get a non-zero entry
    const T __getitem__(int i) const;
    
    /// Python: get a matrix entry
    const T __getitem__(const std::vector<int> &I) const;
    
    /// Python: set a non-zero entry
    void __setitem__(int k, const T& el);
    
    /// Python: set a matrix entry
    void __setitem__(const std::vector<int> &I, const T&  el);

    /** \brief  Make the matrix an dense n-by-m matrix */
    void makeDense(int n, int m, const T& val);

    /** \brief  Make the matrix an empty n-by-m matrix */
    void makeEmpty(int n, int m);

    /** \brief  Unary function */
    void unary(T (*fcn)(const T&), const Matrix<T>& x);
    void binary(T (*fcn)(const T&, const T&), const Matrix<T> &x, const Matrix<T> &y);
    void matrix_matrix(T (*fcn)(const T&, const T&), const Matrix<T>& x, const Matrix<T>& y);
    void matrix_scalar(T (*fcn)(const T&, const T&), const Matrix<T>& x, const T& y);
    void scalar_matrix(T (*fcn)(const T&, const T&), const T& x, const Matrix<T>& y);

    Matrix<T> operator+() const;
    Matrix<T> operator-() const;
    
    Matrix<T> operator+(const Matrix<T> &y) const;
    Matrix<T> operator-(const Matrix<T> &y) const;
    Matrix<T> operator*(const Matrix<T> &y) const;
    Matrix<T> operator/(const Matrix<T> &y) const;

    Matrix<T>& operator+=(const Matrix<T> &y);
    Matrix<T>& operator-=(const Matrix<T> &y);
    Matrix<T>& operator*=(const Matrix<T> &y);
    Matrix<T>& operator/=(const Matrix<T> &y);
    
    
    //@{
    /// Printing
#ifndef SWIG
    virtual void print(std::ostream &stream=std::cout) const; // print default style
#endif
    void printScalar(std::ostream &stream=std::cout) const; // print scalar
    void printVector(std::ostream &stream=std::cout) const; // print vector-style
    void printMatrix(std::ostream &stream=std::cout) const; // print matrix-style
    //@}

    /** \brief Swap the the vector content upon the next copy assignment, saves a copy operation: 
    After the copy operation has taken place, the object can not be used for anything,
    as its vector is no longer of correct length. Only the destructor is allowed to get called.
    Make sure that is the function is only called e.g. for temporary objects at the end of the 
    scope, before copying to a return value. Caution is adviced.
    */
    void swapOnCopy();

    // Get the sparsity pattern
    const std::vector<int>& col() const;
    const std::vector<int>& rowind() const;
    std::vector<int>& col();
    std::vector<int>& rowind();
    int col(int el) const;
    int rowind(int row) const;
    
    std::string __repr__() { return getRepresentation(); }
    
    void clear();
    void resize(int n, int m);
    void reserve(int nnz);
    void reserve(int nnz, int nrow);
    
    // Access the sparsity
    const CRSSparsity& sparsity() const;
    
    // The following need cleaning up

    /** \brief  Set the non-zero elements, scalar */
    void set(T val, Sparsity sp=SPARSE);
    
    /** \brief  Get the non-zero elements, scalar */
    void get(T& val, Sparsity sp=SPARSE) const;

    /** \brief  Set the non-zero elements, vector */
    void set(const std::vector<T>& val, Sparsity sp=SPARSE);

    /** \brief  Get the non-zero elements, vector */
    void get(std::vector<T>& val, Sparsity sp=SPARSE) const;

    /** \brief  Set the non-zero elements, matrix */
    //void set(const Matrix<T>& val);

    /** \brief  Get the non-zero elements, matrix */
    //void get(Matrix<T>& val) const;

    #ifndef SWIG
    /** \brief  Get the non-zero elements, array */
    void get(T* val, Sparsity sp=SPARSE) const;    

    /** \brief  Set the non-zero elements, array */
    void set(const T* val, Sparsity sp=SPARSE);
    #endif

    /** \brief  Get the result */
    void getSparseSym(T *res) const;    // general sparse, symmetric matrix

    /** \brief  Get the result times a vector */
    void getTimesVector(const T *v, T *res) const;

    /** \brief  Save the result to the LAPACK banded format -- see LAPACK documentation 
    kl:    The number of subdiagonals in res 
    ku:    The number of superdiagonals in res 
    ldres: The leading dimension in res 
    res:   The number of superdiagonals */
    void getBand(int kl, int ku, int ldres, T *res) const;

    // Make sure that the number of non-zeros is sz
    void assertNNZ(int sz, Sparsity sp) const;
    void assertNumEl(int sz) const;
    
    
  private:
    /// Sparsity of the matrix in a compressed row storage (CRS) format
    CRSSparsity sparsity_;
    
    /// Swap the content of the vector upon the next copy action
    bool swap_on_copy_;

};

// #ifdef SWIG
// %extend Matrix<T>{
// T __getitem__(const std::vector<PyObject*> &I ) {
//   throw CasadiException("ok!");
// }
// }
// // 
// #endif // SWIG



// #ifdef SWIG
// %extend Matrix<T> {
// std::string __str__() { return $self->getDescription(); }
// std::string __repr__() { return $self->getRepresentation(); }
// }
// #endif // SWIG

#ifndef SWIG
// Implementations

template<class T>
const T Matrix<T>::getElement(int i, int j) const{
  int ind = sparsity_.getNZ(i,j);
  if(ind==-1)
    return 0;
  else
    return std::vector<T>::at(ind);
}

template<class T>
void Matrix<T>::setElement(int i, int j, const T& el){
  getElementRef(i,j) = el;
}

template<class T>
T& Matrix<T>::getElementRef(int i, int j){
  int oldsize = sparsity_.size();
  int ind = sparsity_.getNZ(i,j);
  if(oldsize != sparsity_.size())
    std::vector<T>::insert(std::vector<T>::begin()+ind,0);
  return std::vector<T>::at(ind);
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
void Matrix<T>::makeDense(int n, int m, const T& val){
  if(n*m != numel()){
    sparsity_ = CRSSparsity(n,m,true);
    std::vector<T>::clear();
    std::vector<T>::resize(n*m, val);
  }
}

template<class T>
bool Matrix<T>::empty() const{
  return numel()==0;
}

template<class T>
bool Matrix<T>::scalar() const{
  return size1()==1 && size2()==1;
}

template<class T>
bool Matrix<T>::vector() const{
  return size2()==1;
}

template<class T>
Matrix<T>::Matrix(){
  swap_on_copy_ = false;
  sparsity_ = CRSSparsity(0,0,false);
}

template<class T>
Matrix<T>::Matrix(const Matrix<T>& m){
  swap_on_copy_ = false;
  sparsity_ = m.sparsity_;
  static_cast<std::vector<T>&>(*this) = m;
}

template<class T>
Matrix<T>::Matrix(Matrix<T>& m){
  swap_on_copy_ = false;
  if(m.swap_on_copy_){
    // Swap the vector with m
    m.swap(*this);
    m.sparsity_.swap(sparsity_);
  } else {
    // Copy the content
    static_cast<std::vector<T>&>(*this) = m;
    sparsity_ = m.sparsity_; // shallow copy!
  }
}

template<class T>
Matrix<T>& Matrix<T>::operator=(const Matrix<T>& m){
  sparsity_ = m.sparsity_;
  static_cast<std::vector<T>&>(*this) = m;
}

template<class T>
Matrix<T>& Matrix<T>::operator=(Matrix<T>& m){
  if(m.swap_on_copy_){
    // Swap the vector with m
    m.swap(*this);
    m.sparsity_.swap(sparsity_);
  } else {
    // Copy the content
    static_cast<std::vector<T>&>(*this) = m;
    sparsity_ = m.sparsity_; // shallow copy!
  }
}

template<class T>
void Matrix<T>::swapOnCopy(){
  swap_on_copy_ = true;
}

template<class T>
Matrix<T>::Matrix(int n, int m){
  swap_on_copy_ = false;
  sparsity_ = CRSSparsity(n,m,false);
}

template<class T>
Matrix<T>::Matrix(int n, int m, const T& val){
  swap_on_copy_ = false;
  sparsity_ = CRSSparsity(n,m,true);
  std::vector<T>::resize(n*m, val);
}

template<class T>
void Matrix<T>::makeEmpty(int n, int m){
  sparsity_ = CRSSparsity(n,m,false);
  std::vector<T>::clear();
}

template<class T>
void Matrix<T>::printScalar(std::ostream &stream) const {
  if(numel()!=1) throw CasadiException("operator=: argument not scalar");
  stream << (*this)(0);
}
  
template<class T>
void Matrix<T>::printVector(std::ostream &stream) const {
  // print dimension and first element
  stream << "[" << size1() << "](";
  stream << (*this)(0);
  for(int i=1; i<size1(); i++){
    stream << "," << (*this)(i);
  }
  stream << ")";  
}

template<class T>
void Matrix<T>::printMatrix(std::ostream &stream) const{ 
 // print dimension
 stream << "[" << size1() << "," << size2() << "]";
 // print the first row
 stream << "((";
  stream << (*this)(0,0);
  for(int j=1; j<size2(); j++){
    stream << "," << (*this)(0,j);
 }
 stream << ")";
 // print the rest of the rows
 for(int i=1; i<size1(); i++){  
    stream << ",("; 
    stream << (*this)(i,0);
    for(int j=1; j<size2(); j++){
      stream << "," << (*this)(i,j);
    }
    stream << ")";
  }
  stream << ")";  
}

template<class T>
void Matrix<T>::print(std::ostream &stream) const{
  if (empty())
    stream << "<empty expression>";
  else if (numel()==1)
    printScalar(stream);
  else if (size2()==1)
    printVector(stream);
  else
    printMatrix(stream);
}

template<class T>
const std::vector<int>& Matrix<T>::col() const{
  return sparsity_.col();
}

template<class T>
std::vector<int>& Matrix<T>::col(){
  return sparsity_.col();
}

template<class T>
const std::vector<int>& Matrix<T>::rowind() const{
  return sparsity_.rowind();
}

template<class T>
std::vector<int>& Matrix<T>::rowind(){
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
  std::vector<T>::reserve(nnz);
  sparsity_.reserve(nnz,nrow);
}

template<class T>
void Matrix<T>::resize(int n_, int m_){
  sparsity_.resize(n_,m_);
}

template<class T>
void Matrix<T>::clear(){
  sparsity_ = CRSSparsity(0,0,false);
  std::vector<T>::clear();
}

template<class T>
Matrix<T>::Matrix(const T &val){
  swap_on_copy_ = false;
  sparsity_ = CRSSparsity(1,1,true);
  std::vector<T>::resize(1,val);
}

template<class T>
Matrix<T>::Matrix(int n, int m, const std::vector<int>& col, const std::vector<int>& rowind){
  swap_on_copy_ = false;
  sparsity_ = CRSSparsity(n,m,col,rowind);
  std::vector<T>::resize(sparsity_.size());
}

template<class T>
Matrix<T>::Matrix(const CRSSparsity& sparsity){
  swap_on_copy_ = false;
  sparsity_ = sparsity;
  std::vector<T>::resize(sparsity_.size());
}


template<class T>
const T Matrix<T>::__getitem__(int i) const{
  return std::vector<T>::at(i);
}

template<class T>
const T Matrix<T>::__getitem__(const std::vector<int> &I) const{
  if(I.size()!=2) 
    throw CasADi::CasadiException("__getitem__: not 2D"); 
  return getElement(I[0],I[1]);
}

template<class T>
void Matrix<T>::__setitem__(int k, const T& el){ 
  std::vector<T>::at(k) = el;
}

template<class T>
void Matrix<T>::__setitem__(const std::vector<int> &I, const T&  el){ 
  if(I.size()!=2) 
    throw CasADi::CasadiException("__setitem__: not 2D"); 
  getElementRef(I[0],I[1]) = el;
}

template<class T>
void Matrix<T>::unary(T (*fcn)(const T&), const Matrix<T>& x){
  T temp = fcn(0);
  if(casadi_limits<T>::isZero(temp))
    makeEmpty(x.size1(),x.size2());
  else
    makeDense(x.size1(),x.size2(),temp);
    
  for(int i=0; i<size1(); ++i){ // loop over rows
    for(int el=x.rowind(i); el<x.rowind(i+1); ++el){
      int j = x.col(el);
      getElementRef(i,j) = fcn(x[el]);
    }
  }
}

template<class T>
void Matrix<T>::binary(T (*fcn)(const T&, const T&), const Matrix<T> &x, const Matrix<T> &y){
  if(x.scalar())
    if(y.scalar())
      *this = fcn(x(0),y(0));
    else
      scalar_matrix(fcn,x(0),y);
  else if(y.scalar())
    matrix_scalar(fcn,x,y(0));
  else
    matrix_matrix(fcn,x,y);
}

template<class T>
void Matrix<T>::scalar_matrix(T (*fcn)(const T&, const T&), const T& x, const Matrix<T>& y){
  T temp = fcn(x,0);
  if(casadi_limits<T>::isZero(temp))
    makeEmpty(y.size1(),y.size2());
  else
    makeDense(y.size1(),y.size2(),temp);
  
  for(int i=0; i<size1(); ++i){ // loop over rows
    for(int el=y.rowind(i); el<y.rowind(i+1); ++el){
      int j = y.col(el);
      getElementRef(i,j) = fcn(x,y[el]);
    }
  }
}

template<class T>
void Matrix<T>::matrix_scalar(T (*fcn)(const T&, const T&), const Matrix<T>& x, const T& y){
  T temp = fcn(0,y);
  if(casadi_limits<T>::isZero(temp))
    makeEmpty(x.size1(),x.size2());
  else
    makeDense(x.size1(),x.size2(),temp);
  
  for(int i=0; i<size1(); ++i){ // loop over rows
    for(int el=x.rowind(i); el<x.rowind(i+1); ++el){
      int j = x.col(el);
      getElementRef(i,j) = fcn(x[el],y);
    }
  }
}

template<class T>
void Matrix<T>::matrix_matrix(T (*fcn)(const T&, const T&), const Matrix<T>& x, const Matrix<T>& y){
if(x.size1() != y.size1() || x.size2() != y.size2()) throw CasadiException("matrix_matrix: dimension mismatch");
  T temp = fcn(0,0);
  if(casadi_limits<T>::isZero(temp))
    makeEmpty(x.size1(),x.size2());
  else
    makeDense(x.size1(),x.size2(),temp);
 
  for(int i=0; i<size1(); ++i){ // loop over rows
    int el1 = x.rowind(i);
    int el2 = y.rowind(i);
    int k1 = x.rowind(i+1);
    int k2 = y.rowind(i+1);
    while(el1 < k1 || el2 < k2){
      int j1 = (el1 < k1) ? x.col(el1) : numel() ;
      int j2 = (el2 < k2) ? y.col(el2) : numel() ;
      
      if(j1==j2)
        getElementRef(i,j1) = fcn(x[el1++],y[el2++]); 
      else if(j1>j2)
        getElementRef(i,j2) = fcn(0,y[el2++]);
      else
        getElementRef(i,j1) = fcn(x[el1++],0);
      }
    }
}

template<class T>
Matrix<T> Matrix<T>::operator-() const{
  Matrix<T> temp;
  temp.unary(casadi_operators<T>::neg,*this);
  return temp;
}

template<class T>
Matrix<T> Matrix<T>::operator+() const{
  return *this;
}

template<class T>
Matrix<T> Matrix<T>::operator+(const Matrix<T> &y) const{
  Matrix<T> r;
  r.binary(casadi_operators<T>::add,*this,y);
  return r;
}

template<class T>
Matrix<T> Matrix<T>::operator-(const Matrix<T> &y) const{
  Matrix<T> r;
  r.binary(casadi_operators<T>::sub,*this,y);
  return r;
}

template<class T>
Matrix<T> Matrix<T>::operator*(const Matrix<T> &y) const{
  Matrix<T> r;
  r.binary(casadi_operators<T>::mul,*this,y);
  return r;
}

template<class T>
Matrix<T> Matrix<T>::operator/(const Matrix<T> &y) const{
  Matrix<T> r;
  r.binary(casadi_operators<T>::div,*this,y);
  return r;
}

template<class T>
Matrix<T>& Matrix<T>::operator+=(const Matrix<T> &y){
  Matrix<T> x = *this;
  binary(casadi_operators<T>::add,x,y);
  return *this;
}

template<class T>
Matrix<T>& Matrix<T>::operator-=(const Matrix<T> &y){
  Matrix<T> x = *this;
  binary(casadi_operators<T>::sub,x,y);
  return *this;
}

template<class T>
Matrix<T>& Matrix<T>::operator*=(const Matrix<T> &y){
  Matrix<T> x = *this;
  binary(casadi_operators<T>::mul,x,y);
  return *this;
}

template<class T>
Matrix<T>& Matrix<T>::operator/=(const Matrix<T> &y){
  Matrix<T> x = *this;
  binary(casadi_operators<T>::div,x,y);
  return *this;
}

template<class T>
const CRSSparsity& Matrix<T>::sparsity() const{
  return sparsity_;
}

template<class T>
void Matrix<T>::getSparseSym(T *res) const{
  // copy to the result vector
  int nz = 0;
  for(int row=0; row<size1(); ++row)
  {
    // Loop over the elements in the row
    for(int el=rowind(row); el<rowind(row+1); ++el){ // loop over the non-zero elements
      if(col(el) > row) break; // break inner loop (only lower triangular part is used)
      res[nz] = (*this)[el];
      nz++;
    }
  }
}

template<class T>
void Matrix<T>::getTimesVector(const T *v, T *res) const{
  // copy the result
  for(int i=0; i<size1(); ++i){ // loop over rows
    res[i] = 0;
    for(int el=rowind(i); el<rowind(i+1); ++el){ // loop over the non-zero elements
      int j=col(el);  // column

      // Multiply with the vector
      res[i] += v[j]*(*this)[el];
    }
  }
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
      res[s + ldres*j] = (*this)[el];
    }
  }
}

template<class T>
void Matrix<T>::set(T val, Sparsity sp){
  assertNumEl(1);
  assertNNZ(1,sp);
  (*this)[0] = val;
}
    
template<class T>
void Matrix<T>::get(T& val, Sparsity sp) const{
  assertNumEl(1);
  assertNNZ(1,sp);
  val = (*this)[0];
}

template<class T>
void Matrix<T>::set(const std::vector<T>& val, Sparsity sp){
  assertNNZ(val.size(),sp);
  set(&val[0],sp);
}

template<class T>
void Matrix<T>::get(std::vector<T>& val, Sparsity sp) const{
  assertNNZ(val.size(),sp);
  get(&val[0],sp);
}

template<class T>
void Matrix<T>::set(const T* val, Sparsity sp){
  std::vector<T> &v = *this;
  if(sp==SPARSE || (sp==DENSE && numel()==std::vector<T>::size())){
    copy(val,val+v.size(),v.begin());
  } else if(sp==DENSE){
    for(int i=0; i<size1(); ++i) // loop over rows
      for(int el=rowind(i); el<rowind(i+1); ++el){ // loop over the non-zero elements
        // column
        int j=col(el);
        
        // Set the element
        v[el] = val[i*size2()+j];
    }
  } else {
    throw CasadiException("Matrix<T>::set: not SPARSE or DENSE");
  }
}

template<class T>
void Matrix<T>::get(T* val, Sparsity sp) const{
  const std::vector<T> &v = *this;
  if(sp==SPARSE || (sp==DENSE && numel()==v.size())){
    copy(v.begin(),v.end(),val);
  } else if(sp==DENSE){
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
  } else {
    throw CasadiException("Matrix<T>::get: not SPARSE or DENSE");
  }
}

template<class T>
void Matrix<T>::assertNNZ(int sz, Sparsity sp) const{
  int nnz_correct = -1;
  switch(sp){
    case SPARSE:
      nnz_correct = std::vector<T>::size();
      break;
    case DENSE:
      nnz_correct = numel();
      break;
    default:
      throw CasadiException("Matrix<T>::assertNNZ: unknown sparsity");
  }

  if(nnz_correct!=sz){
    std::stringstream ss;
    ss << "Matrix<T>::assertNNZ: wrong number of elements (" << sz << "), but should be " << nnz_correct << std::flush;
    throw CasadiException(ss.str());
  }
}

template<class T>
void Matrix<T>::assertNumEl(int sz) const{
  if(numel()!=sz){
    std::stringstream ss;
    ss << "Matrix<T>::assertNumEl: wrong number of elements (" << sz << " ), but should be " << numel();
    throw CasadiException(ss.str());
  }
}









#endif // SWIG

} // namespace CasADi


#endif // MATRIX_HPP

