/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A minimalistic computer algebra system with automatic differentiation 
 *    and framework for dynamic optimization.
 *    Copyright (C) 2010 by Joel Andersson, Moritz Diehl et al., K.U.Leuven. All rights reserved.
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
#include "casadi_exception.hpp"
#include "printable_object.hpp"
#include "element.hpp"

namespace CasADi{

  
/** \brief General sparse matrix class
  \author Joel Andersson 
  \date 2010	
*/
template<class T>
class Matrix : public std::vector<T>, public PrintableObject{

  public:
    /** \brief  constructors */
    /// empty 0-by-0 matrix constructor
    Matrix();                               // 
    /// empty n-by-m matrix constructor
    Matrix(int n, int m);                   
    /// dense n-by-m matrix filled with val constructor
    Matrix(int n, int m, const T& val);    

    /** \brief  This constructor enables implicit type conversion from a scalar type */
    template<typename A>
    Matrix(const A &val){
      makeEmpty(1,1);
      getElementRef()=val;
    }

    /** \brief  Create an expression from an stl vector  */
    template<typename A>
    Matrix(const std::vector<A>& x){
      makeDense(x.size(),1,1);
      copy(x.begin(),x.end(),std::vector<T>::begin());
    }

    /** \brief  Create a non-vector expression from an stl vector */
    template<typename A>
    Matrix(const std::vector<A>& x,  int n, int m){
      if(x.size() != n*m) throw CasadiException("Matrix::Matrix(const std::vector<double>& x,  int n, int m): dimension mismatch");
      makeDense(n,m,1);
      copy(x.begin(),x.end(),std::vector<T>::begin());
    }
    
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

    /// get an element
    const T getElement(int i=0, int j=0) const;
    
    /// get a reference to an element
    T& getElementRef(int i=0, int j=0);
  
    /// Access an element 
    Element<Matrix<T>,T> operator()(int i, int j=0){ return Element<Matrix<T>,T>(*this,i,j); }

    /// Const access an element 
    const T operator()(int i, int j=0) const{ return getElement(i,j); }

    /** \brief  Make the matrix an dense n-by-m matrix */
    void makeDense(int n, int m, const T& val);

    /** \brief  Make the matrix an empty n-by-m matrix */
    void makeEmpty(int n, int m);

    //@{
    /// Printing
    virtual void print(std::ostream &stream=std::cout) const; // print default style
    void printScalar(std::ostream &stream=std::cout) const; // print scalar
    void printVector(std::ostream &stream=std::cout) const; // print vector-style
    void printMatrix(std::ostream &stream=std::cout) const; // print matrix-style
    //@}

    // Get the sparsity pattern
    const std::vector<int>& getCol() const;
    const std::vector<int>& getRowInd() const;
    
    // make protected?
  public:
    // For efficiency, the following elements can be moved to a reference counted object, in which case the functions that actually change the sparsity or size have to be implemented a bit more carefully. Changes should only be permitted if only one instance exists, otherwise a deep copy must be made. This is probably simular to STL's std::string implementation.
    
    /** \brief  Elements of the matrix in a compressed row storage (CRS) format */
    std::vector<int> col;          // vector of length nnz containing the columns for all the indices of the non-zero elements
    std::vector<int> rowind;       // vector of length n+1 containing the index of the last non-zero element up till each row 
    int nrow;
    int ncol;

};

// Implementations
template<class T>
const T Matrix<T>::getElement(int i, int j) const{
  if(i >= size1() || j>=size2()) throw CasadiException("Matrix::getElement: out of bounds");
  for(int ind=rowind.at(i); ind<rowind.at(i+1); ++ind){
    if(col[ind] == j)
      return std::vector<T>::at(ind);     // quick return if element exists
    else if(col[ind] > j)
      break;                // break at the place where the element should be added
  }
  return 0;
}

template<class T>
T& Matrix<T>::getElementRef(int i, int j){
  if(i >= size1() || j>=size2()) throw CasadiException("Matrix::getElementRef: out of bounds");

  // go to the place where the element should be
  int ind;
  for(ind=rowind[i]; ind<rowind[i+1]; ++ind){ // better: loop from the back to the front
    if(col[ind] == j){
      return std::vector<T>::at(ind); // element exists
    } else if(col[ind] > j)
      break;                // break at the place where the element should be added
  }
  
  // insert the element
  std::vector<T>::insert(std::vector<T>::begin()+ind,0);
  col.insert(col.begin()+ind,j);

  for(int row=i+1; row<size1()+1; ++row)
    rowind[row]++;
    
  return std::vector<T>::at(ind);
}

template<class T>
int Matrix<T>::size1() const{
  return nrow;
}

template<class T>
int Matrix<T>::size2() const{
  return ncol;
}

template<class T>
int Matrix<T>::numel() const{
  return size1()*size2();
}

        
template<class T>
void Matrix<T>::makeDense(int n, int m, const T& val){
  nrow = n;
  ncol = m;
  std::vector<T>::clear();
  std::vector<T>::resize(n*m, val);
  col.resize(n*m);
  rowind.resize(n+1);
  
  int el =0;
  for(int i=0; i<size1(); ++i){
    rowind[i] = el;
    for(int j=0; j<size2(); ++j){
      col[el] = j;
      el++;
    }
  }
  rowind.back() = el;
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
  makeEmpty(0,0);
}

template<class T>
Matrix<T>::Matrix(int n, int m){
  makeEmpty(n,m);
}

template<class T>
Matrix<T>::Matrix(int n, int m, const T& val){
  makeDense(n,m,val);
}

template<class T>
void Matrix<T>::makeEmpty(int n, int m){
  nrow = n;
  ncol = m;
  std::vector<T>::clear();
  col.clear();
  rowind.clear();
  rowind.resize(n+1, 0);
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
const std::vector<int>& Matrix<T>::getCol() const{
  return col;
}

template<class T>
const std::vector<int>& Matrix<T>::getRowInd() const{
  return rowind;
}



} // namespace CasADi


#endif // MATRIX_HPP

