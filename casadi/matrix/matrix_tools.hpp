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

#ifndef MATRIX_TOOLS_HPP
#define MATRIX_TOOLS_HPP

#include "matrix.hpp"

namespace CasADi{
  
/// Transpose of a matrix
template<class T>
Matrix<T> trans(const Matrix<T> &x);

/// Product of two matrices
template<class T>
Matrix<T> prod(const Matrix<T> &x, const Matrix<T> &y);

template<class T>
void append(Matrix<T>& expr, const Matrix<T>& add);

/// ... = A(i:ki:i+ni,j:kj:j+nj)
template<class T>
void getSub(Matrix<T> &res, const Matrix<T> &expr, int i, int j=0, int ni=1, int nj=1, int ki=1, int kj=1); 

/** \brief  A(i:ki:i+ni,j:kj:j+nj) = expr */
template<class T>
void setSub(const Matrix<T> &expr, Matrix<T> &res, int i, int j=0);

/// ... = A(i:ki:i+ni,:)
template<class T>
void getRow(Matrix<T> &res, const Matrix<T> &expr, int i, int ni=1, int ki=1);

/** \brief  A(i:ki:i+ni,:) = expr */
template<class T>
void setRow(const Matrix<T>& expr, Matrix<T> &res, int i, int ni=1, int ki=1);

/// ... = A(:,j:kj:j+nj)
template<class T>
void getColumn(Matrix<T> &res, const Matrix<T> &expr, int j, int nj=1, int kj=1);

/** \brief  A(:,j:kj:j+nj) = expr */
template<class T>
void setColumn(const Matrix<T>& expr, Matrix<T> &res, int j, int nj=1, int kj=1);

#ifndef SWIG
// Implementations
template<class T>
Matrix<T> trans(const Matrix<T> &x){
  // quick return if empty or scalar
  if(x.empty() || x.scalar()) return x;

  // We do matrix transpose by the (linear time) "bucket sort" algorithm
  std::vector<std::vector<int> > buckets(x.size2()); // one bucket for each column
  
  // Create a vector with the rows for each non-zero element
  std::vector<int> row(x.size());
  
  // Loop over the rows of the original matrix
  for(int i=0; i<x.size1(); ++i)
  {
    // Loop over the elements in the row
    for(int el=x.rowind(i); el<x.rowind(i+1); ++el){ // loop over the non-zero elements
      int j=x.col(el);  // column
      
     // put the element into the right bucket
     buckets[j].push_back(el);
     
     // save the row index
     row[el] = i;
    }
  }

  Matrix<T> ret(x.size2(),x.size1()); // create the return matrix

  // reserve space (to make the calculations quicker)
  ret.reserve(x.capacity());
  ret.col().reserve(x.col().size());
  ret.rowind().reserve(x.rowind().size());

  for(int j=0; j<x.size2(); ++j)   // loop over the columns
    for(int r=0; r<buckets[j].size(); ++r){ // loop over the bucket content
     int el =  buckets[j][r]; // the index of the non-zero element
     int i = row[el]; // the row of the element
     ret.getElementRef(j,i) = x[el]; // add the element
    }

    return ret;
}

template<class T>
Matrix<T> prod(const Matrix<T> &x, const Matrix<T> &y){
  if(x.size2() != y.size1()) throw CasadiException("prod: dimension mismatch");

  Matrix<T> ret(x.size1(),y.size2());
  Matrix<T> b = trans(y); // take the transpose of the second matrix: linear time operation

  for(int i=0; i<x.size1(); ++i) // loop over the row of the resulting matrix)
    for(int j=0; j<b.size1(); ++j){ // loop over the column of the resulting matrix
      int el1 = x.rowind(i);
      int el2 = b.rowind(j);
      while(el1 < x.rowind(i+1) && el2 < b.rowind(j+1)){ // loop over non-zero elements
        int j1 = x.col(el1);
        int i2 = b.col(el2);      
        if(j1==i2){
          T temp = x[el1++] * b[el2++];
          if(!casadi_limits<T>::isZero(temp))
            ret(i,j) += temp;
        } else if(j1<i2) {
          el1++;
        } else {
          el2++;
        }
      }
    }
    return ret;
}

template<class T>
void append(Matrix<T>& expr, const Matrix<T>& add){
  // Quick return if we are adding an empty expression
  if(add.empty()) return;

  // Likewise if expr is empty
  if(expr.empty()){
    expr=add;
    return;
  }

  // Check dimensions
  if(expr.size2() != add.size2()) throw CasadiException("append: dimensions do not match");

  // Resize the expression
  int oldn = expr.size1();
  int n    = expr.size1() + add.size1();  
  int m    = expr.size2();
  expr.resize(n,m);

  // Copy the lower expression to the end
  setSub(add, expr, oldn, 0);
}

template<class T>
void getSub(Matrix<T> &res, const Matrix<T> &expr, int i, int j, int ni, int nj, int ki, int kj){
  if(!(ki==1 && kj==1)) throw CasadiException("getSub: ki!=1 and kj!=1 not implemented");
  if(!(i+ni <= expr.size1() && j+nj <= expr.size2())) throw CasadiException("getSub: dimension mismatch");
  res = Matrix<T>(ni,nj);
  for(int r=0; r<ni; ++r)
    for(int c=0; c<nj; ++c)
      if(!casadi_limits<T>::isZero(expr(i+r,j+c)))
        res(r,c) = expr(i+r,j+c);
}

template<class T>
void setSub(const Matrix<T> &expr, Matrix<T> &res, int i, int j){
  for(int r=0; r<expr.size1(); ++r)
    for(int c=0; c<expr.size2(); ++c)
      if(!casadi_limits<T>::isZero(expr(r,c)))
        res(i+r,j+c) = expr(r,c);
}

template<class T>
void getRow(Matrix<T> &res, const Matrix<T> &expr, int i, int ni, int ki){
  if(!(i<expr.size1())) throw CasadiException("getRow: dimension mismatch");
  res = Matrix<T>(ni,expr.size2());
  for(int ii=0; ii<ni; ++ii)
    for(int j=0; j<expr.size2(); ++j){
      T temp = expr(i+ii,j);
      if(!casadi_limits<T>::isZero(temp))
        res(ii,j) = temp;
    }
}

template<class T>
Matrix<T> getRow( const Matrix<T> &expr, int i, int ni, int ki){
  Matrix<T> res(ni,expr.size2());
  getRow(res,expr,i,ni,ki);
  return res;
}

template<class T>
void getColumn(Matrix<T> &res, const Matrix<T> &expr, int j, int nj, int kj){
  if(!(j<expr.size2())) throw CasadiException("getColumn: dimension mismatch");
  res = Matrix<T>(expr.size1(),nj);
  for(int i=0; i<expr.size1(); ++i)
    for(int jj=0; jj<nj; ++jj){
      T temp = expr(i,j+jj);
      if(!casadi_limits<T>::isZero(temp))
        res(i,jj) = temp;
    }
}

template<class T>
Matrix<T> getColumn( const Matrix<T> &expr, int j, int nj, int kj){
  Matrix<T> res(expr.size1(),nj);
  getColumn(res,expr,j,nj,kj);
  return res;
}

template<class T>
void setRow(const Matrix<T>& expr, Matrix<T> &res, int i, int ni, int ki){
  if(!(i<res.size1())) throw CasadiException("setRow: dimension mismatch");
  for(int j=0; j<res.size2(); ++j){
    if(!casadi_limits<T>::isZero(expr(0,j)))
      res(i,j) = expr(0,j);
    }
}

template<class T>
void setColumn(const Matrix<T>& expr, Matrix<T> &res, int j, int nj, int kj){
  if(!(j<res.size2())) throw CasadiException("setColumn: dimension mismatch");
  for(int i=0; i<res.size1(); ++i){
    if(!casadi_limits<T>::isZero(expr(i,0)))
      res(i,j) = expr(i,0);
    }
}

# else // SWIG

// map the template name to the instantiated name
#define MTT_INST(T,function_name) \
%template(function_name) function_name<T>;

// Define template instanciations
#define MATRIX_TOOLS_TEMPLATES(T) \
MTT_INST(T,trans) \
MTT_INST(T,prod) \
MTT_INST(T,append) \
MTT_INST(T,getSub) \
MTT_INST(T,setSub) \
MTT_INST(T,getRow) \
MTT_INST(T,setRow) \
MTT_INST(T,getColumn) \
MTT_INST(T,setColumn) \

#endif //SWIG


} // namespace CasADi

#endif // MATRIX_TOOLS_HPP
