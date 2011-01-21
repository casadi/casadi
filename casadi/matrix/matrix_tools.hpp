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

#ifndef SWIG
/// Matrix product of two matrices - not available in Python since prod in numpy means elementwise multiplication
template<class T>
Matrix<T> prod(const Matrix<T> &x, const Matrix<T> &y);
#endif // SWIG

/// Product of two matrices - Python naming
template<class T>
Matrix<T> dot(const Matrix<T> &x, const Matrix<T> &y);

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

template<class T>
Matrix<T> getRow( const Matrix<T> &expr, int i, int ni=1, int ki=1);

/** \brief  A(i:ki:i+ni,:) = expr */
template<class T>
void setRow(const Matrix<T>& expr, Matrix<T> &res, int i, int ni=1, int ki=1);

/// ... = A(:,j:kj:j+nj)
template<class T>
void getColumn(Matrix<T> &res, const Matrix<T> &expr, int j, int nj=1, int kj=1);

template<class T>
Matrix<T> getColumn( const Matrix<T> &expr, int j, int nj=1, int kj=1);

/** \brief  A(:,j:kj:j+nj) = expr */
template<class T>
void setColumn(const Matrix<T>& expr, Matrix<T> &res, int j, int nj=1, int kj=1);

/** \brief  check if the matrix has certain properties */
template<class T>
bool isConstant(const Matrix<T>& ex);

template<class T>
bool isDense(const Matrix<T>& ex);

template<class T>
bool isEmpty(const Matrix<T>& ex);

template<class T>
bool isInteger(const Matrix<T>& ex);

template<class T>
bool isScalar(const Matrix<T>& ex);

template<class T>
bool isVector(const Matrix<T>& ex);

/** \brief  Check if a matrix is lower triangular (complexity ~ A.size1()) */
template<class T>
bool isTril(const Matrix<T> &A);

/** \brief  Check if a matrix is upper triangular (complexity ~ A.size1()) */
template<class T>
bool isTriu(const Matrix<T> &A);

template<class T>
T det(const Matrix<T>& a);

template<class T>
T getMinor(const Matrix<T> &x, int i, int j);

template<class T>
T cofactor(const Matrix<T> &x, int i, int j);

template<class T>
Matrix<T> adj(const Matrix<T>& a);

template<class T>
Matrix<T> inv(const Matrix<T>& a);

template<class T>
Matrix<T> reshape(const Matrix<T>& a, int n, int m);

/** \brief  make a vector
  Reshapes/flattens the Matrix<T> such that the shape becomes (expr.numel(),1).
  Columns are stacked on top of each other.
 */
template<class T>
Matrix<T> vec(const Matrix<T>& a);

template<class T>
Matrix<T> vertcat(const std::vector<Matrix<T> > &v);

template<class T>
Matrix<T> horzcat(const std::vector<Matrix<T> > &v);

#ifndef SWIG
template<class T>
Matrix<T> vertcat(const Matrix<T> &x, const Matrix<T> &y);

template<class T>
Matrix<T> horzcat(const Matrix<T> &x, const Matrix<T> &y);
#endif // SWIG

/** \brief Inner product of two vectors
        Equals
        \code
        trans(x)*y
        \endcode
        with x and y vectors
*/
template<class T>
Matrix<T> inner_prod(const Matrix<T> &x, const Matrix<T> &y); // inner product

/** \brief Outer product of two vectors
        Equals
        \code
        x*trans(y)
        \endcode
        with x and y vectors
*/
template<class T>
Matrix<T> outer_prod(const Matrix<T> &x, const Matrix<T> &y);

/** \brief  QR factorization using the modified Gram-Schmidt algorithm 
More stable than the classical Gram-Schmidt, but may break down if the columns of A are nearly linearly dependent
See J. Demmel: Applied Numerical Linear Algebra (algorithm 3.1.) */
template<class T>
void qr(const Matrix<T>& A, Matrix<T>& Q, Matrix<T> &R);

template<class T>
std::vector<Matrix<T> > qr(const Matrix<T>& A);

/** \brief  Solve a system of equations: A*x = b */
template<class T>
Matrix<T> solve(const Matrix<T>& A, const Matrix<T>& b);

/** \brief  Matlab's linspace function */
template<class T>
Matrix<T> linspace(const Matrix<T> &a, const Matrix<T> &b, int nsteps);

template<class T>
bool isZero(const Matrix<T>& ex);

template<class T>
int nnz(const Matrix<T>& ex);

template<class T>
int nnz_sym(const Matrix<T>& ex);

template<class T>
bool isEqual(const Matrix<T>& ex1,const Matrix<T> &ex2);

/** \brief  z += x*y: with x and y given (currently only if dense) */
template<class T>
void matrix_matrix_mult(const Matrix<T>& x, const Matrix<T>& y, Matrix<T>& z);

/** \brief  z += x*trans(y): with x and y and (currently only if dense) */
template<class T>
void matrix_matrix_trans_mult(const Matrix<T>& x, const Matrix<T>& y_trans, Matrix<T>& z);

/** \brief  z += trans(x)*y with x and y givem (currently only if dense)*/
template<class T>
void matrix_trans_matrix_mult(const Matrix<T>& x_trans, const Matrix<T>& y, Matrix<T>& z);

/// Make the vector 1-norm of an Matrix<T>
template<class T>
Matrix<T> norm_1(const Matrix<T>& x);
  
/// Make the vector 2-norm of an Matrix<T>
template<class T>
Matrix<T> norm_2(const Matrix<T>& x);

/// Return summation of all elements
template<class T>
T sum_all(const Matrix<T> &x); 

/** \brief Return summation of elements along specific axis
    \param axis either 0 or 1 */
template<class T>
Matrix<T> sum(const Matrix<T> &x, int axis=0);  

} // namespace CasADi

#ifndef SWIG

// Implementations of the functions in standard namespace

#include <iterator>

namespace CasADi{
// Implementations

template<class T>
Matrix<T> trans(const Matrix<T> &x){
  // quick return if empty or scalar
  if(x.empty() || x.scalar()) return x;

  // Create the new sparsity pattern and the mapping
  std::vector<int> mapping;
  CRSSparsity sparsity = x.sparsity().transpose(mapping);

  // create the return matrix
  Matrix<T> ret(sparsity);
  
  // Copy the content
  for(int i=0; i<mapping.size(); ++i)
    ret[i] = x[mapping[i]];
  
  return ret;
}

template<class T>
Matrix<T> dot(const Matrix<T> &x, const Matrix<T> &y){
  return prod<T>(x,y);
}

template<class T>
Matrix<T> prod(const Matrix<T> &x, const Matrix<T> &y){
  casadi_assert_message(x.size2() == y.size1(),"prod: dimension mismatch");
  
  // Find the mapping corresponding to the transpose of y (no need to form the transpose explicitly)
  std::vector<int> y_trans_map;
  CRSSparsity y_trans_sparsity = y.sparsity().transpose(y_trans_map);

  // Create the return object
  Matrix<T> ret(x.size1(),y.size2());
  std::vector<int>& rowind = ret.sparsityRef().rowind();

#ifdef LISTS_IN_PROD
  // ALTERNATIVE 1: GROWING LISTS, THEN COPY: BETTER FOR SYMBOLIC TYPES WHEN COPYING IS EXPENSIVE?
  std::list<int> col;
  std::list<T> val;
  std::vector<int>& ret_col = ret.col();
  std::vector<T>& ret_val = ret;
#else // LISTS_IN_PROD
  // ALTERNATIVE 2: GROWING VECTORS: APPEARS FASTER
  std::vector<int>& col = ret.sparsityRef().col();
  std::vector<T>& val = ret;

#ifdef ESTIMATE_NNZ_IN_PROD
  // Sparsity densitity in x and y
  double x_dens = x.size()/double(x.numel());
  double y_dens = y.size()/double(y.numel());
  
  // Probability that a product of an element of x and an element of y is zero
  double zero_prob = (1-x_dens)*(1-y_dens);
  
  // Probability that all of the products that make up a matrix entry is zero
  double element_zero_prob = std::pow(zero_prob,x.size2());

  // Make sure between zero and one (you never know how crude pow might be)
  element_zero_prob = fmin(fmax(element_zero_prob,0.0),1.0);

  // We obtain the estimated number of non-zeros, with 20 % safety margin
  int nnz_est = int(1.20*ret.numel()*(1-element_zero_prob));
  
  // Make sure not larger than the number of elements
  if(nnz_est>ret.numel())
    nnz_est = ret.numel();
    
  // Pass the guess to the vectors
  col.reserve(nnz_est);
  val.reserve(nnz_est);
#endif // ESTIMATE_NNZ_IN_PROD
#endif // LISTS_IN_PROD

  // Direct access to the arrays
  const std::vector<int> &x_col = x.col();
  const std::vector<int> &y_row = y_trans_sparsity.col();
  const std::vector<int> &x_rowind = x.rowind();
  const std::vector<int> &y_colind = y_trans_sparsity.rowind();
  
  // loop over the row of the resulting matrix)
  for(int i=0; i<x.size1(); ++i){
    for(int j=0; j<y.size2(); ++j){ // loop over the column of the resulting matrix
      int el1 = x_rowind[i];
      int el2 = y_colind[j];
      T d = 0; // the entry of the matrix to be calculated
      bool added = false;
      while(el1 < x_rowind[i+1] && el2 < y_colind[j+1]){ // loop over non-zero elements
        int j1 = x_col[el1];
        int i2 = y_row[el2];      
        if(j1==i2){
          added = true;
          d += x[el1++] * y[y_trans_map[el2++]];
        } else if(j1<i2) {
          el1++;
        } else {
          el2++;
        }
      }
      if(added){
        col.push_back(j);
        val.push_back(d);
      }
    }
    rowind[i+1] = col.size();
  }
  
#ifdef LISTS_IN_PROD
  // Save the column indices to the return matrix
  ret_col.resize(col.size());
  copy(col.begin(),col.end(),ret_col.begin());
  
  // Save the non-zero entries to the return matrix
  ret_val.resize(val.size());
  copy(val.begin(),val.end(),ret_val.begin());
#endif // LISTS_IN_PROD
  
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
  casadi_assert_message(expr.size2() == add.size2(),"append: dimensions do not match");

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
  casadi_assert_message(ki==1 && kj==1, "getSub: ki!=1 and kj!=1 not implemented");
  casadi_assert_message(i+ni <= expr.size1() && j+nj <= expr.size2(),"getSub: dimension mismatch");
  res = Matrix<T>(ni,nj);
  for(int r=0; r<ni; ++r)
    for(int c=0; c<nj; ++c)
      if(!casadi_limits<T>::isZero(expr(i+r,j+c)))
        res(r,c) = expr(i+r,j+c);
}

template<class T>
void setSub(const Matrix<T> &expr, Matrix<T> &res, int i, int j){
  // TODO: Very inefficient, exploit sparsity
  for(int r=0; r<expr.size1(); ++r)
    for(int c=0; c<expr.size2(); ++c)
      if(!casadi_limits<T>::isZero(expr(r,c)))
        res(i+r,j+c) = expr(r,c);
}

template<class T>
void getRow(Matrix<T> &res, const Matrix<T> &expr, int i, int ni, int ki){
  // TODO: Very inefficient, exploit sparsity
  casadi_assert_message(i<expr.size1(),"dimension mismatch");
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
  // TODO: Very inefficient, exploit sparsity
  casadi_assert_message(j<expr.size2(),"dimension mismatch");
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
  // TODO: Very inefficient, exploit sparsity
  casadi_assert_message(i<res.size1(),"dimension mismatch");
  for(int j=0; j<res.size2(); ++j){
    if(!casadi_limits<T>::isZero(expr(0,j)))
      res(i,j) = expr(0,j);
    }
}

template<class T>
void setColumn(const Matrix<T>& expr, Matrix<T> &res, int j, int nj, int kj){
  // TODO: Very inefficient, exploit sparsity
  casadi_assert_message(j<res.size2(),"dimension mismatch");
  for(int i=0; i<res.size1(); ++i){
    if(!casadi_limits<T>::isZero(expr(i,0)))
      res(i,j) = expr(i,0);
    }
}

template<class T>
bool isScalar(const Matrix<T>& ex){
  return ex.size1()==1 && ex.size2()==1;
}

template<class T>
bool isVector(const Matrix<T>& ex){
  return ex.size2()==1;
}

template<class T>
bool isInteger(const Matrix<T>& ex){
  // loop over non-zero elements
  for(int k=0; k<ex.size(); ++k) 
    if(!casadi_limits<T>::isInteger(ex[k])) // if an element is not integer
      return false;
    
  // Integer if reached this point
  return true;
}

template<class T>
bool isConstant(const Matrix<T>& ex){
  // loop over non-zero elements
  for(int k=0; k<ex.size(); ++k) 
    if(!casadi_limits<T>::isConstant(ex[k])) // if an element is not constant
      return false;
    
  // Constant if we reach this point
  return true;
}

template<class T>
bool isDense(const Matrix<T>& ex){
  return ex.size() == ex.numel();
}

template<class T>
bool isEmpty(const Matrix<T>& ex){
  return ex.empty();
}

template<class T>
bool isTril(const Matrix<T> &A){
  // TODO: Move implementation to CRSSparsity as it does not depend on the matrix entries 
  // loop over rows
  for(int i=0; i<A.size1(); ++i){
    if(A.rowind(i) != A.rowind(i+1)){ // if there are any elements of the row
      // check column of the right-most element of the row
      int col = A.col(A.rowind(i+1)-1);

      // not lower triangular if col>i
      if(col>i) return false;
    }
  }
  // all rows ok
  return true;
}

template<class T>
bool isTriu(const Matrix<T> &A){
  // TODO: Move implementation to CRSSparsity as it does not depend on the matrix entries 
  // loop over rows
  for(int i=0; i<A.size1(); ++i){
    if(A.rowind(i) != A.rowind(i+1)){ // if there are any elements of the row
      // check column of the left-most element of the row
      int col = A.col(A.rowind(i));

      // not lower triangular if col>i
      if(col<i) return false;
    }
  }
  // all rows ok
  return true;
}

template<class T>
T det(const Matrix<T>& a){
  int n = a.size1();
  casadi_assert_message(n == a.size2(),"matrix must be square");

  // Trivial return if scalar
  if(isScalar(a)) return a(0);

  // Return expression
  T ret = 0;

  // We expand the matrix along the first column
  for(int i=0; i<n; ++i){

    // Sum up the cofactors
    ret += a(i,0)*cofactor(a,i,0);

  }
  return ret;
}

template<class T>
T getMinor(const Matrix<T> &x, int i, int j){
  int n = x.size1();
  casadi_assert_message(n == x.size2(), "getMinor: matrix must be square");

  // Trivial return if scalar
  if(n==1) return 1;

  // Remove row i and column j
  Matrix<T> M(n-1,n-1);

   for(int i1=0; i1<n; ++i1)
       for(int j1=0; j1<n; ++j1){
           if(i1 == i || j1 == j)
              continue;

            int i2 = (i1<i)?i1:i1-1;
            int j2 = (j1<j)?j1:j1-1;
    
            M(i2,j2) = x(i1,j1);
       }
  return det(M);
}

template<class T>
T cofactor(const Matrix<T> &x, int i, int j){

    // Calculate the i,j minor
    T minor_ij = getMinor(x,i,j);

    // Calculate the cofactor
    int sign_i = 1-2*((i+j) % 2);

    return sign_i * minor_ij;
}

template<class T>
Matrix<T> adj(const Matrix<T>& a){
  int n = a.size1();
  casadi_assert_message(n == a.size2(),"adj: matrix must be square");

  // Cofactor matrix
  Matrix<T> C(n,n);
  for(int i=0; i<n; ++i)
    for(int j=0; j<n; ++j)
      C(i,j) = cofactor(a,i,j);
  
  return trans(C);
}

template<class T>
Matrix<T> inv(const Matrix<T>& a){
  // laplace formula
  return adj(a)/det(a);
}

template<class T>
Matrix<T> reshape(const Matrix<T>& a, int n, int m){
  casadi_assert_message(a.numel() == n*m, "resize: number of elements must remain the same");

  Matrix<T> ret(n,m);
  for(int i=0; i<a.size1(); ++i){
    for(int el=a.rowind(i); el<a.rowind(i+1); ++el){
      int j = a.col(el);
      int k = j+i*a.size1();
      ret(k/m,k%m) = a[el];
    }
  }
  return ret;
}

template<class T>
Matrix<T> vec(const Matrix<T>& a){
  Matrix<T> ret = reshape(trans(a),a.numel(),1);
  return ret;
}

template<class T>
Matrix<T> vertcat(const std::vector<Matrix<T> > &v){
  Matrix<T> ret;
  for(int i=0; i<v.size(); ++i)
    append(ret,v[i]);
  return ret;
}

template<class T>
Matrix<T> horzcat(const std::vector<Matrix<T> > &v){
  Matrix<T> ret;
  for(int i=0; i<v.size(); ++i)
    append(ret,trans(v[i]));
  return trans(ret);  
}

template<class T>
Matrix<T> vertcat(const Matrix<T> &x, const Matrix<T> &y){
  Matrix<T> xy = x;
  append(xy,y);
  return xy;
}

template<class T>
Matrix<T> horzcat(const Matrix<T> &x, const Matrix<T> &y){
  return trans(vertcat(trans(x),trans(y)));
}

template<class T>
Matrix<T> inner_prod(const Matrix<T> &x, const Matrix<T> &y){
  casadi_assert_message(x.vector() && y.vector(), "inner_prod: arguments must be vectors");
  return prod(trans(x),y);
}

template<class T>
Matrix<T> outer_prod(const Matrix<T> &x, const Matrix<T> &y){
  casadi_assert_message(x.vector() && y.vector(), "outer_prod: arguments must be vectors");
  return prod(x,trans(y));  
}

template<class T>
T sum_all(const Matrix<T> &x) {
  // Sum non-zero elements
  T res=0;
  for(int k=0; k<x.size(); k++){
    res += x[k];
  }
  return res;
}

template<class T>
Matrix<T> sum(const Matrix<T> &x, int axis) {
  casadi_assert_message(axis==0 || axis==1,"axis argument should be zero or one");
  if (axis==1){
    return trans(sum(trans(x),0));
  } else {
    Matrix<T> res(1,x.size2());
    const std::vector<int> &col = x.col();
    for(int k=0; k<col.size(); k++){
      res(0,col[k]) += x[k];
    }
    return res;
  }
}

template<class T>
Matrix<T> norm_1(const Matrix<T>& x){
  return sum_all(fabs(x));
}

template<class T>
Matrix<T> norm_2(const Matrix<T>& x){
  return sqrt(inner_prod(x,x));
}

template<class T>
std::vector<Matrix<T> > qr(const Matrix<T>& A){
  // Create return value
  std::vector<Matrix<T> > QR(2);
  Matrix<T> &Q = QR[0];
  Matrix<T> &R = QR[1];
  
  qr(A,Q,R);
  return QR;
}

template<class T>
void qr(const Matrix<T>& A, Matrix<T>& Q, Matrix<T> &R){
  // The following algorithm is taken from J. Demmel: Applied Numerical Linear Algebra (algorithm 3.1.)
  int m = A.size1();
  int n = A.size2();
  casadi_assert_message(m>=n, "qr: fewer rows than columns");

  // Transpose of A
  Matrix<T> AT = trans(A);

  // Transposes of the output matrices
  Matrix<T> QT, RT;

  // compute Q and R column by column
  for(int i=0; i<n; ++i){
    // Initialize qi to be the i-th column of A
    Matrix<T> ai = getRow(AT,i);
    Matrix<T> qi = ai;
    // The i-th column of R
    Matrix<T> ri(1,n);
  
    // subtract the projection of of qi in the previous directions from ai
    for(int j=0; j<i; ++j){
      // Get the j-th column of Q
      Matrix<T> qj = getRow(QT,j);

      ri(0,j) = prod(qj,trans(qi))[0]; // Modified Gram-Schmidt
      // ri[j] = inner_prod(qj,ai); // Classical Gram-Schmidt

      // Remove projection in direction j
      qi -= ri(0,j) * qj;
    }

    // Normalize qi
    ri(0,i) = norm_2(trans(qi))[0];
    qi /= ri(0,i);

    // Update RT and QT
    append(QT,qi);
    append(RT,ri);
  }

  // Save to output
  Q = trans(QT);
  R = trans(RT);
}

template<class T>
Matrix<T> solve(const Matrix<T>& A, const Matrix<T>& b){
  // check dimensions
  casadi_assert_message(A.size1() == b.size1(),"solve: dimension mismatch");
  casadi_assert_message(A.size1() == A.size2(),"solve: A not square");
  if(isTril(A)){
    // forward substiution if lower triangular
    Matrix<T> x = b;
    for(int i=0; i<A.size1(); ++i){ // loop over rows
      for(int k=0; k<b.size2(); ++k){ // for every right hand side
        for(int j=0; j<i; ++j){
          x(i,k) -= A(i,j)*x(j,k);
        }
        x(i,k) /= A(i,i);
      }
    }
    return x;
  } else if(isTriu(A)){
    // backward substiution if upper triangular
    Matrix<T> x = b;
    for(int i=A.size1()-1; i>=0; --i){ // loop over rows from the back
      for(int k=0; k<b.size2(); ++k){ // for every right hand side
        for(int j=A.size1()-1; j>i; --j){
          x(i,k) -= A(i,j)*x(j,k);
        }
        x(i,k) /= A(i,i);
      }
    }
    return x;
  } else {
    // Make a QR factorization
    Matrix<T> Q,R;
    qr(A,Q,R);

    // Solve the factorized system
    return solve(R,prod(trans(Q),b));
  }
}

template<class T>
Matrix<T> linspace(const Matrix<T> &a_, const Matrix<T> &b_, int nsteps){
  casadi_assert_message(isScalar(a_) && isScalar(b_), "linspace: a and b must be scalar");
  T a = a_(0);
  T b = b_(0);
  Matrix<T> ret(nsteps,1);
  ret(0) = a;
  T step = (b-a)/(nsteps-1);

  for(int i=1; i<nsteps-1; ++i)
    ret(i) = ret(i-1) + step;
  
  ret(nsteps-1) = b;
  return ret;
}

template<class T>
bool isZero(const Matrix<T>& ex) {  

  // loop over (potentially) non-zero elements
  for(int el=0; el<ex.size(); ++el)
    if(!casadi_limits<T>::isZero(ex[el]))
      return false;
  
  return true;
}

template<class T>
int nnz(const Matrix<T>& ex) {
  return ex.size();
}

template<class T>
int nnz_sym(const Matrix<T>& ex) {
  int nz = 0; // number of non-zeros  
  for(int row=0; row<ex.size1(); ++row)
  {
    // Loop over the elements in the row
    for(int el=ex.rowind(row); el<ex.rowind(row+1); ++el){ // loop over the non-zero elements
      if(ex.col(el) > row) break; // break inner loop (only lower triangular part is used)
      nz++;
    }
  }
  return nz;
}

template<class T>
bool isEqual(const Matrix<T>& ex1,const Matrix<T> &ex2){
  if ((nnz(ex1)!=0 || nnz(ex2)!=0) && (ex1.size1()!=ex2.size1() || ex1.size2()!=ex2.size2())) return false;
  Matrix<T> difference = ex1 - ex2;  
  return isZero(difference);
}

template<class T>
void matrix_matrix_mult(const Matrix<T>& x, const Matrix<T>& y, Matrix<T>& z){
  //dimensions
  int ni = x.size1();
  int nk = x.size2();
  int nj = y.size2();
  if(x.dense() && y.dense() && z.dense()){
    for(int j=0; j<nj; ++j)
      for(int k=0; k<nk; ++k)
        for(int i=0; i<ni; ++i)
          z[i + ni*j] += x[i + ni*k]*y[k + nk*j];
  } else {
    throw CasadiException("matrix_matrix_mult: not dense");
  }
}

template<class T>
void matrix_matrix_trans_mult(const Matrix<T>& x, const Matrix<T>& y_trans, Matrix<T>& z){
  //dimensions
  int ni = z.size1();
  int nk = z.size2();
  int nj = y_trans.size2();
  if(x.dense() && y_trans.dense() && z.dense()){
    for(int j=0; j<nj; ++j)
      for(int k=0; k<nk; ++k)
        for(int i=0; i<ni; ++i)
          z[i + ni*k] += x[i + ni*j]*y_trans[k + nk*j];
  } else {
    throw CasadiException("matrix_matrix_trans_mult: not dense");
  }
}

template<class T>
void matrix_trans_matrix_mult(const Matrix<T>& x_trans, const Matrix<T>& y, Matrix<T>& z){
  //dimensions
  int ni = x_trans.size1();
  int nk = x_trans.size2();
  int nj = z.size2();
  if(x_trans.dense() && y.dense() && z.dense()){
    for(int j=0; j<nj; ++j)
      for(int k=0; k<nk; ++k)
        for(int i=0; i<ni; ++i)
          z[k + nk*j] += x_trans[i + ni*k]*y[i + ni*j];
  } else {
    throw CasadiException("matrix_trans_matrix_mult: not dense");
  }
}



} // namespace CasADi

#endif //SWIG




#ifdef SWIG

// map the template name to the instantiated name
#define MTT_INST(T,function_name) \
%template(function_name) CasADi::function_name<T>;

// Define template instanciations
#define MATRIX_TOOLS_TEMPLATES(T) \
MTT_INST(T,trans) \
MTT_INST(T,dot) \
MTT_INST(T,append) \
MTT_INST(T,getSub) \
MTT_INST(T,setSub) \
MTT_INST(T,getRow) \
MTT_INST(T,setRow) \
MTT_INST(T,getColumn) \
MTT_INST(T,setColumn) \
MTT_INST(T,isConstant) \
MTT_INST(T,isDense) \
MTT_INST(T,isEmpty) \
MTT_INST(T,isInteger) \
MTT_INST(T,isScalar) \
MTT_INST(T,isVector) \
MTT_INST(T,isTril) \
MTT_INST(T,isTriu) \
MTT_INST(T,det) \
MTT_INST(T,getMinor) \
MTT_INST(T,cofactor) \
MTT_INST(T,adj) \
MTT_INST(T,inv) \
MTT_INST(T,reshape) \
MTT_INST(T,vec) \
MTT_INST(T,horzcat) \
MTT_INST(T,vertcat) \
MTT_INST(T,inner_prod) \
MTT_INST(T,outer_prod) \
MTT_INST(T,norm_2) \
MTT_INST(T,qr) \
MTT_INST(T,solve) \
MTT_INST(T,linspace) \
MTT_INST(T,isZero) \
MTT_INST(T,nnz) \
MTT_INST(T,nnz_sym) \
MTT_INST(T,isEqual) \


#endif //SWIG



#endif // MATRIX_TOOLS_HPP
