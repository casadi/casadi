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

// The following functions must be placed in the standard namespace so that the old ones are not shadowed when CasADi namespace is used
namespace std{

template<class T>
CasADi::Matrix<T> sin(const CasADi::Matrix<T>& x);

template<class T>
CasADi::Matrix<T> cos(const CasADi::Matrix<T>& x);

template<class T>
CasADi::Matrix<T> tan(const CasADi::Matrix<T>& x);

template<class T>
CasADi::Matrix<T> asin(const CasADi::Matrix<T>& x);

template<class T>
CasADi::Matrix<T> acos(const CasADi::Matrix<T>& x);

template<class T>
CasADi::Matrix<T> atan(const CasADi::Matrix<T>& x);

template<class T>
CasADi::Matrix<T> exp(const CasADi::Matrix<T>& x);

template<class T>
CasADi::Matrix<T> log(const CasADi::Matrix<T>& x);

template<class T>
CasADi::Matrix<T> sqrt(const CasADi::Matrix<T>& x);

template<class T>
CasADi::Matrix<T> floor(const CasADi::Matrix<T>& x);

template<class T>
CasADi::Matrix<T> ceil(const CasADi::Matrix<T>& x);

template<class T>
CasADi::Matrix<T> fabs(const CasADi::Matrix<T>& x);

template<class T>
CasADi::Matrix<T> fmin(const CasADi::Matrix<T>& x, const CasADi::Matrix<T>& y);

template<class T>
CasADi::Matrix<T> fmax(const CasADi::Matrix<T>& x, const CasADi::Matrix<T>& y);

template<class T>
CasADi::Matrix<T> pow(const CasADi::Matrix<T>& x, const CasADi::Matrix<T>& y);
  
} // namespace std

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

/// Get the 2-norm of a vector
template<class T>
Matrix<T> norm_2(const Matrix<T>& x);

/** \brief  QR factorization using the modified Gram-Schmidt algorithm 
More stable than the classical Gram-Schmidt, but may break down if the columns of A are nearly linearly dependent
See J. Demmel: Applied Numerical Linear Algebra (algorithm 3.1.) */
template<class T>
void qr(const Matrix<T>& A, Matrix<T>& Q, Matrix<T> &R);

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

} // namespace CasADi

#ifndef SWIG

// Implementations of the functions in standard namespace

#include <iterator>

namespace std{
  
#define UNOP_DEF(fname,opname) \
template<class T> \
CasADi::Matrix<T> fname(const CasADi::Matrix<T>& x){ \
  CasADi::Matrix<T> temp; \
  temp.unary(CasADi::casadi_operators<T>::opname,x); \
  return temp;\
}

#define BINOP_DEF(fname,opname) \
template<class T> \
CasADi::Matrix<T> fname(const CasADi::Matrix<T>& x, const CasADi::Matrix<T>& y){ \
  CasADi::Matrix<T> temp; \
  temp.binary(CasADi::casadi_operators<T>::opname,x,y); \
  return temp;\
}

UNOP_DEF(sin,sin)
UNOP_DEF(cos,cos)
UNOP_DEF(tan,tan)
UNOP_DEF(asin,asin)
UNOP_DEF(acos,acos)
UNOP_DEF(atan,atan)
UNOP_DEF(exp,exp)
UNOP_DEF(log,log)
UNOP_DEF(sqrt,sqrt)
UNOP_DEF(floor,floor)
UNOP_DEF(ceil,ceil)
UNOP_DEF(fabs,fabs)
UNOP_DEF(fmin,fmin)
UNOP_DEF(fmax,fmax)
UNOP_DEF(pow,pow)
  
} // namespace std

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
Matrix<T> prod(const Matrix<T> &x, const Matrix<T> &y){
  if(x.size2() != y.size1()) throw CasadiException("prod: dimension mismatch");
  
  // Find the mapping corresponding to the transpose of y (no need to form the transpose explicitly)
  std::vector<int> y_trans_map;
  CRSSparsity y_trans_sparsity = y.sparsity().transpose(y_trans_map);

  // Create the return object
  Matrix<T> ret(x.size1(),y.size2());
  std::vector<int>& rowind = ret.rowind();

#ifdef LISTS_IN_PROD
  // ALTERNATIVE 1: GROWING LISTS, THEN COPY: BETTER FOR SYMBOLIC TYPES WHEN COPYING IS EXPENSIVE?
  std::list<int> col;
  std::list<T> val;
  std::vector<int>& ret_col = ret.col();
  std::vector<T>& ret_val = ret;
#else // LISTS_IN_PROD
  // ALTERNATIVE 2: GROWING VECTORS: APPEARS FASTER
  std::vector<int>& col = ret.col();
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
  // TODO: Very inefficient, exploit sparsity
  for(int r=0; r<expr.size1(); ++r)
    for(int c=0; c<expr.size2(); ++c)
      if(!casadi_limits<T>::isZero(expr(r,c)))
        res(i+r,j+c) = expr(r,c);
}

template<class T>
void getRow(Matrix<T> &res, const Matrix<T> &expr, int i, int ni, int ki){
  // TODO: Very inefficient, exploit sparsity
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
  // TODO: Very inefficient, exploit sparsity
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
  // TODO: Very inefficient, exploit sparsity
  if(!(i<res.size1())) throw CasadiException("setRow: dimension mismatch");
  for(int j=0; j<res.size2(); ++j){
    if(!casadi_limits<T>::isZero(expr(0,j)))
      res(i,j) = expr(0,j);
    }
}

template<class T>
void setColumn(const Matrix<T>& expr, Matrix<T> &res, int j, int nj, int kj){
  // TODO: Very inefficient, exploit sparsity
  if(!(j<res.size2())) throw CasadiException("setColumn: dimension mismatch");
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
  if(n != a.size2()) throw CasadiException("det: matrix must be square");

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
  if(n != x.size2()) throw CasadiException("getMinor: matrix must be square");

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
  if(n != a.size2()) throw CasadiException("adj: matrix must be square");

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
  if(a.numel() != n*m)
    throw CasadiException("resize: number of elements must remain the same");
  
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
  if(!x.vector() || !y.vector()) throw CasadiException("inner_prod: arguments must be vectors");
  return prod(trans(x),y);
}

template<class T>
Matrix<T> outer_prod(const Matrix<T> &x, const Matrix<T> &y){
  if(!x.vector() || !y.vector()) throw CasadiException("outer_prod: arguments must be vectors");
  return prod(x,trans(y));  
}

template<class T>
Matrix<T> norm_2(const Matrix<T>& x){
  return sqrt(inner_prod(x,x));
}

template<class T>
void qr(const Matrix<T>& A, Matrix<T>& Q, Matrix<T> &R){
  // The following algorithm is taken from J. Demmel: Applied Numerical Linear Algebra (algorithm 3.1.)
  int m = A.size1();
  int n = A.size2();
  if(m<n) throw CasadiException("qr: fewer rows than columns");

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
      qi -= T(ri(0,j)) * qj;
    }

    // Normalize qi
    ri(0,i) = norm_2(trans(qi))[0];
    qi /= T(ri(0,i));

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
  if(A.size1() != b.size1()) throw CasadiException("solve: dimension mismatch");
  if(A.size1() != A.size2()) throw CasadiException("solve: A not square");
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
  if(!isScalar(a_) || !isScalar(b_)) throw CasadiException("Matrix<T>::linspace: a and b must be scalar");
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


} // namespace CasADi

#endif //SWIG




#ifdef SWIG

// map the template name to the instantiated name
#define MTT_INST(T,function_name,ns) \
%template(function_name) ns::function_name<T>;

// Define template instanciations
#define MATRIX_TOOLS_TEMPLATES(T) \
MTT_INST(T,trans,CasADi) \
MTT_INST(T,prod,CasADi) \
MTT_INST(T,append,CasADi) \
MTT_INST(T,getSub,CasADi) \
MTT_INST(T,setSub,CasADi) \
MTT_INST(T,getRow,CasADi) \
MTT_INST(T,setRow,CasADi) \
MTT_INST(T,getColumn,CasADi) \
MTT_INST(T,setColumn,CasADi) \
MTT_INST(T,isConstant,CasADi) \
MTT_INST(T,isDense,CasADi) \
MTT_INST(T,isEmpty,CasADi) \
MTT_INST(T,isInteger,CasADi) \
MTT_INST(T,isScalar,CasADi) \
MTT_INST(T,isVector,CasADi) \
MTT_INST(T,isTril,CasADi) \
MTT_INST(T,isTriu,CasADi) \
MTT_INST(T,sin,std) \
MTT_INST(T,cos,std) \
MTT_INST(T,tan,std) \
MTT_INST(T,asin,std) \
MTT_INST(T,acos,std) \
MTT_INST(T,atan,std) \
MTT_INST(T,exp,std) \
MTT_INST(T,log,std) \
MTT_INST(T,sqrt,std) \
MTT_INST(T,floor,std) \
MTT_INST(T,ceil,std) \
MTT_INST(T,fabs,std) \
MTT_INST(T,det,CasADi) \
MTT_INST(T,getMinor,CasADi) \
MTT_INST(T,cofactor,CasADi) \
MTT_INST(T,adj,CasADi) \
MTT_INST(T,inv,CasADi) \
MTT_INST(T,reshape,CasADi) \
MTT_INST(T,vec,CasADi) \
MTT_INST(T,horzcat,CasADi) \
MTT_INST(T,vertcat,CasADi) \
MTT_INST(T,inner_prod,CasADi) \
MTT_INST(T,outer_prod,CasADi) \
MTT_INST(T,norm_2,CasADi) \
MTT_INST(T,qr,CasADi) \
MTT_INST(T,solve,CasADi) \
MTT_INST(T,linspace,CasADi) \
MTT_INST(T,isZero,CasADi) \
MTT_INST(T,nnz,CasADi) \
MTT_INST(T,nnz_sym,CasADi) \
MTT_INST(T,isEqual,CasADi) \





#endif //SWIG



#endif // MATRIX_TOOLS_HPP
