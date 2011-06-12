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
/// Python's range function
std::vector<int> range(int start, int stop, int step=1, int len=std::numeric_limits<int>::max());

/// Python's range function, start = 0
std::vector<int> range(int stop);

/// Matrix product of two matrices - not available in Python since prod in numpy means elementwise multiplication
template<class T>
Matrix<T> prod(const Matrix<T> &x, const Matrix<T> &y);
#endif // SWIG

/// Product of two matrices - Python naming
template<class T>
Matrix<T> dot(const Matrix<T> &x, const Matrix<T> &y);

template<class T>
void append(Matrix<T>& expr, const Matrix<T>& add);

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
Matrix<T> reshape(const Matrix<T>& a, const std::vector<int> sz);

template<class T>
Matrix<T> reshape(const Matrix<T>& a, const CRSSparsity& sp);

template<class T>
T trace(const Matrix<T>& a);

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
bool isIdentity(const Matrix<T>& ex);

template<class T>
int nnz(const Matrix<T>& ex);

template<class T>
int nnz_sym(const Matrix<T>& ex);

template<class T>
bool isEqual(const Matrix<T>& ex1,const Matrix<T> &ex2);

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

/** \brief Repeat matrix A n times vertically and m times horizontally */
template<class T>
Matrix<T> repmat(const Matrix<T> &A, int n, int m);

/** \brief  create an n-by-n identity matrix */
template<class T>
Matrix<T> eye(int n);

/** \brief  create a matrix with all ones */
template<class T>
Matrix<T> ones(int n, int m=1);

/** \brief  create a matrix with all zeros */
template<class T>
Matrix<T> zeros(int n, int m=1);

#ifndef SWIG
/** \brief  Get the sparsity in sparse triplet format */
template<class T>
void getSparseTriplet(const Matrix<T>& A, std::vector<int>& row, std::vector<int>& col);
#endif // SWIG

/** \brief  Get the sparsity in sparse triplet format - Python style: [row,col] =  getSparseTriplet(A) */
template<class T>
std::vector<std::vector<int> > getSparseTriplet(const Matrix<T>& A);

/** \brief  Unite two matrices no overlapping sparsity */
template<class T>
Matrix<T> unite(const Matrix<T>& A, const Matrix<T>& B);

/** \brief  Make a matrix dense */
template<class T>
void makeDense(Matrix<T>& A);

} // namespace CasADi

#ifndef SWIG

// Implementations of the functions in standard namespace

#include <iterator>

namespace CasADi{
// Implementations

template<class T>
Matrix<T> trans(const Matrix<T> &x){
  return x.trans();
}

template<class T>
Matrix<T> dot(const Matrix<T> &x, const Matrix<T> &y){
  return prod<T>(x,y);
}

template<class T>
Matrix<T> prod(const Matrix<T> &x, const Matrix<T> &y){
  return x.prod(y);
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
  
  // Append the sparsity pattern
  expr.sparsityRef().append(add.sparsity());
  
  // Add the non-zeros
  expr.data().insert(expr.end(),add.begin(),add.end());
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
  if(isScalar(a)) return a.toScalar();

  // Return expression
  Matrix<T> ret = 0;

  // We expand the matrix along the first column
  for(int i=0; i<n; ++i){

    // Sum up the cofactors
    ret += a(i,0)*cofactor(a,i,0);

  }
  return ret.toScalar();
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
  casadi_assert_message(a.numel() == n*m, "reshape: number of elements must remain the same");

  Matrix<T> ret(n,m);
  for(int i=0; i<a.size1(); ++i){
    for(int el=a.rowind(i); el<a.rowind(i+1); ++el){
      int j = a.col(el);
      int k = j+i*a.size2();
      ret(k/m,k%m) = a[el];
    }
  }
  return ret;
}

template<class T>
Matrix<T> reshape(const Matrix<T>& a, const std::vector<int> sz){
  casadi_assert_message(sz.size() == 2, "reshape: must be two dimensional");
  return reshape(a,sz[0],sz[1]);
}

template<class T>
Matrix<T> reshape(const Matrix<T>& x, const CRSSparsity& sp){
  // quick return if already the right shape
  if(sp==x.sparsity())
    return x;
  
  // make sure that the number of zeros agree
  casadi_assert(x.size()==sp.size());
  
  return Matrix<T>(sp,x.data());
}

template<class T>
T trace(const Matrix<T>& a){
  casadi_assert_message(a.size1() == a.size2(), "trace: must be square");
  T res=0;
  for (int i=0; i< a.size1(); i ++) {
    res+=a(i,i);
  }
  return res;
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
  return std::sqrt(inner_prod(x,x));
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
    Matrix<T> ai = AT(i,ALL);
    Matrix<T> qi = ai;
    // The i-th column of R
    Matrix<T> ri(1,n);
  
    // subtract the projection of of qi in the previous directions from ai
    for(int j=0; j<i; ++j){
      // Get the j-th column of Q
      Matrix<T> qj = QT(j,ALL);

      ri(0,j) = prod(qj,trans(qi))(0,0); // Modified Gram-Schmidt
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
  T a = a_.toScalar();
  T b = b_.toScalar();
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
bool isIdentity(const Matrix<T>& ex) {
  
  // Make sure that the matrix is diagonal
  if(!ex.sparsity().diagonal())
    return false;
  
  // Make sure that all entries are one
  for(typename Matrix<T>::const_iterator it=ex.begin(); it!=ex.end(); ++it){
    if(!casadi_limits<T>::isOne(*it))
      return false;
  }
      
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
Matrix<T> repmat(const Matrix<T> &A, int n, int m){
  // First concatenate horizontally
  Matrix<T> row = horzcat(std::vector<Matrix<T> >(m, A));
  
  // Then vertically
  return vertcat(std::vector<Matrix<T> >(n, row));
}

template<class T>
Matrix<T> ones(int n, int m){
  return Matrix<T>(n,m,1);
}

template<class T>
Matrix<T> zeros(int n, int m){
  return Matrix<T>(n,m);
}

template<class T>
Matrix<T> eye(int n){
  return Matrix<T>(CRSSparsity::createDiagonal(n),1);
}

template<class T>
void getSparseTriplet(const Matrix<T>& A, std::vector<int>& row, std::vector<int>& col){
  col = A.sparsity().col();
  row = A.sparsity().getRow();
}

template<class T>
std::vector<std::vector<int> > getSparseTriplet(const Matrix<T>& A){
  std::vector<std::vector<int> > ret(2);
  getSparseTriplet(A,ret[0],ret[1]);
  return ret;
}

template<class T>
Matrix<T> unite(const Matrix<T>& A, const Matrix<T>& B){
  // Join the sparsity patterns
  std::vector<int> mapping;
  CRSSparsity sp = A.sparsity().patternUnion(B.sparsity(),mapping);
  
  // Create return matrix
  Matrix<T> ret(sp);
  
  // Copy sparsity
  int elA=0, elB=0;
  for(int k=0; k<mapping.size(); ++k){
    if(mapping[k]<0){
      ret[k] = A[elA++];
    } else if(mapping[k]>0){
      ret[k] = B[elB++];
    } else {
      throw CasadiException("Pattern intersection not empty");
    }
  }
  
  casadi_assert(A.size()==elA);
  casadi_assert(B.size()==elB);
  
  return ret;
}

template<class T>
void makeDense(Matrix<T>& A){
  A.makeDense(A.size1(),A.size2(),0);
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
MTT_INST(T,isIdentity) \
MTT_INST(T,nnz) \
MTT_INST(T,nnz_sym) \
MTT_INST(T,isEqual) \
MTT_INST(T,repmat) \
MTT_INST(T,getSparseTriplet) \
MTT_INST(T,unite) \
MTT_INST(T,sum) \
MTT_INST(T,sum_all) \
MTT_INST(T,trace) \
MTT_INST(T,makeDense) \

#endif //SWIG

#endif // MATRIX_TOOLS_HPP
