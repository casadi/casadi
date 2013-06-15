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
#include <algorithm>

namespace CasADi{

/// Transpose of a matrix
template<class T>
Matrix<T> trans(const Matrix<T> &x);

/// Matrix product of two matrices
template<class T>
Matrix<T> mul(const Matrix<T> &x, const Matrix<T> &y);

/// Matrix product of n matrices
template<class T>
Matrix<T> mul(const std::vector< Matrix<T> > &args);

/** \brief  check if the matrix is constant (note that false negative answers are possible)*/
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

//{@
/// Checks if vector does not contain NaN or Inf
bool isRegular(const Matrix<int>& ex);
bool isRegular(const Matrix<double>& ex);
//@}

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
  Reshapes/vectorizes the Matrix<T> such that the shape becomes (expr.numel(),1).
  Columns are stacked on top of each other.
  Same as reshape(trans(expr), expr.numel(),1)
  
    a b \n
    c d \n
    
    turns into
    
    a \n
    c \n
    b \n
    d \n
    
 */
template<class T>
Matrix<T> vec(const Matrix<T>& a);

/** \brief  make a vector
  Flattens the Matrix<T> such that the shape becomes (expr.numel(),1).
  Transposed rows are stacked on top of each other.
  Same as reshape(expr, expr.numel(),1)
  
    a b \n
    c d \n
    
    turns into
    
    a \n
    b \n
    c \n
    d \n
    
 */
template<class T>
Matrix<T> flatten(const Matrix<T>& a);

/** \brief Returns a flattened version of the Matrix, preserving only nonzeros
 */
template<class T>
Matrix<T> vecNZ(const Matrix<T>& a);

/** \brief Construct a matrix from a list of list of blocks.
*/
template<class T>
Matrix<T> blockcat(const std::vector< std::vector<Matrix<T> > > &v);

#ifndef SWIG
/** \brief Construct a matrix from 4 blocks
*/
template<class T>
Matrix<T> blockcat(const Matrix<T> &A,const Matrix<T> &B,const Matrix<T> &C,const Matrix<T> &D);
#endif // SWIG

/** \brief Concatenate a list of matrices vertically
* Alternative terminology: vertical stack, vstack, vertical append, [a;b]
*/
template<class T>
Matrix<T> vertcat(const std::vector<Matrix<T> > &v);

/** \brief Concatenate a list of matrices horizontally
* Alternative terminology: horizontal stack, hstack, horizontal append, [a b]
*/
template<class T>
Matrix<T> horzcat(const std::vector<Matrix<T> > &v);

#ifndef SWIG
template<class T>
Matrix<T> vertcat(const Matrix<T> &x, const Matrix<T> &y);

template<class T>
Matrix<T> horzcat(const Matrix<T> &x, const Matrix<T> &y);
#endif // SWIG

template<class T>
/** \brief  concatenate vertically while vectorizing all arguments with vec */
Matrix<T> veccat(const std::vector< Matrix<T> >& comp);

template<class T>
/** \brief  concatenate vertically while vectorizing all arguments with vecNZ */
Matrix<T> vecNZcat(const std::vector< Matrix<T> >& comp);

/** \brief Inner product of two matrices
        Equals
        \code
        sumAll(x*y)
        \endcode
        with x and y matrices of the same dimension
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
 * More stable than the classical Gram-Schmidt, but may break down if the columns of A are nearly linearly dependent
 * See J. Demmel: Applied Numerical Linear Algebra (algorithm 3.1.). 
 * Note that in SWIG, Q and R are returned by value. */
template<class T>
#ifndef SWIG
void qr(const Matrix<T>& A, Matrix<T>& Q, Matrix<T>& R);
#else // SWIG
void qr(const Matrix<T>& A, Matrix<T>& OUTPUT, Matrix<T>& OUTPUT);
#endif

/** \brief  Solve a system of equations: A*x = b 
The solve routine works similar to Matlab's backslash when A is square and nonsingular. The algorithm
used is the following:
  1. A simple forward or backward substitution if A is upper or lower triangular
  2. If the linear system is at most 3-by-3, form the inverse via minor expansion and multiply
  3. Permute the variables and equations as to get a (structurally) nonzero diagonal, then perform
     a QR factorization without pivoting and solve the factorized system.

Note 1: If there are entries of the linear system known to be zero, these will be removed. Elements
  that are very small, or will evaluate to be zero, can still cause numerical errors, due to the lack
  of pivoting (which is not possible since cannot compare the size of entries)

Note 2: When permuting the linear system, a BLT (block lower triangular) transformation is formed.
  Only the permutation part of this is however used. An improvement would be to solve block-by-block
  if there are multiple BLT blocks.
  
*/
template<class T>
Matrix<T> solve(const Matrix<T>& A, const Matrix<T>& b);

/** \brief  check if the matrix is 0 (note that false negative answers are possible)*/
template<class T>
bool isZero(const Matrix<T>& ex);

/** \brief  check if the matrix is 1 (note that false negative answers are possible)*/
template<class T>
bool isOne(const Matrix<T>& ex);

/** \brief  check if the matrix is -1 (note that false negative answers are possible)*/
template<class T>
bool isMinusOne(const Matrix<T>& ex);

/** \brief  check if the matrix is an identity matrix (note that false negative answers are possible)*/
template<class T>
bool isIdentity(const Matrix<T>& ex);

template<class T>
int nnz(const Matrix<T>& ex);

template<class T>
int nnz_sym(const Matrix<T>& ex);

/** \brief Check if two expressions are equal
*
*  Might very well give false negatives
*
*   Note: does not work when CasadiOptions.setSimplificationOnTheFly(False) was called
*/
template<class T>
bool isEqual(const Matrix<T>& ex1,const Matrix<T> &ex2);

/// Make the vector 1-norm of an Matrix<T>
template<class T>
Matrix<T> norm_1(const Matrix<T>& x);
  
/// Make the vector 2-norm (Frobenius Norm) of an Matrix<T>
template<class T>
Matrix<T> norm_2(const Matrix<T>& x);

/// Make the vector 2-norm (Frobenius Norm) squared of an Matrix<T>
template<class T>
Matrix<T> norm_22(const Matrix<T>& x);

/// Return summation of all elements
template<class T>
T sumAll(const Matrix<T> &x); 

/** \brief Return a row-wise summation of elements */
template<class T>
Matrix<T> sumRows(const Matrix<T> &x);

/** \brief Return a column-wise summation of elements */
template<class T>
Matrix<T> sumCols(const Matrix<T> &x);

#ifdef SWIG
/// Returns true only if every element in the matrix is true
template<class T>
T all(const Matrix<T> &x); 

/// Returns true if any element in the matrix is true
template<class T>
T any(const Matrix<T> &x); 
#endif //SWIG

/** \brief Repeat matrix A n times vertically and m times horizontally */
template<class T>
Matrix<T> repmat(const Matrix<T> &A, int n, int m);

/** \brief  Evaluate a polynomial with coefficeints p in x */
template<class T>
Matrix<T> polyval(const Matrix<T>& p, const Matrix<T>& x);

/** \brief   Get the diagonal of a matrix or construct a diagonal
When the input is square, the diagonal elements are returned.
If the input is vector-like, a diagonal matrix is constructed with it. */
template<class T>
Matrix<T> diag(const Matrix<T> &A);

/** \brief   Construct a matrix with given block on the diagonal */
template<class T>
Matrix<T> blkdiag(const std::vector< Matrix<T> > &A);

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

/** \brief  Make a matrix dense */
template<class T>
Matrix<T> densify(const Matrix<T>& A);

#ifndef SWIGOCTAVE
/** \brief  Make a matrix dense */
template<class T>
Matrix<T> full(const Matrix<T>& A);
#endif // SWIGOCTAVE

/** \brief  Make a matrix sparse by removing numerical */
template<class T>
void makeSparse(Matrix<T>& A);

#ifndef SWIGOCTAVE
/** \brief  Make a matrix sparse by removing numerical */
template<class T>
Matrix<T> sparse(const Matrix<T>& A);
#endif // SWIGOCTAVE

/** \brief  Check if the matrix has any zero entries which are not structural zeros */
template<class T>
bool hasNonStructuralZeros(const Matrix<T>& A);

/// same as: res += mul(A,v)
template<typename T>
void addMultiple(const Matrix<T>& A, const std::vector<T>& v, std::vector<T>& res, bool trans_A=false);

/// Get a pointer to the data contained in the vector
template<typename T>
T* getPtr(Matrix<T> &v);
  
/// Get a pointer to the data contained in the vector
template<typename T>
const T* getPtr(const Matrix<T> &v);

/** \brief Create a new matrix with a given sparsity pattern but with the nonzeros taken from an existing matrix */
template<typename T>
Matrix<T> project(const Matrix<T>& A, const CRSSparsity& sparsity);

/// Obtain the structural rank of a sparsity-pattern
template<typename T>
int sprank(const Matrix<T>& A);

} // namespace CasADi

// Global namespace

#ifndef SWIG
#include <iterator>

namespace CasADi{
// Implementations

template<class T>
Matrix<T> trans(const Matrix<T> &x){
  return x.trans();
}

template<class T>
Matrix<T> mul(const Matrix<T> &x, const Matrix<T> &y){
  return x.mul(y);
}

template<class T>
Matrix<T> mul(const std::vector< Matrix<T> > &args){
  casadi_assert_message(args.size()>=1,"mul(std::vector< Matrix<T> > &args): supplied list must not be empty.");
  if (args.size()==1) return args[0];
  Matrix<T> ret = args[0].mul(args[1]);
  for (int i=2;i<args.size();++i) {
    ret = ret.mul(args[i]);
  }
  return ret;
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
    if(!casadi_limits<T>::isInteger(ex.data()[k])) // if an element is not integer
      return false;
    
  // Integer if reached this point
  return true;
}

template<class T>
bool isConstant(const Matrix<T>& ex){
  // loop over non-zero elements
  for(int k=0; k<ex.size(); ++k) 
    if(!casadi_limits<T>::isConstant(ex.data()[k])) // if an element is not constant
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

  // Trivial case 2 x 2
  if(n==2) return a.elem(0,0) * a.elem(1,1) - a.elem(1,0) * a.elem(0,1);
  
  // Return expression
  Matrix<T> ret = 0;
  
  // Find out which is the best direction to expand along

  // Build up an IMatrix with ones on the non-zeros
  Matrix<int> sp = IMatrix(a.sparsity(),1);
  
  // Have a count of the nonzeros for each column
  Matrix<int> col_count = sumRows(sp);
  
  // A blank column? determinant is structurally zero
  if (!col_count.dense()) return 0;

  // Have a count of the nonzeros for each row
  Matrix<int> row_count = trans(sumCols(sp));
  
  // A blank row? determinant is structurally zero
  if (!col_count.dense()) return 0;
  
  int min_col = std::distance(col_count.data().begin(), std::min_element(col_count.data().begin(),col_count.data().end()));
  int min_row = std::distance(row_count.data().begin(), std::min_element(row_count.data().begin(),row_count.data().end()));
  
  if (min_col <= min_row) {
    // Expand along column j
    int j = col_count.sparsity().col(min_col);
    
    Matrix<T> col = a(range(n),j);

    std::vector< int > row_i = col.sparsity().getRow();

    for(int k=0; k<col.size(); ++k) {
      // Sum up the cofactors
      ret += col.at(k)*cofactor(a,row_i.at(k),j);
    }
    return ret.toScalar();
  } else {
    // Expand along row i
    int i = row_count.sparsity().col(min_row);

    Matrix<T> row = a(i,range(n));
    
    const std::vector< int > &col_i = row.sparsity().col();

    for(int k=0; k<row.size(); ++k) {
      // Sum up the cofactors
      ret += row.at(k)*cofactor(a,i,col_i.at(k));
    }
    return ret.toScalar();
  }
 
}

template<class T>
T getMinor(const Matrix<T> &x, int i, int j){
  int n = x.size1();
  casadi_assert_message(n == x.size2(), "getMinor: matrix must be square");

  // Trivial return if scalar
  if(n==1) return 1;

  // Remove row i and column j
  Matrix<T> M(n-1,n-1);
  
  std::vector<int> row = x.sparsity().getRow();
  const std::vector<int> &col = x.sparsity().col();

  for(int k=0;k<x.size();++k) {
    int i1 = row[k];
    int j1 = col[k];

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

  // Temporary placeholder
  T temp;
  
  // Cofactor matrix
  Matrix<T> C(n,n);
  for(int i=0; i<n; ++i)
    for(int j=0; j<n; ++j) {
      temp = cofactor(a,i,j);
      if (!casadi_limits<T>::isZero(temp))
        C(i,j) = temp;
    }
  
  return trans(C);
}

template<class T>
Matrix<T> inv(const Matrix<T>& a){
  // laplace formula
  return adj(a)/det(a);
}

template<class T>
Matrix<T> reshape(const Matrix<T>& a, int n, int m){
  CRSSparsity sp = a.sparsity().reshape(n,m);
  return Matrix<T>(sp,a.data());
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
    res+=a.elem(i,i);
  }
  return res;
}

template<class T>
Matrix<T> vec(const Matrix<T>& a){
  Matrix<T> ret = reshape(trans(a),a.numel(),1);
  return ret;
}

template<class T>
Matrix<T> flatten(const Matrix<T>& a){
  Matrix<T> ret = reshape(a,a.numel(),1);
  return ret;
}


template<class T>
Matrix<T> vecNZ(const Matrix<T>& a){
  return Matrix<T>(vec(a).data());
}

template<class T>
Matrix<T> blockcat(const std::vector< std::vector<Matrix<T> > > &v) {
  std::vector< Matrix<T> > ret;
  for(int i=0; i<v.size(); ++i)
    ret.push_back(horzcat(v[i]));
  return vertcat(ret);
}

template<class T>
Matrix<T> blockcat(const Matrix<T> &A,const Matrix<T> &B,const Matrix<T> &C,const Matrix<T> &D) {
  return vertcat(horzcat(A,B),horzcat(C,D));
}

template<class T>
Matrix<T> vertcat(const std::vector<Matrix<T> > &v){
  Matrix<T> ret;
  for(int i=0; i<v.size(); ++i)
    ret.append(v[i]);
  return ret;
}

template<class T>
Matrix<T> horzcat(const std::vector<Matrix<T> > &v){
  Matrix<T> ret;
  for(int i=0; i<v.size(); ++i)
    ret.append(trans(v[i]));
  return trans(ret);  
}

template<class T>
Matrix<T> vertcat(const Matrix<T> &x, const Matrix<T> &y){
  Matrix<T> xy = x;
  xy.append(y);
  return xy;
}

template<class T>
Matrix<T> horzcat(const Matrix<T> &x, const Matrix<T> &y){
  return trans(vertcat(trans(x),trans(y)));
}

template<class T>
Matrix<T> veccat(const std::vector< Matrix<T> >& comp) {
  Matrix<T> (&f)(const Matrix<T>&) = vec;
  return vertcat(applymap(f,comp));
}

template<class T>
Matrix<T> vecNZcat(const std::vector< Matrix<T> >& comp) {
  Matrix<T> (&f)(const Matrix<T>&) = vecNZ;
  return vertcat(applymap(f,comp));
}

template<class T>
Matrix<T> inner_prod(const Matrix<T> &x, const Matrix<T> &y){
  casadi_assert_message(x.shape()==y.shape(), "inner_prod: Dimension mismatch");
  return Matrix<T>(sumAll(x*y));
}

template<class T>
Matrix<T> outer_prod(const Matrix<T> &x, const Matrix<T> &y){
  casadi_assert_message(x.vector() && y.vector(), "outer_prod: arguments must be vectors");
  return mul(x,trans(y));  
}

template<class T>
T sumAll(const Matrix<T> &x) {
  // Sum non-zero elements
  T res=0;
  for(int k=0; k<x.size(); k++){
    res += x.data()[k];
  }
  return res;
}

template<class T>
Matrix<T> sumRows(const Matrix<T> &x) {
  return mul(Matrix<T>::ones(1,x.size1()),x);
}

template<class T>
Matrix<T> sumCols(const Matrix<T> &x) {
  return mul(x,Matrix<T>::ones(x.size2(),1));
}

template<class T>
T all(const Matrix<T> &x) {
  if (!x.dense()) return false;
  T ret=1;
  for (int i=0;i<x.size();++i) {
  ret = ret && x.at(i)==1;
  }
  return ret;
}

template<class T>
T any(const Matrix<T> &x) {
  if (!x.dense()) return false;
  T ret=0;
  for (int i=0;i<x.size();++i) {
  ret = ret || x.at(i)==1;
  }
  return ret;
}


template<class T>
Matrix<T> norm_1(const Matrix<T>& x){
  return sumAll(fabs(x));
}

template<class T>
Matrix<T> norm_2(const Matrix<T>& x){
  return sqrt(1.0*sumAll(x*x));
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
  
    // subtract the projection of qi in the previous directions from ai
    for(int j=0; j<i; ++j){
      
      // Get the j-th column of Q
      Matrix<T> qj = QT(j,ALL);

      ri(0,j) = mul(qj,trans(qi)); // Modified Gram-Schmidt
      // ri[j] = inner_prod(qj,ai); // Classical Gram-Schmidt
     
      // Remove projection in direction j
      qi -= ri(0,j) * qj;
    }

    // Normalize qi
    ri(0,i) = norm_2(trans(qi));
    qi /= ri(0,i);

    // Update RT and QT
    QT.append(qi);
    RT.append(ri);
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
    // forward substitution if lower triangular
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
    // backward substitution if upper triangular
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
  } else if(hasNonStructuralZeros(A)){
    
    // If there are structurally nonzero entries that are known to be zero, remove these and rerun the algorithm
    Matrix<T> A_sparse = A;
    makeSparse(A_sparse);
    return solve(A_sparse,b);

  } else {
    
    // Make a BLT transformation of A
    std::vector<int> rowperm, colperm, rowblock, colblock, coarse_rowblock, coarse_colblock;
    A.sparsity().dulmageMendelsohn(rowperm, colperm, rowblock, colblock, coarse_rowblock, coarse_colblock);

    // Get the inverted column permutation
    std::vector<int> inv_colperm(colperm.size());
    for(int k=0; k<colperm.size(); ++k)
      inv_colperm[colperm[k]] = k;
    
    // Permute the right hand side
    Matrix<T> bperm(0,b.size2());
    for(int i=0; i<b.size1(); ++i){
      bperm.resize(i+1,b.size2());
      for(int el=b.rowind(rowperm[i]); el<b.rowind(rowperm[i]+1); ++el){
        bperm(i,b.col(el)) = b[el];
      }
    }

    // Permute the linear system
    Matrix<T> Aperm(0,A.size2());
    for(int i=0; i<A.size1(); ++i){
      Aperm.resize(i+1,A.size2());
      for(int el=A.rowind(rowperm[i]); el<A.rowind(rowperm[i]+1); ++el){
        Aperm(i,inv_colperm[A.col(el)]) = A[el];
      }
    }
    
    // Permuted solution
    Matrix<T> xperm;
    
    // Solve permuted system
    if(isTril(Aperm)){
      
      // Forward substitution if lower triangular after sorting the equations
      xperm = solve(Aperm,bperm);
      
    } else if(A.size1()<=3){
      
      // Form inverse by minor expansion and multiply if very small (up to 3-by-3)
      xperm = mul(inv(A),bperm);

    } else {
      
      // Make a QR factorization
      Matrix<T> Q,R;
      qr(Aperm,Q,R);

      // Solve the factorized system (note that solve will now be fast since it is triangular)
      xperm = solve(R,mul(trans(Q),bperm));
    }
    
    // Permute back the solution
    Matrix<T> x(0,xperm.size2());
    for(int i=0; i<xperm.size1(); ++i){
      x.resize(i+1,xperm.size2());
      for(int el=xperm.rowind(inv_colperm[i]); el<xperm.rowind(inv_colperm[i]+1); ++el){
        x(i,xperm.col(el)) = xperm[el];
      }
    }
    return x;
  }
}

template<class T>
bool isOne(const Matrix<T>& ex){  
  if(!ex.dense()){
    return false;
  }
  
  // loop over non-zero elements
  for(int el=0; el<ex.size(); ++el)
    if(!casadi_limits<T>::isOne(ex.at(el)))
      return false;
  
  return true;
}

template<class T>
bool isMinusOne(const Matrix<T>& ex){  
  if(!ex.dense()){
    return false;
  }
  
  // loop over non-zero elements
  for(int el=0; el<ex.size(); ++el)
    if(!casadi_limits<T>::isMinusOne(ex.at(el)))
      return false;
  
  return true;
}

template<class T>
bool isZero(const Matrix<T>& ex) {  

  // loop over (potentially) non-zero elements
  for(int el=0; el<ex.size(); ++el)
    if(!casadi_limits<T>::isZero(ex.at(el)))
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
Matrix<T> diag(const Matrix<T>&A){
  // Nonzero mapping
  std::vector<int> mapping;
  // Get the sparsity
  CRSSparsity sp = A.sparsity().diag(mapping);
  
  Matrix<T> ret = Matrix<T>(sp);
  
  for (int k=0;k<mapping.size();k++) ret[k] = A[mapping[k]];
  return ret;
}

/** \brief   Construct a matrix with given block on the diagonal */
template<class T>
Matrix<T> blkdiag(const std::vector< Matrix<T> > &A) {
  int n = 0;
  int m = 0;
  
  std::vector<int> rowind(1,0);
  std::vector<int> col;
  std::vector<T> data;
  
  int nz = 0;
  for (int i=0;i<A.size();++i) {
    data.insert(data.end(),A[i].data().begin(),A[i].data().end());
    const std::vector<int> &rowind_ = A[i].rowind();
    const std::vector<int> &col_ = A[i].col();
    for (int k=1;k<rowind_.size();++k) {
      rowind.push_back(rowind_[k]+nz);
    }
    for (int k=0;k<col_.size();++k) {
      col.push_back(col_[k]+m);
    }
    n+= A[i].size1();
    m+= A[i].size2();
    nz+= A[i].size();
  }
  
  return Matrix<T>(n,m,col,rowind,data);
  
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
  std::vector<unsigned char> mapping;
  CRSSparsity sp = A.sparsity().patternUnion(B.sparsity(),mapping);
  
  // Create return matrix
  Matrix<T> ret(sp);
  
  // Copy sparsity
  int elA=0, elB=0;
  for(int k=0; k<mapping.size(); ++k){
    if(mapping[k]==1){
      ret.data()[k] = A.data()[elA++];
    } else if(mapping[k]==2){
      ret.data()[k] = B.data()[elB++];
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

template<class T>
Matrix<T> densify(const Matrix<T>& A){
  Matrix<T> ret(A);
  makeDense(ret);
  return ret;
}

template<class T>
Matrix<T> full(const Matrix<T>& A){
  return densify(A);
}

template<class T>
void makeSparse(Matrix<T>& A){
  // Quick return if there are no structurally zero entries
  if(!hasNonStructuralZeros(A))
    return;
  
  // Start with a matrix with no rows
  Matrix<T> Asp(0,A.size2());

  // Loop over the rows
  for(int i=0; i<A.size1(); ++i){
    // Resize the matrix to accomodate the row
    Asp.resize(i+1,A.size2());
    
    // Loop over the existing, possible nonzeros
    for(int el=A.rowind(i); el<A.rowind(i+1); ++el){
      
      // If it is not known to be a zero
      if(!casadi_limits<T>::isZero(A.at(el))){
        
        // Get the column
        int j=A.col(el);
      
        // Save to new, sparse matrix
        Asp(i,j) = A.at(el);
      }
    }
  }
  
  // Save to A
  A = Asp;
}

template<class T>
Matrix<T> sparse(const Matrix<T>& A){
  Matrix<T> ret(A);
  makeSparse(ret);
  return ret;
}

template<class T>
bool hasNonStructuralZeros(const Matrix<T>& A){
  // Loop over the structural nonzeros
  for(int el=0; el<A.size(); ++el){
    
    // Check if the structural nonzero is known to be zero
    if(casadi_limits<T>::isZero(A.at(el)))
      return true;
  }
  
  // No known zeros amongst the structurally nonzero entries
  return false;
}
  
template<class T>
Matrix<T> polyval(const Matrix<T>& p, const Matrix<T>& x){
  casadi_assert_message(isDense(p),"polynomial coefficients vector must be a vector");
  casadi_assert_message(isVector(p) && p.size()>0,"polynomial coefficients must be a vector");
  Matrix<T> ret = p[0];
  for(int i=1; i<p.size(); ++i){
    ret = ret*x + p[i];
  }
  return ret;
}

template<typename T>
void addMultiple(const Matrix<T>& A, const std::vector<T>& v, std::vector<T>& res, bool trans_A){
  // Get dimension and sparsity
  int d1=A.size1(), d2=A.size2();
  const std::vector<int> &rowind=A.rowind();
  const std::vector<int> &col=A.col();
  const std::vector<T>& data = A.data();

  // Assert consistent dimensions
  if(trans_A){
    casadi_assert(v.size()==d1);
    casadi_assert(res.size()==d2);
  } else {
    casadi_assert(v.size()==d2);
    casadi_assert(res.size()==d1);
  }

  // Carry out multiplication
  for(int i=0; i<d1; ++i){ // loop over rows
    for(int el=rowind[i]; el<rowind[i+1]; ++el){ // loop over the non-zero elements
      int j=col[el];  // column
      // Add scalar product
      if(trans_A){
        res[j] += v[i]*data[el];
      } else {
        res[i] += v[j]*data[el];
      }
    }
  }
}

  template<typename T>
  T* getPtr(Matrix<T> &v){
    if(v.empty())
      return 0;
    else
      return &v.front();
  }
  
  template<typename T>
  const T* getPtr(const Matrix<T> &v){
    if(v.empty())
      return 0;
    else
      return &v.front();
  }

  template<typename T>
  Matrix<T> project(const Matrix<T>& A, const CRSSparsity& sparsity){
    // Check dimensions
    if(!(A.empty() && sparsity.numel()==0)){
      casadi_assert_message(A.size1()==sparsity.size1() && A.size2()==sparsity.size2(),
                            "Shape mismatch. Expecting " << A.dimString() << ", but got " << 
                            sparsity.dimString() << " instead.");
    }
    
    // Return value
    Matrix<T> ret(sparsity,0);
    
    // Get the elements of the known matrix
    std::vector<int> known_ind = A.sparsity().getElements(false);
      
    // Find the corresponding nonzeros in the return matrix
    sparsity.getNZInplace(known_ind);
      
    // Set the element values
    const std::vector<T>& A_data = A.data();
    std::vector<T>& ret_data = ret.data();
    for(int k=0; k<known_ind.size(); ++k){
      if(known_ind[k]!=-1){
        ret_data[known_ind[k]] = A_data[k];
      }
    }
    return ret;
  }

  template<typename T>
  int sprank(const Matrix<T>& A) {
    return rank(A.sparsity());
  }
  
} // namespace CasADi


#endif //SWIG




#ifdef SWIG

// map the template name to the instantiated name
#define MTT_INST(T,function_name) \
%template(function_name) CasADi::function_name < T >;

// Define template instanciations
#define MATRIX_TOOLS_TEMPLATES_COMMON(T) \
MTT_INST(T,trans) \
MTT_INST(T,mul) \
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
MTT_INST(T,flatten) \
MTT_INST(T,vecNZ) \
MTT_INST(T,blockcat) \
MTT_INST(T,horzcat) \
MTT_INST(T,vertcat) \
MTT_INST(T,inner_prod) \
MTT_INST(T,outer_prod) \
MTT_INST(T,norm_2) \
MTT_INST(T,qr) \
MTT_INST(T,solve) \
MTT_INST(T,isZero) \
MTT_INST(T,isOne) \
MTT_INST(T,isMinusOne) \
MTT_INST(T,isIdentity) \
MTT_INST(T,nnz) \
MTT_INST(T,nnz_sym) \
MTT_INST(T,isEqual) \
MTT_INST(T,repmat) \
MTT_INST(T,getSparseTriplet) \
MTT_INST(T,unite) \
MTT_INST(T,sumCols) \
MTT_INST(T,sumRows) \
MTT_INST(T,sumAll) \
MTT_INST(T,trace) \
MTT_INST(T,makeDense) \
MTT_INST(T,densify) \
MTT_INST(T,makeSparse) \
MTT_INST(T,hasNonStructuralZeros) \
MTT_INST(T,diag) \
MTT_INST(T,blkdiag) \
MTT_INST(T,polyval) \
MTT_INST(T,addMultiple) \
MTT_INST(T,veccat) \
MTT_INST(T,vecNZcat) \
MTT_INST(T,project) \
MTT_INST(T,sprank) 
#endif //SWIG

#ifdef SWIGOCTAVE
#define MATRIX_TOOLS_TEMPLATES(T) MATRIX_TOOLS_TEMPLATES_COMMON(T)
#else
#define MATRIX_TOOLS_TEMPLATES(T) \
MATRIX_TOOLS_TEMPLATES_COMMON(T) \
MTT_INST(T,sparse) \
MTT_INST(T,full)
#endif //SWIGOCTAVE

#endif // MATRIX_TOOLS_HPP
