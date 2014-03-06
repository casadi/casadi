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

#include "../options_functionality.hpp"

#include "sparsity_tools.hpp"

namespace CasADi{

  /// Transpose of a matrix
  template<class T>
  Matrix<T> trans(const Matrix<T> &x);

  /** \brief  Matrix product of two matrices
   *
   * With optional sp_z you can specify the sparsity of the result
   * A typical use case might be where the product is only constructed to 
   * inspect the trace of it. sp_z diagonal will be more efficient then. 
   */
  template<class T>
  Matrix<T> mul(const Matrix<T> &x, const Matrix<T> &y, const Sparsity& sp_z=Sparsity());

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

  /** \brief  Check if a matrix is lower triangular (complexity ~ A.size2()) */
  template<class T>
  bool isTril(const Matrix<T> &A);

  /** \brief  Check if a matrix is upper triangular (complexity ~ A.size2()) */
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
  Matrix<T> reshape(const Matrix<T>& a, int nrow, int ncol);

  template<class T>
  Matrix<T> reshape(const Matrix<T>& a, std::pair<int,int> rc);

  template<class T>
  Matrix<T> reshape(const Matrix<T>& a, const Sparsity& sp);

  template<class T>
  T trace(const Matrix<T>& a);

  /** \brief  make a vector
      Reshapes/vectorizes the Matrix<T> such that the shape becomes (expr.numel(),1).
      Columns are stacked on top of each other.
      Same as reshape(expr, expr.numel(),1)
  
      a c \n
      b d \n
    
      turns into
    
      a \n
      b \n
      c \n
      d \n
    
  */
  template<class T>
  Matrix<T> vec(const Matrix<T>& a);

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
   *
   *   horzcat(horzsplit(x,...)) = x
   */
  template<class T>
  Matrix<T> horzcat(const std::vector<Matrix<T> > &v);

  /** \brief  split vertically, retaining groups of cols
   * \param offset List of all start cols for each group
   *      the last col group will run to the end.
   *
   *   horzcat(horzsplit(x,...)) = x
   */
  template<class T>
  std::vector<Matrix<T> > horzsplit(const Matrix<T> &v, const std::vector<int>& offset);

  /** \brief  split vertically, retaining fixed-sized groups of cols
   * \param incr Size of each group of cols
   *
   *   horzcat(horzsplit(x,...)) = x
   */
  template<class T>
  std::vector<Matrix<T> > horzsplit(const Matrix<T> &v, int incr=1);

  /** \brief Concatenate a list of matrices horizontally
   * Alternative terminology: horizontal stack, hstack, horizontal append, [a b]
   *
   *   vertcat(vertsplit(x,...)) = x
   */
  template<class T>
  Matrix<T> vertcat(const std::vector<Matrix<T> > &v);

  /** \brief  split horizontally, retaining groups of rows
   * \param output_offset List of all start rows for each group
   *      the last row group will run to the end.
   *
   *   vertcat(vertsplit(x,...)) = x
   */
  template<class T>
  std::vector<Matrix<T> > vertsplit(const Matrix<T> &v, const std::vector<int>& offset);

  /** \brief  split horizontally, retaining fixed-sized groups of rows
   * \param incr Size of each group of rows
   *
   *   vertcat(vertsplit(x,...)) = x
   */
  template<class T>
  std::vector<Matrix<T> > vertsplit(const Matrix<T> &v, int incr=1);


  /** \brief  chop up into blocks
   * \brief vert_offset Defines the boundaries of the block rows
   * \brief horz_offset Defines the boundaries of the block columns
   *
   *   blockcat(blocksplit(x,...,...)) = x
   */
  template<class T>
  std::vector< std::vector< Matrix<T> > > blocksplit(const Matrix<T>& x, const std::vector<int>& vert_offset, const std::vector<int>& horz_offset);

  /** \brief  chop up into blocks
   * \brief vert_incr Defines the increment for block boundaries in row dimension
   * \brief horz_incr Defines the increment for block boundaries in column dimension
   *
   *   blockcat(blocksplit(x,...,...)) = x
   */
  template<class T>
  std::vector< std::vector< Matrix<T> > > blocksplit(const Matrix<T>& x, int vert_incr = 1, int horz_incr = 1);

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
   * More stable than the classical Gram-Schmidt, but may break down if the rows of A are nearly linearly dependent
   * See J. Demmel: Applied Numerical Linear Algebra (algorithm 3.1.). 
   * Note that in SWIG, Q and R are returned by value. */
  template<class T>
#ifndef SWIG
  void qr(const Matrix<T>& A, Matrix<T>& Q, Matrix<T>& R);
#else // SWIG
  void qr(const Matrix<T>& A, Matrix<T>& OUTPUT, Matrix<T>& OUTPUT);
#endif

  /** \brief Computes the nullspace of a matrix A
  *
  * Finds Z m-by-(m-n) such that AZ = 0 
  * with A n-by-m with m > n
  *
  * Assumes A is full rank
  *
  * Inspired by Numerical Methods in Scientific Computing by Ake Bjorck
  */
  template<class T>
  Matrix<T> nullspace(const Matrix<T>& A);

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
  
  /** \brief Computes the Moore-Penrose pseudo-inverse
  * 
  * If the matrix A is fat (size2>size1), mul(A,pinv(A)) is unity.
  * If the matrix A is slender (size1<size2), mul(pinv(A),A) is unity.
  *
  */
  template<class T>
  Matrix<T> pinv(const Matrix<T>& A);
  
  /** \brief Solve a system of equations: A*x = b 
  */
  Matrix<double> solve(const Matrix<double>& A, const Matrix<double>& b, linearSolverCreator lsolver, const Dictionary& dict = Dictionary());
  
  
  /** \brief Computes the Moore-Penrose pseudo-inverse
  * 
  * If the matrix A is fat (size1>size2), mul(A,pinv(A)) is unity.
  * If the matrix A is slender (size2<size1), mul(pinv(A),A) is unity.
  *
  */
  Matrix<double> pinv(const Matrix<double>& A,linearSolverCreator lsolver, const Dictionary& dict = Dictionary());
  
  /** \brief Kronecker tensor product
  *
  * Creates a block matrix in which each element (i,j) is a_ij*b 
  */
  template<class T>
  Matrix<T> kron(const Matrix<T>& a, const Matrix<T>& b);

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

  /** \brief  Frobenius norm  */
  template<class T>
  Matrix<T> norm_F(const Matrix<T> &x);

  /** \brief  2-norm  */
  template<class T>
  Matrix<T> norm_2(const Matrix<T> &x);

  /** \brief 1-norm  */
  template<class T>
  Matrix<T> norm_1(const Matrix<T> &x);

  /** \brief Infinity-norm */
  template<class T>
  Matrix<T> norm_inf(const Matrix<T> &x);

  /// Return summation of all elements
  template<class T>
  Matrix<T> sumAll(const Matrix<T> &x); 

  /** \brief Return a col-wise summation of elements */
  template<class T>
  Matrix<T> sumCols(const Matrix<T> &x);

  /** \brief Return a row-wise summation of elements */
  template<class T>
  Matrix<T> sumRows(const Matrix<T> &x);

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
  void getSparseTriplet(const Matrix<T>& A, std::vector<int>& col, std::vector<int>& row);
#endif // SWIG

  /** \brief  Get the sparsity in sparse triplet format - Python style: [col,row] =  getSparseTriplet(A) */
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

  /** \brief  Make a matrix sparse by removing numerical zeros */
  template<class T>
  void makeSparse(Matrix<T>& A, double tol=0);

#ifndef SWIGOCTAVE
  /** \brief  Make a matrix sparse by removing numerical zeros*/
  template<class T>
  Matrix<T> sparse(const Matrix<T>& A, double tol=0);
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
  Matrix<T> project(const Matrix<T>& A, const Sparsity& sparsity);

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
  Matrix<T> mul(const Matrix<T> &x, const Matrix<T> &y, const Sparsity &sp_z){
    return x.mul(y,sp_z);
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
    return ex.size2()==1 && ex.size1()==1;
  }

  template<class T>
  bool isVector(const Matrix<T>& ex){
    return ex.size1()==1;
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
    return A.sparsity().tril();
  }

  template<class T>
  bool isTriu(const Matrix<T> &A){
    return A.sparsity().triu();
  }

  template<class T>
  T det(const Matrix<T>& a){
    int n = a.size2();
    casadi_assert_message(n == a.size1(),"matrix must be square");

    // Trivial return if scalar
    if(isScalar(a)) return a.toScalar();

    // Trivial case 2 x 2
    if(n==2) return a.elem(0,0) * a.elem(1,1) - a.elem(0,1) * a.elem(1,0);
  
    // Return expression
    Matrix<T> ret = 0;
  
    // Find out which is the best direction to expand along

    // Build up an IMatrix with ones on the non-zeros
    Matrix<int> sp = IMatrix(a.sparsity(),1);
  
    // Have a count of the nonzeros for each row
    Matrix<int> row_count = sumCols(sp);
  
    // A blank row? determinant is structurally zero
    if (!row_count.dense()) return 0;

    // Have a count of the nonzeros for each col
    Matrix<int> col_count = trans(sumRows(sp));
  
    // A blank col? determinant is structurally zero
    if (!row_count.dense()) return 0;
  
    int min_row = std::distance(row_count.data().begin(), std::min_element(row_count.data().begin(),row_count.data().end()));
    int min_col = std::distance(col_count.data().begin(), std::min_element(col_count.data().begin(),col_count.data().end()));
  
    if (min_row <= min_col) {
      // Expand along row j
      int j = row_count.sparsity().row(min_row);
    
      Matrix<T> row = a(j,range(n));

      std::vector< int > col_i = row.sparsity().getCol();

      for(int k=0; k<row.size(); ++k) {
        // Sum up the cofactors
        ret += row.at(k)*cofactor(a,col_i.at(k),j);
      }
      return ret.toScalar();
    } else {
      // Expand along col i
      int i = col_count.sparsity().row(min_col);

      Matrix<T> col = a(range(n),i);
    
      const std::vector< int > &row_i = col.sparsity().row();

      for(int k=0; k<col.size(); ++k) {
        // Sum up the cofactors
        ret += col.at(k)*cofactor(a,i,row_i.at(k));
      }
      return ret.toScalar();
    }
 
  }

  template<class T>
  T getMinor(const Matrix<T> &x, int i, int j){
    int n = x.size2();
    casadi_assert_message(n == x.size1(), "getMinor: matrix must be square");

    // Trivial return if scalar
    if(n==1) return 1;

    // Remove col i and row j
    Matrix<T> M = Matrix<T>::sparse(n-1,n-1);
  
    std::vector<int> col = x.sparsity().getCol();
    const std::vector<int> &row = x.sparsity().row();

    for(int k=0;k<x.size();++k) {
      int i1 = col[k];
      int j1 = row[k];

      if(i1 == i || j1 == j)
        continue;

      int i2 = (i1<i)?i1:i1-1;
      int j2 = (j1<j)?j1:j1-1;

      M(j2,i2) = x(j1,i1);
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
    int n = a.size2();
    casadi_assert_message(n == a.size1(),"adj: matrix must be square");

    // Temporary placeholder
    T temp;
  
    // Cofactor matrix
    Matrix<T> C = Matrix<T>::sparse(n,n);
    for(int i=0; i<n; ++i)
      for(int j=0; j<n; ++j) {
        temp = cofactor(a,i,j);
        if (!casadi_limits<T>::isZero(temp))
          C(j,i) = temp;
      }
  
    return trans(C);
  }

  template<class T>
  Matrix<T> inv(const Matrix<T>& a){
    // laplace formula
    return adj(a)/det(a);
  }

  template<class T>
  Matrix<T> reshape(const Matrix<T>& a, int nrow, int ncol){
    Sparsity sp = a.sparsity().reshape(nrow,ncol);
    return Matrix<T>(sp,a.data());
  }

  template<class T>
  Matrix<T> reshape(const Matrix<T>& a, std::pair<int,int> rc){
    return reshape(a,rc.first,rc.second);
  }

  template<class T>
  Matrix<T> reshape(const Matrix<T>& x, const Sparsity& sp){
    // quick return if already the right shape
    if(sp==x.sparsity())
      return x;
  
    // make sure that the number of zeros agree
    casadi_assert(x.size()==sp.size());
  
    return Matrix<T>(sp,x.data());
  }

  template<class T>
  T trace(const Matrix<T>& a){
    casadi_assert_message(a.size2() == a.size1(), "trace: must be square");
    T res=0;
    for (int i=0; i< a.size2(); i ++) {
      res+=a.elem(i,i);
    }
    return res;
  }

  template<class T>
  Matrix<T> vec(const Matrix<T>& a){
    Matrix<T> ret = reshape(a,a.numel(),1);
    return ret;
  }

  template<class T>
  Matrix<T> vecNZ(const Matrix<T>& a){
    return Matrix<T>(a.data());
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
  Matrix<T> horzcat(const std::vector<Matrix<T> > &v){
    Matrix<T> ret;
    for(int i=0; i<v.size(); ++i)
      ret.appendColumns(v[i]);
    return ret;
  }

  template<class T>
  std::vector<Matrix<T> > horzsplit(const Matrix<T> &v, const std::vector<int>& offset) {
    // Consistency check
    casadi_assert(offset.size()>=1);
    casadi_assert(offset.front()==0);
    casadi_assert_message(offset.back()<=v.size2(),"horzsplit(const Matrix<T> &v, const std::vector<int>& offset): Last elements of offset (" << offset.back() << ") must be at maximum the number of cols in v (" << v.size2() << ")");
    casadi_assert(isMonotone(offset));
  
    std::vector<Matrix<T> > ret;
  
    // Obtain sparsity pattern
    const std::vector<int> & colind = v.sparsity().colind();
    const std::vector<int> & row = v.sparsity().row();
  
    for(int i=0; i<offset.size(); ++i) {
      int start = offset[i];
      int stop = i+1 < offset.size() ? offset[i+1] : v.size2(); 
  
      // colind for the submatrix: a portion of the original colind, 
      // but with a common offset substracted such that colind_s[0]==0
      std::vector<int> colind_s(stop-start+1,-colind[start]);
      std::transform(colind.begin()+start,colind.begin()+stop+1,colind_s.begin(),colind_s.begin(),std::plus<int>());
    
      // row for the submatrix: a portion of the original row
      std::vector<int> row_s(colind[stop]-colind[start]);
      std::copy(row.begin()+colind[start],row.begin()+colind[stop],row_s.begin());
    
      Sparsity s = Sparsity(v.size1(),stop-start,colind_s,row_s);
      Matrix<T> r(s);
    
      // data for the submatrix: a portion of the original data
      std::copy(v.begin()+colind[start],v.begin()+colind[stop],r.begin());
    
      // Append submatrix to list
      ret.push_back(r);
    }
    return ret;
  }

  template<class T>
  std::vector<Matrix<T> > horzsplit(const Matrix<T> &v, int incr) {
    casadi_assert(incr>=1);
    return horzsplit(v,range(0,v.size2(),incr));
  }


  template<class T>
  Matrix<T> vertcat(const std::vector<Matrix<T> > &v){
    Matrix<T> ret;
    for(int i=0; i<v.size(); ++i)
      ret.appendColumns(trans(v[i]));
    return trans(ret);  
  }

  template<class T>
  std::vector< Matrix<T> > vertsplit(const Matrix<T>& x, const std::vector<int>& offset){
    std::vector< Matrix<T> > ret = horzsplit(trans(x),offset);
    Matrix<T> (*transT)(const Matrix<T>& x) = trans;
    std::transform(ret.begin(),ret.end(),ret.begin(),transT);
    return ret;
  }
  
  template<class T>
  std::vector< Matrix<T> > vertsplit(const Matrix<T>& x, int incr){
    casadi_assert(incr>=1);
    return vertsplit(x,range(0,x.size1(),incr));
  }

  template<class T>
  std::vector< std::vector< Matrix<T> > > blocksplit(const Matrix<T>& x, const std::vector<int>& vert_offset, const std::vector<int>& horz_offset) {
    std::vector< Matrix<T> > rows = vertsplit(x,vert_offset);
    std::vector< std::vector< Matrix<T> > > ret;
    for (int i=0;i<rows.size();++i) {
      ret.push_back(horzsplit(rows[i],horz_offset));
    }
    return ret;
  }

  template<class T>
  std::vector< std::vector< Matrix<T> > > blocksplit(const Matrix<T>& x, int vert_incr, int horz_incr) {
    casadi_assert(horz_incr>=1);
    casadi_assert(vert_incr>=1);
    return blocksplit(x,range(0,x.size1(),vert_incr),range(0,x.size2(),horz_incr));
  }

  template<class T>
  Matrix<T> horzcat(const Matrix<T> &x, const Matrix<T> &y){
    Matrix<T> xy = x;
    xy.appendColumns(y);
    return xy;
  }

  template<class T>
  Matrix<T> vertcat(const Matrix<T> &x, const Matrix<T> &y){
    return trans(horzcat(trans(x),trans(y)));
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
    return sumAll(x*y);
  }

  template<class T>
  Matrix<T> outer_prod(const Matrix<T> &x, const Matrix<T> &y){
    return mul(x,trans(y));  
  }

  template<class T>
  Matrix<T> sumAll(const Matrix<T> &x) {
    // Quick return if empty
    if (x.empty()) return Matrix<T>::sparse(1,1);
    // Sum non-zero elements
    T res=0;
    for(int k=0; k<x.size(); k++){
      res += x.data()[k];
    }
    return res;
  }

  template<class T>
  Matrix<T> sumCols(const Matrix<T> &x) {
    return mul(x,Matrix<T>::ones(x.size2(),1));
  }

  template<class T>
  Matrix<T> sumRows(const Matrix<T> &x) {
    return mul(Matrix<T>::ones(1,x.size1()),x);
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
    if(x.vector()){
      return norm_F(x);
    } else {
      casadi_error("2-norms currently only supported for vectors. Did you intend to calculate a Frobenius norms (norm_F)?");
    }
  }

  template<class T>
  Matrix<T> norm_F(const Matrix<T>& x){
    return sqrt(1.0*sumAll(x*x));
  }

  template<class T>
  Matrix<T> norm_inf(const Matrix<T>& x){
    // Get largest element by absolute value
    T s = 0;
    for(typename std::vector<T>::const_iterator i=x.begin(); i!=x.end(); ++i){
      s = fmax(s,T(abs(*i)));
    }
    
    return s;
  }

  template<class T>
  void qr(const Matrix<T>& A, Matrix<T>& Q, Matrix<T> &R){
    // The following algorithm is taken from J. Demmel: Applied Numerical Linear Algebra (algorithm 3.1.)
    casadi_assert_message(A.size1()>=A.size2(), "qr: fewer rows than columns");

    // compute Q and R column by column
    Q = R = Matrix<T>();
    for(int i=0; i<A.size2(); ++i){
      // Initialize qi to be the i-th column of A
      Matrix<T> ai = A(ALL,i);
      Matrix<T> qi = ai;
      // The i-th column of R
      Matrix<T> ri = Matrix<T>::sparse(A.size2(),1);
  
      // subtract the projection of qi in the previous directions from ai
      for(int j=0; j<i; ++j){
      
        // Get the j-th column of Q
        Matrix<T> qj = Q(ALL,j);

        ri(j,0) = mul(trans(qi),qj); // Modified Gram-Schmidt
        // ri[j] = inner_prod(qj,ai); // Classical Gram-Schmidt
     
        // Remove projection in direction j
        if (ri.hasNZ(j,0))
          qi -= ri(j,0) * qj;
      }

      // Normalize qi
      ri(i,0) = norm_2(qi);
      qi /= ri(i,0);

      // Update R and Q
      Q.appendColumns(qi);
      R.appendColumns(ri);
    }
  }
  
  template<class T>
  Matrix<T> nullspace(const Matrix<T>& A) {
    int n = A.size1();
    int m = A.size2();
    
    Matrix<T> X = A;
    
    casadi_assert_message(m>=n,"nullspace(A): expecting a flat matrix (more columns than rows), but got " << A.dimString() << ".");
    
    Matrix<T> seed = DMatrix::eye(m)(Slice(0,m),Slice(n,m));

    std::vector< Matrix<T> > us;
    std::vector< Matrix<T> > betas;
    
    Matrix<T> beta;
    
    for (int i=0;i<n;++i) {
      Matrix<T> x = X(i,Slice(i,m));
      Matrix<T> u = Matrix<T>(x);
      Matrix<T> sigma = sqrt(sumCols(x*x));
      const Matrix<T>& x0 = x(0,0);
      u(0,0) = 1;
      
      Matrix<T> b = -copysign(sigma,x0);
      
      u(Slice(0),Slice(1,m-i))*= 1/(x0-b);
      beta = 1-x0/b;
      
      X(Slice(i,n),Slice(i,m))-= beta*mul(mul(X(Slice(i,n),Slice(i,m)),trans(u)),u);
      us.push_back(u);
      betas.push_back(beta);
    }
    
    for (int i=n-1;i>=0;--i) {
      seed(Slice(i,m),Slice(0,m-n)) -= betas[i]*mul(trans(us[i]),mul(us[i],seed(Slice(i,m),Slice(0,m-n))));
    }
    
    return seed;

  }

  template<class T>
  Matrix<T> solve(const Matrix<T>& A, const Matrix<T>& b){
    // check dimensions
    casadi_assert_message(A.size1() == b.size1(),"solve Ax=b: dimension mismatch: b has " << b.size1() << " rows while A has " << A.size1() << ".");
    casadi_assert_message(A.size1() == A.size2(),"solve: A not square but " << A.dimString());

    if(isTril(A)){
      // forward substitution if lower triangular
      Matrix<T> x = b;
      const std::vector<int> & Arow = A.row();
      const std::vector<int> & Acolind = A.colind();
      const std::vector<T> & Adata = A.data();
      for(int i=0; i<A.size2(); ++i){ // loop over columns forwards
        for(int k=0; k<b.size2(); ++k){ // for every right hand side
          if(!x.hasNZ(i,k)) continue;
          x(i,k) /= A(i,i);
          for(int kk=Acolind[i+1]-1; kk>=Acolind[i] && Arow[kk]>i; --kk){
            int j = Arow[kk]; 
            x(j,k) -= Adata[kk]*x(i,k);
          }
        }
      }
      return x;
    } else if(isTriu(A)){
      // backward substitution if upper triangular
      Matrix<T> x = b;
      const std::vector<int> & Arow = A.row();
      const std::vector<int> & Acolind = A.colind();
      const std::vector<T> & Adata = A.data();
      for(int i=A.size2()-1; i>=0; --i){ // loop over columns backwards
        for(int k=0; k<b.size2(); ++k){ // for every right hand side
          if(!x.hasNZ(i,k)) continue;
          x(i,k) /= A(i,i);
          for(int kk=Acolind[i]; kk<Acolind[i+1] && Arow[kk]<i; ++kk){ 
            int j = Arow[kk];
            x(j,k) -= Adata[kk]*x(i,k);
          }
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

      // Permute the right hand side
      Matrix<T> bperm = b(rowperm,ALL);

      // Permute the linear system
      Matrix<T> Aperm = A(rowperm,colperm);

      // Solution
      Matrix<T> xperm;

      // Solve permuted system
      if(isTril(Aperm)){
      
        // Forward substitution if lower triangular
        xperm = solve(Aperm,bperm);
      
      } else if(A.size2()<=3){
      
        // Form inverse by minor expansion and multiply if very small (up to 3-by-3)
        xperm = mul(inv(Aperm),bperm);

      } else {
      
        // Make a QR factorization
        Matrix<T> Q,R;
        qr(Aperm,Q,R);

        // Solve the factorized system (note that solve will now be fast since it is triangular)
        xperm = solve(R,mul(trans(Q),bperm));
      }

      // get the inverted column permutation
      std::vector<int> inv_colperm(colperm.size());
      for(int k=0; k<colperm.size(); ++k)
        inv_colperm[colperm[k]] = k;

      // Permute back the solution and return
      Matrix<T> x = xperm(inv_colperm,ALL);
      return x;
    }
  }
  
  template<class T>
  Matrix<T> pinv(const Matrix<T>& A) {
    if (A.size2()>=A.size1()) {
      return trans(solve(mul(A,trans(A)),A));
    } else {
      return solve(mul(trans(A),A),trans(A));
    }
  }
  
  template<class T>
  Matrix<T> kron(const Matrix<T>& a, const Matrix<T>& b) {
    const Sparsity &a_sp = a.sparsity();
    Matrix<T> filler = Matrix<T>::sparse(b.shape());
    std::vector< std::vector< Matrix<T> > > blocks(a.size1(),std::vector< Matrix<T> >(a.size2(),filler));
    for (int i=0;i<a.size1();++i) {
      for (int j=0;j<a.size2();++j) {
        int k = a_sp.getNZ(i,j);
        if (k!=-1) {
          blocks[i][j] = a[k]*b;
        }
      }
    }
    return blockcat(blocks);
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
    for(int col=0; col<ex.size2(); ++col)
      {
        // Loop over the elements in the col
        for(int el=ex.colind(col); el<ex.colind(col+1); ++el){ // loop over the non-zero elements
          if(ex.row(el) > col) break; // break inner loop (only lower triangular part is used)
          nz++;
        }
      }
    return nz;
  }

  template<class T>
  bool isEqual(const Matrix<T>& ex1,const Matrix<T> &ex2){
    if ((nnz(ex1)!=0 || nnz(ex2)!=0) && (ex1.size2()!=ex2.size2() || ex1.size1()!=ex2.size1())) return false;
    Matrix<T> difference = ex1 - ex2;  
    return isZero(difference);
  }

  template<class T>
  Matrix<T> repmat(const Matrix<T> &A, int n, int m){
    // First concatenate horizontally
    Matrix<T> col = horzcat(std::vector<Matrix<T> >(m, A));
  
    // Then vertically
    return vertcat(std::vector<Matrix<T> >(n, col));
  }

  template<class T>
  Matrix<T> diag(const Matrix<T>&A){
    // Nonzero mapping
    std::vector<int> mapping;
    // Get the sparsity
    Sparsity sp = A.sparsity().diag(mapping);
  
    Matrix<T> ret = Matrix<T>(sp);
  
    for (int k=0;k<mapping.size();k++) ret[k] = A[mapping[k]];
    return ret;
  }

  /** \brief   Construct a matrix with given block on the diagonal */
  template<class T>
  Matrix<T> blkdiag(const std::vector< Matrix<T> > &A) {
    std::vector<T> data;
  
    std::vector<Sparsity> sp;
    for (int i=0;i<A.size();++i) {
      data.insert(data.end(),A[i].data().begin(),A[i].data().end());
      sp.push_back(A[i].sparsity());
    }
  
    return Matrix<T>(blkdiag(sp),data);
  
  }

  template<class T>
  void getSparseTriplet(const Matrix<T>& A, std::vector<int>& col, std::vector<int>& row){
    row = A.sparsity().row();
    col = A.sparsity().getCol();
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
    Sparsity sp = A.sparsity().patternUnion(B.sparsity(),mapping);
  
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
  void makeSparse(Matrix<T>& A, double tol){
    // Quick return if there are no structurally zero entries
    if(!hasNonStructuralZeros(A) && tol==0)
      return;
  
    // Start with a matrix with no cols
    Matrix<T> Asp = Matrix<T>::sparse(A.size1(),0);

    // Loop over the cols
    for(int i=0; i<A.size2(); ++i){
      // Resize the matrix to accomodate the col
      Asp.resize(A.size1(),i+1);
    
      // Loop over the existing, possible nonzeros
      for(int el=A.colind(i); el<A.colind(i+1); ++el){
      
        // If it is not known to be a zero
        if(!casadi_limits<T>::isAlmostZero(A.at(el),tol)){
        
          // Get the row
          int j=A.row(el);
      
          // Save to new, sparse matrix
          Asp(j,i) = A.at(el);
        }
      }
    }
  
    // Save to A
    A = Asp;
  }

  template<class T>
  Matrix<T> sparse(const Matrix<T>& A, double tol){
    Matrix<T> ret(A);
    makeSparse(ret,tol);
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
    int d1=A.size2(), d2=A.size1();
    const std::vector<int> &colind=A.colind();
    const std::vector<int> &row=A.row();
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
    for(int i=0; i<d1; ++i){ // loop over cols
      for(int el=colind[i]; el<colind[i+1]; ++el){ // loop over the non-zero elements
        int j=row[el];  // row
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
  Matrix<T> project(const Matrix<T>& A, const Sparsity& sparsity){
    // Check dimensions
    if(!(A.empty() && sparsity.numel()==0)){
      casadi_assert_message(A.size2()==sparsity.size2() && A.size1()==sparsity.size1(),
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
#define MTT_INST(T,function_name)                       \
  %template(function_name) CasADi::function_name < T >;

// Define template instanciations
#define MATRIX_TOOLS_TEMPLATES_COMMON(T)        \
  MTT_INST(T,trans)                             \
  MTT_INST(T,mul)                               \
  MTT_INST(T,isConstant)                        \
  MTT_INST(T,isDense)                           \
  MTT_INST(T,isEmpty)                           \
  MTT_INST(T,isInteger)                         \
  MTT_INST(T,isScalar)                          \
  MTT_INST(T,isVector)                          \
  MTT_INST(T,isTril)                            \
  MTT_INST(T,isTriu)                            \
  MTT_INST(T,det)                               \
  MTT_INST(T,getMinor)                          \
  MTT_INST(T,cofactor)                          \
  MTT_INST(T,adj)                               \
  MTT_INST(T,inv)                               \
  MTT_INST(T,reshape)                           \
  MTT_INST(T,vec)                           \
  MTT_INST(T,vecNZ)                         \
  MTT_INST(T,blockcat)                          \
  MTT_INST(T,blocksplit)                        \
  MTT_INST(T,vertcat)                           \
  MTT_INST(T,vertsplit)                         \
  MTT_INST(T,horzcat)                           \
  MTT_INST(T,horzsplit)                         \
  MTT_INST(T,inner_prod)                        \
  MTT_INST(T,outer_prod)                        \
  MTT_INST(T,norm_1)                            \
  MTT_INST(T,norm_2)                            \
  MTT_INST(T,norm_inf)                          \
  MTT_INST(T,norm_F)                            \
  MTT_INST(T,qr)                                \
  MTT_INST(T,nullspace)                         \
  MTT_INST(T,solve)                             \
  MTT_INST(T,pinv)                              \
  MTT_INST(T,isZero)                            \
  MTT_INST(T,isOne)                             \
  MTT_INST(T,isMinusOne)                        \
  MTT_INST(T,isIdentity)                        \
  MTT_INST(T,nnz)                               \
  MTT_INST(T,nnz_sym)                           \
  MTT_INST(T,isEqual)                           \
  MTT_INST(T,repmat)                            \
  MTT_INST(T,getSparseTriplet)                  \
  MTT_INST(T,unite)                             \
  MTT_INST(T,sumRows)                           \
  MTT_INST(T,sumCols)                           \
  MTT_INST(T,sumAll)                            \
  MTT_INST(T,trace)                             \
  MTT_INST(T,makeDense)                         \
  MTT_INST(T,densify)                           \
  MTT_INST(T,makeSparse)                        \
  MTT_INST(T,hasNonStructuralZeros)             \
  MTT_INST(T,diag)                              \
  MTT_INST(T,blkdiag)                           \
  MTT_INST(T,polyval)                           \
  MTT_INST(T,addMultiple)                       \
  MTT_INST(T,veccat)                        \
  MTT_INST(T,vecNZcat)                      \
  MTT_INST(T,project)                           \
  MTT_INST(T,sprank)                            \
  MTT_INST(T,kron) 
#endif //SWIG

#ifdef SWIGOCTAVE
#define MATRIX_TOOLS_TEMPLATES(T) MATRIX_TOOLS_TEMPLATES_COMMON(T)
#else
#define MATRIX_TOOLS_TEMPLATES(T)               \
  MATRIX_TOOLS_TEMPLATES_COMMON(T)              \
  MTT_INST(T,sparse)                            \
  MTT_INST(T,full)
#endif //SWIGOCTAVE

#endif // MATRIX_TOOLS_HPP
