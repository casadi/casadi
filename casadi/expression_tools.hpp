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

#ifndef EXPRESSION_TOOLS_HPP
#define EXPRESSION_TOOLS_HPP

#include "sx/sx_matrix.hpp"
#include "fx/sx_function.hpp"

#ifdef WITH_UBLAS
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
namespace ublas = boost::numeric::ublas;
#endif

namespace CasADi{

/** \brief  Expand the expression as a weighted sum (with constant weights)  */
void expand(const SXMatrix& ex, SXMatrix &weights, SXMatrix& terms);

/** \brief  Simplify the expression: formulates the expression as and eliminates terms */
void simplify(SX& ex);

/** \brief Make a vector/matrix of symbolic variables */
template <typename iter_type, typename func_type>
void make_symbolic(iter_type first, func_type last, const std::string& name){
  int count = 0;
  for(iter_type cur=first; cur!=last; ++cur){
    std::stringstream ss;
    ss << name << "_" << count++;
    *cur = SX(ss.str());
  }
}

  
  
/// Make the 2-norm of an SXMatrix
SXMatrix norm_2(const SXMatrix& x);

/// Matrix product iof two SXMatrix es
SXMatrix prod(const SXMatrix &x, const SXMatrix &y);
/** \brief Inner product of two vectors
	Equals
	\code
	trans(x)*y
	\endcode
	with x and y vectors
*/
SXMatrix inner_prod(const SXMatrix &x, const SXMatrix &y); // inner product
/** \brief Outer product of two vectors
	Equals
	\code
	x*trans(y)
	\endcode
	with x and y vectors
*/
SXMatrix outer_prod(const SXMatrix &x, const SXMatrix &y); // outer product

//@{
/** \brief  Concatenate */
SXMatrix& operator<<(SXMatrix& expr, const SXMatrix& add); // remove when C++0X becomes available
SXMatrix vertcat(const std::vector<SXMatrix> &v);
SXMatrix horzcat(const std::vector<SXMatrix> &v);
//@}

/** \brief  Matlab's linspace function */
SXMatrix linspace(const SXMatrix &a, const SXMatrix &b, int nsteps); 

//@{
/** \brief  Calculating the inverse of a (very small) matrix: for larger matrices, make a QR factorization and solve R*x = trans(Q) for x */
SX det(const SXMatrix &x);  // determinant
SXMatrix adj(const SXMatrix &x);  // adjungate
SXMatrix inv(const SXMatrix &x);  // inverse
SX getMinor(const SXMatrix &x, int i, int j); // the (i,j) minor
SX cofactor(const SXMatrix &x, int i, int j); // the (i,j) cofactor
//@}

/** \brief  Transpose */
SXMatrix trans(const SXMatrix &expr); // transpose

/** \brief  create an n-by-n identity matrix */
SXMatrix eye(int n); 

/** \brief  create a matrix with all infinities */
SXMatrix inf(int n=1,int m=1);

/** \brief  create a matrix with all ones */
SXMatrix ones(int n, int m=1);

/** \brief  create a matrix with all ones */
SXMatrix zeros(int n, int m=1);


/** \brief Create a piecewise constant function 
  Create a piecewise constant function with n=val.size() intervals

  Inputs:
  \param t a scalar variable (e.g. time)
  \param tval vector with the discrete values of t at the interval transitions (length n-1)
  \param val vector with the value of the function for each interval (length n)
*/
SXMatrix pw_const(const SXMatrix &t, const SXMatrix &tval, const SXMatrix &val);

/** Create a piecewise linear function 
  Create a piecewise linear function:

  Inputs:
  \brief t a scalar variable (e.g. time)
  \brief tval vector with the the discrete values of t (monotonically increasing)
  \brief val vector with the corresponding function values (same length as tval)
*/
SXMatrix pw_lin(const SX &t, const SXMatrix &tval, const SXMatrix &val);

SXMatrix if_else(const SXMatrix &cond, const SXMatrix &if_true, const SXMatrix &if_false);
/// Heaviside step function
SXMatrix heaviside(const SXMatrix &x); // heaviside step function

/// sign function
SXMatrix sign(const SXMatrix &x);     // sign function

/** \brief  Integrate f from a to b using Gaussian quadrature with n points */
SXMatrix gauss_quadrature(SXMatrix f, const SXMatrix &x, const SXMatrix &a, const SXMatrix &b, int order=5, const SXMatrix& w=SXMatrix());

/** \brief  make a vector */
SXMatrix vec(const SXMatrix &expr); 

/// ... = A(i:ki:i+ni,j:kj:j+nj)
void getSub(SXMatrix &res, const SXMatrix &expr, int i, int j=0, int ni=1, int nj=1, int ki=1, int kj=1); 
/** \brief  A(i:ki:i+ni,j:kj:j+nj) = expr */
void setSub(const SXMatrix &expr, SXMatrix &res, int i, int j=0);
/// ... = A(i:ki:i+ni,:)
void getRow(SXMatrix &res, const SXMatrix &expr, int i, int ni=1, int ki=1);
/** \brief  A(i:ki:i+ni,:) = expr */
void setRow(const SXMatrix& expr, SXMatrix &res, int i, int ni=1, int ki=1);
/// ... = A(:,j:kj:j+nj)
void getColumn(SXMatrix &res, const SXMatrix &expr, int j, int nj=1, int kj=1);
/** \brief  A(:,j:kj:j+nj) = expr */
void setColumn(const SXMatrix& expr, SXMatrix &res, int j, int nj=1, int kj=1);


SXMatrix getRow(const SXMatrix &expr, int i, int ni=1, int ki=1);
SXMatrix getColumn(const SXMatrix &expr, int j, int nj=1, int kj=1);

/** \brief  Convert stl vector to expression */
template<typename T>
SXMatrix toSXMatrix(const std::vector<T> &v){
  SXMatrix ret;
  for(int i=0; i<v.size(); ++i) ret << SXMatrix(v[i]);
  return ret;
}

/** \brief  concatenate */
template<typename T>
SXMatrix vertcat(const std::vector<T> &v){
  SXMatrix ret;
  for(int i=0; i<v.size(); ++i) ret << SXMatrix(v[i]);
  return ret;
}

// SXMatrix vertcat(const std::vector<SX>& comp);
SXMatrix vertcat(const SXMatrix& a, const SXMatrix& b);

/** \brief  the following functions ensures ublas compability */
#ifdef HAVE_UBLAS
template<typename T>
SXMatrix toSXMatrix(const ublas::vector<T> &v){
  SXMatrix ret;
  for(int i=0; i<v.size(); ++i) ret << SXMatrix(v[i]);
  return ret;
}

template<typename T>
SXMatrix toSXMatrix(const ublas::matrix<T> &v){
  SXMatrix ret(x.size1(),x.size2());
  ret.reserve(x.size1() * x.size2());
  for(int i=0; i< sz.n; ++i)
    for(int j=0; j< sz.m; ++j)
      ret(i,j) = x(i,j);  
  return ret;
}
#endif


/**
Returns true if at least one element in list contains the scalar e.
*/
bool contains(const SXMatrix &list, const SXMatrix &e);

/** \brief  Simplify an expression */
void simplify(SXMatrix &ex);
/// remove identical calculations
void compress(SXMatrix &ex, int level=5); 
/// substitute variable var with expression expr
void substitute(SXMatrix &ex, const SXMatrix &var, const SXMatrix &expr); 

/** \brief  Make the expression smooth by replacing non-smooth nodes with binary variables */
void makeSmooth(SXMatrix &ex, SXMatrix &bvar, SXMatrix &bexpr);

/** \brief  Substitute derivatives with variables */
/** \brief void replaceDerivatives(SXMatrix &ex, const SXMatrix &var, const SXMatrix &dvar); */

// "operator?:" can not be overloaded
template<typename T>
T if_else(const SX& cond, const T& if_true, const T &if_false){
  return if_false + (if_true-if_false)*cond;
}

/** \brief  QR factorization using the modified Gram-Schmidt algorithm */
/** \brief  More stable than the classical Gram-Schmidt, but may break down if the columns of A are nearly linearly dependent */
/** \brief  See J. Demmel: Applied Numerical Linear Algebra (algorithm 3.1.) */
void qr(const SXMatrix& A, SXMatrix& Q, SXMatrix &R);

/** \brief  Check if a matrix is lower triangular (complexity ~ A.size1()) */
bool isTril(const SXMatrix &A);

/** \brief  Check if a matrix is upper triangular (complexity ~ A.size1()) */
bool isTriu(const SXMatrix &A);

/** \brief  Solve a system of equations: A*x = b */
SXMatrix solve(const SXMatrix& A, const SXMatrix& b);

/** \brief  Get the sparsity pattern of a matrix */
SXMatrix spy(const SXMatrix& A);

/** \brief  Create a block matrix */
template<int n, int m>
SXMatrix blockmatrix(SXMatrix array[n][m]){
/** \brief  Return matrix */
  SXMatrix ret;

/** \brief  loop over rows */
  for(int i=0; i<n; ++i){
/** \brief  Create a row */
    SXMatrix row;
    
/** \brief  append components to the row */
    for(int j=0; j<m; ++j){
      row << array[i][j];
    }
    
/** \brief  append row to matrix */
    ret << trans(row);
  }

  return ret;
}

/** \brief  Create a block matrix (vector) */
template<int n>
SXMatrix blockmatrix(SXMatrix array[n]){
/** \brief  Return matrix */
  SXMatrix ret;

/** \brief  loop over rows */
  for(int i=0; i<n; ++i){
/** \brief  append components */
    ret << array[i];
  }

  return ret;
}

/** \brief  Get the number of nodes of the tree of a SXMatrix */
int numNodes(const SXMatrix& A);

/// Check dependency: very inefficient algorithm
bool dependsOn(const SXMatrix& f, const SXMatrix &arg);

//@{
/** \brief  check if the matrix has a certain properties */
bool isConstant(const SXMatrix& ex);
bool isSymbolic(const SXMatrix& ex);
bool isDense(const SXMatrix& ex);
bool isEmpty(const SXMatrix& ex);
bool isInteger(const SXMatrix& ex);
bool isScalar(const SXMatrix& ex);
bool isVector(const SXMatrix& ex);
bool isSmooth(const SXMatrix& ex);
//@}

/** \brief  check if two expressions are the same */
bool isEqual(const SXMatrix &ex1, const SXMatrix &ex2);

//@{
/** \brief  Automatic differentiation */
SXMatrix jacobian(const SXMatrix &ex, const SXMatrix &arg);
SXMatrix gradient(const SXMatrix &ex, const SXMatrix &arg);
SXMatrix hessian(const SXMatrix &ex, const SXMatrix &arg);
void hessian(const SXMatrix &ex, const SXMatrix &arg, SXMatrix &H, SXMatrix &g); // hessian and gradient
//@}

/** \brief  Obtain the values of a constant expression */
double getValue(const SXMatrix &ex, int i=0, int j=0);          // for constant expressions only
int getIntValue(const SXMatrix &ex, int i=0, int j=0);          // integer version
void getValue(const SXMatrix &ex, double *res); // for all constant expressions
void getIntValue(const SXMatrix &ex, int *res); // integer version
const std::string& getName(const SXMatrix &ex); // get the name (only for scalar variables)

/// Check if zero
bool isZero(const SXMatrix &ex);

/** \brief  number of non-zeros */
int nnz(const SXMatrix &ex) ;
/** \brief  number of non-zeros, assuming the matrix is diagonal */
int nnz_sym(const SXMatrix &ex) ;

/** \brief  To and from string */
std::istream& operator>>(std::istream &stream, SXMatrix &expr);

SXMatrix operator<=(const SXMatrix &a, const SXMatrix &b);
SXMatrix operator>=(const SXMatrix &a, const SXMatrix &b);
SXMatrix operator<(const SXMatrix &a, const SXMatrix &b);
SXMatrix operator>(const SXMatrix &a, const SXMatrix &b);
SXMatrix operator&&(const SXMatrix &a, const SXMatrix &b);
SXMatrix operator||(const SXMatrix &a, const SXMatrix &b);
SXMatrix operator==(const SXMatrix &a, const SXMatrix &b);
SXMatrix operator!=(const SXMatrix &a, const SXMatrix &b);
SXMatrix operator!(const SXMatrix &a);


} // namespace CasADi

#endif // EXPRESSION_TOOLS_HPP
