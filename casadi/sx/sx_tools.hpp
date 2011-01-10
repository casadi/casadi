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

#ifndef SX_TOOLS_HPP
#define SX_TOOLS_HPP

#include "sx_matrix.hpp"
#include "../matrix/matrix_tools.hpp"

#ifdef WITH_UBLAS
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
namespace ublas = boost::numeric::ublas;
#endif

/** NOTE: This file needs some cleaning up
    Basically everything which deals with SXMatrix should be turned into template functions and moved to matrix_tools.hpp
*/

namespace CasADi{

  #ifndef SWIG
/** \brief Make a vector/matrix of symbolic variables - dimension 0 */
void make_symbolic(SX& v, const std::string& name);

/** \brief Make a vector/matrix of symbolic variables - higher dimension recursively */
template<typename A>
void make_symbolic(std::vector< A >& v, const std::string& name){
  for(int i=0; i<v.size(); ++i){
    std::stringstream ss;
    ss << name << "_" << i;
    make_symbolic(v[i],ss.str());
  }
}
#endif

/** \brief Create a one-dimensional stl vector of length n with symbolic variables */
std::vector<SX> create_symbolic(const std::string& name, int n);

/** \brief Create a two-dimensional stl vector of length n-by-m with symbolic variables */
std::vector< std::vector<SX> > create_symbolic(const std::string& name, int n, int m);

/** \brief Create a three-dimensional stl vector of length n-by-m-by-p with symbolic variables */
std::vector< std::vector< std::vector< SX> > > create_symbolic(const std::string& name, int n, int m, int p);

  
/** \brief  Expand the expression as a weighted sum (with constant weights)  */
void expand(const SXMatrix& ex, SXMatrix &weights, SXMatrix& terms);

/** \brief  Simplify the expression: formulates the expression as and eliminates terms */
void simplify(SX& ex);

#ifndef SWIG
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
#endif
  
  
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
SXMatrix vertcat(const SXMatrix &x, const SXMatrix &y);
SXMatrix horzcat(const SXMatrix &x, const SXMatrix &y);
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

/** \brief  make a vector
  Reshapes/flattens the SXMatrix such that the shape becomes (expr.numel(),1)
 */
SXMatrix vec(const SXMatrix &expr); 

SXMatrix getRow(const SXMatrix &expr, int i, int ni=1, int ki=1);
SXMatrix getColumn(const SXMatrix &expr, int j, int nj=1, int kj=1);

#ifndef SWIG
/** \brief  Convert stl vector to expression */
template<typename T>
SXMatrix toSXMatrix(const std::vector<T> &v){
  SXMatrix ret;
  for(int i=0; i<v.size(); ++i) ret << SXMatrix(v[i]);
  return ret;
}
#endif

#ifndef SWIG
/** \brief  concatenate */
template<typename T>
SXMatrix vertcat(const std::vector<T> &v){
  SXMatrix ret;
  for(int i=0; i<v.size(); ++i) ret << SXMatrix(v[i]);
  return ret;
}
#endif
// SXMatrix vertcat(const std::vector<SX>& comp);
SXMatrix vertcat(const SXMatrix& a, const SXMatrix& b);

#ifndef SWIG
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
#endif

// For some mysterious reason, activating this function for SWIG will throw an 
// ImportError: /home/jg/programs/casadi/build/swig/python/casadi/_casadi.so: undefined symbol: _ZN6CasADi8containsERKNS_8SXMatrixES2_
// when importing casadi in python

#ifndef SWIG
/**
Returns true if at least one element in list contains the scalar e.
*/
bool contains(const SXMatrix &list, const SXMatrix &e);
#endif

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

#ifndef SWIG
// "operator?:" can not be overloaded
template<typename T>
T if_else(const SX& cond, const T& if_true, const T &if_false){
  return if_false + (if_true-if_false)*cond;
}
#endif

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

#ifndef SWIG
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
#endif

/** \brief  Get the number of nodes of the tree of a SXMatrix */
int numNodes(const SXMatrix& A);

/// Check dependency: very inefficient algorithm
bool dependsOn(const SXMatrix& f, const SXMatrix &arg);

//@{
/** \brief  check if the matrix has certain properties */
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
#ifndef SWIG
SXMatrix operator&&(const SXMatrix &a, const SXMatrix &b);
SXMatrix operator||(const SXMatrix &a, const SXMatrix &b);
SXMatrix operator!(const SXMatrix &a);
#endif
SXMatrix operator==(const SXMatrix &a, const SXMatrix &b);
SXMatrix operator!=(const SXMatrix &a, const SXMatrix &b);

/** \brief  Fill the matrix with the value val, make empty sparse if zero */
void fill(SXMatrix& mat, const SX& val);

  /** \brief  Perform operations by ID */
SXMatrix binary(int op, const SXMatrix &x, const SXMatrix &y);
SXMatrix unary(int op, const SXMatrix &x);
SXMatrix scalar_matrix(int op, const SX &x, const SXMatrix &y);
SXMatrix matrix_scalar(int op, const SXMatrix &x, const SX &y);
SXMatrix matrix_matrix(int op, const SXMatrix &x, const SXMatrix &y);


} // namespace CasADi

#ifndef SWIG

/** \brief  Global functions with c equivalents: The implementation and syntax mirrors the standard c functions in math.h */
namespace std{
#define SXMatrix CasADi::SXMatrix
SXMatrix sin(const SXMatrix &x);
SXMatrix cos(const SXMatrix &x);
SXMatrix tan(const SXMatrix &x);
SXMatrix atan(const SXMatrix &x);
SXMatrix asin(const SXMatrix &x);
SXMatrix acos(const SXMatrix &x);
SXMatrix exp(const SXMatrix &x); // natural expontial
SXMatrix log(const SXMatrix &x); // natuaral logarithm
SXMatrix pow(const SXMatrix &x, const SXMatrix &n); // power function
SXMatrix sqrt(const SXMatrix &x); // square root
SXMatrix fmin(const SXMatrix &x, const SXMatrix &y);
SXMatrix fmax(const SXMatrix &x, const SXMatrix &y);
SXMatrix floor(const SXMatrix &x);
SXMatrix ceil(const SXMatrix &x); 
SXMatrix erf(const SXMatrix &x); 
#undef SXMatrix
} // namespace std

#endif

#endif // SX_TOOLS_HPP
