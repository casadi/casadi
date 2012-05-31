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

#include "sx.hpp"
#include "../matrix/matrix_tools.hpp"
#include "../matrix/generic_matrix_tools.hpp"

#ifdef WITH_UBLAS
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
namespace ublas = boost::numeric::ublas;
#endif

/** Functions using SX */

namespace CasADi{

/** \brief  Construct symbolic arrays and variables using CasADi's SX expression graph representation
The "ssym" function is intended to work in a similar way as "sym" used in the Symbolic Toolbox for Matlab but instead creating an SXMatrix object.
The SX expression graph has much less overhead, but is also more restricted than the alternative MX expression graph.
*/
//@{

/** \brief  Construct symbolic arrays and variables using CasADi's more restricted, but more efficient SX expression graph
*/
//@{

  /** \brief Create an matrix with symbolic variables, with the dimension given by the string */
  Matrix<SX> ssym(const std::string& name);

  /** \brief Create an n-by-m matrix with symbolic variables */
  Matrix<SX> ssym(const std::string& name, int n, int m=1);

  /** \brief Create an n-by-m matrix with symbolic variables */
  Matrix<SX> ssym(const std::string& name, const std::pair<int,int> & nm); 

  /** \brief Create a vector of length p with with matrices with symbolic variables of given sparsity */
  std::vector<Matrix<SX> > ssym(const std::string& name, const CRSSparsity& sp, int p);

  /** \brief Create a vector of length p with n-by-m matrices with symbolic variables */
  std::vector<Matrix<SX> > ssym(const std::string& name, int n, int m, int p);

  /** \brief Create a vector of length r of vectors of length p with matrices with symbolic variables with given sparsity */
  std::vector<std::vector<Matrix<SX> > > ssym(const std::string& name, const CRSSparsity& sp, int p, int r);
  
  /** \brief Create a vector of length r of vectors of length p with n-by-m matrices with symbolic variables */
  std::vector<std::vector<Matrix<SX> > > ssym(const std::string& name, int n, int m, int p, int r);

  /** \brief Create an matrix with symbolic variables, given a sparsity pattern */
  Matrix<SX> ssym(const std::string& name, const CRSSparsity& sp);

  /** \brief Create a symbolic matrix out of a numeric one */
  Matrix<SX> ssym(const Matrix<double>& x);

//@}

/** \brief  Expand the expression as a weighted sum (with constant weights)  */
void expand(const Matrix<SX>& ex, Matrix<SX> &weights, Matrix<SX>& terms);

/** \brief  Simplify the expression: formulates the expression as and eliminates terms */
void simplify(SX& ex);

/** \brief Create a piecewise constant function 
  Create a piecewise constant function with n=val.size() intervals

  Inputs:
  \param t a scalar variable (e.g. time)
  \param tval vector with the discrete values of t at the interval transitions (length n-1)
  \param val vector with the value of the function for each interval (length n)
*/
Matrix<SX> pw_const(const Matrix<SX> &t, const Matrix<SX> &tval, const Matrix<SX> &val);

/** Create a piecewise linear function 
  Create a piecewise linear function:

  Inputs:
  \brief t a scalar variable (e.g. time)
  \brief tval vector with the the discrete values of t (monotonically increasing)
  \brief val vector with the corresponding function values (same length as tval)
*/
Matrix<SX> pw_lin(const SX &t, const Matrix<SX> &tval, const Matrix<SX> &val);

Matrix<SX> if_else(const Matrix<SX> &cond, const Matrix<SX> &if_true, const Matrix<SX> &if_false);
/**  \brief Heaviside function
*
* \f[
* \begin{cases}
* H(x) = 0 & x<0 \\
* H(x) = 1/2 & x=0 \\
* H(x) = 1 & x>0 \\
* \end{cases}
* \f]
*/
Matrix<SX> heaviside(const Matrix<SX> &x);

/** 
* \brief rectangle function
*
* \f[
* \begin{cases}
* \Pi(x) = 1     & |x| < 1/2 \\ 
* \Pi(x) = 1/2   & |x| = 1/2  \\
* \Pi(x) = 0     & |x| > 1/2  \\
* \end{cases}
* \f]
*
* Also called: gate function, block function, band function, pulse function, window function
*/
Matrix<SX> rectangle(const Matrix<SX> &x);

/** 
* \brief triangle function
*
* \f[
* \begin{cases}
* \Lambda(x) = 0 &    |x| >= 1  \\
* \Lambda(x) = 1-|x| &  |x| < 1 
* \end{cases}
* \f]
*
*/
Matrix<SX> triangle(const Matrix<SX> &x);

/** 
* \brief ramp function
*
* 
* \f[
* \begin{cases}
*  R(x) = 0   & x <= 1 \\
*  R(x) = x   & x > 1 \\
* \end{cases}
* \f]
*
* Also called: slope function
*/
Matrix<SX> ramp(const Matrix<SX> &x);

/** \brief  Integrate f from a to b using Gaussian quadrature with n points */
Matrix<SX> gauss_quadrature(Matrix<SX> f, const Matrix<SX> &x, const Matrix<SX> &a, const Matrix<SX> &b, int order=5, const Matrix<SX>& w=Matrix<SX>());

#ifndef SWIG
/**
Returns true if at least one element in list contains the scalar e.
*/
bool contains(const Matrix<SX> &list, const Matrix<SX> &e);
#endif

/** \brief  Simplify an expression */
void simplify(Matrix<SX> &ex);

/// Remove identical calculations
void compress(Matrix<SX> &ex, int level=5); 

/// Substitute variable v with expression vdef in an expression ex
Matrix<SX> substitute(const Matrix<SX> &ex, const Matrix<SX> &v, const Matrix<SX> &vdef);

/// Substitute variable var with expression expr in multiple expressions
std::vector<Matrix<SX> > substitute(const std::vector<Matrix<SX> > &ex, const Matrix<SX> &var, const Matrix<SX> &expr);

/// Substitute variable var out of or into an expression expr
void substituteInPlace(const Matrix<SX> &v, Matrix<SX> &vdef, bool reverse=false);

/// Substitute variable var out of or into an expression expr, with an arbitrary number of other expressions piggyback
void substituteInPlace(const Matrix<SX> &v, Matrix<SX> &vdef, std::vector<Matrix<SX> >& ex, bool reverse=false);

// /** \brief  Make the expression smooth by replacing non-smooth nodes with binary variables */
//void makeSmooth(Matrix<SX> &ex, Matrix<SX> &bvar, Matrix<SX> &bexpr);

/** \brief  Substitute derivatives with variables */
/** \brief void replaceDerivatives(Matrix<SX> &ex, const Matrix<SX> &var, const Matrix<SX> &dvar); */

#ifndef SWIG
// "operator?:" can not be overloaded
template<typename T>
T if_else(const SX& cond, const T& if_true, const T &if_false){
  return if_false + (if_true-if_false)*cond;
}
#endif

/** \brief  Get the sparsity pattern of a matrix */
Matrix<SX> spy(const Matrix<SX>& A);

#ifndef SWIG
/** \brief  Create a block matrix */
template<int n, int m>
Matrix<SX> blockmatrix(Matrix<SX> array[n][m]){
/** \brief  Return matrix */
  Matrix<SX> ret;

/** \brief  loop over rows */
  for(int i=0; i<n; ++i){
/** \brief  Create a row */
    Matrix<SX> row;
    
/** \brief  append components to the row */
    for(int j=0; j<m; ++j){
      row.append(array[i][j]);
    }
    
/** \brief  append row to matrix */
    ret.append(trans(row));
  }

  return ret;
}

/** \brief  Create a block matrix (vector) */
template<int n>
Matrix<SX> blockmatrix(Matrix<SX> array[n]){
/** \brief  Return matrix */
  Matrix<SX> ret;

/** \brief  loop over rows */
  for(int i=0; i<n; ++i){
/** \brief  append components */
    ret.append(array[i]);
  }

  return ret;
}
template<>
inline void sym(Matrix<SX>& ret, const std::string& name, int n, int m) {
  ret = ssym(name,n,m);
}

#endif

/// Check dependency: very inefficient algorithm
bool dependsOn(const Matrix<SX>& f, const Matrix<SX> &arg);

/** \brief  check if smooth */
bool isSmooth(const Matrix<SX>& ex);

/** \brief  check if symbolic (Dense)

 Sparse matrices invariable return false 
  */
bool isSymbolic(const Matrix<SX>& ex);

/** \brief  check if symbolic

 Sparse matrices can return true if all non-zero elements are symbolic
 */
bool isSymbolicSparse(const Matrix<SX>& ex);

//@{
/** \brief Calculate jacobian via source code transformation

Uses CasADi::SXFunction::jac
 */
Matrix<SX> jacobian(const Matrix<SX> &ex, const Matrix<SX> &arg);
Matrix<SX> gradient(const Matrix<SX> &ex, const Matrix<SX> &arg);
Matrix<SX> hessian(const Matrix<SX> &ex, const Matrix<SX> &arg);
void hessian(const Matrix<SX> &ex, const Matrix<SX> &arg, Matrix<SX> &H, Matrix<SX> &g); // hessian and gradient
//@}

/** \brief  Obtain the values of a constant expression */
double getValue(const Matrix<SX> &ex, int i=0, int j=0);          // for constant expressions only
int getIntValue(const Matrix<SX> &ex, int i=0, int j=0);          // integer version
void getValue(const Matrix<SX> &ex, double *res); // for all constant expressions
void getIntValue(const Matrix<SX> &ex, int *res); // integer version
const std::string& getName(const Matrix<SX> &ex); // get the name (only for scalar variables)

Matrix<SX> operator<=(const Matrix<SX> &a, const Matrix<SX> &b);
Matrix<SX> operator>=(const Matrix<SX> &a, const Matrix<SX> &b);
Matrix<SX> operator<(const Matrix<SX> &a, const Matrix<SX> &b);
Matrix<SX> operator>(const Matrix<SX> &a, const Matrix<SX> &b);
#ifndef SWIG
Matrix<SX> operator&&(const Matrix<SX> &a, const Matrix<SX> &b);
Matrix<SX> operator||(const Matrix<SX> &a, const Matrix<SX> &b);
Matrix<SX> operator!(const Matrix<SX> &a);
#endif
Matrix<SX> operator==(const Matrix<SX> &a, const Matrix<SX> &b);
Matrix<SX> operator!=(const Matrix<SX> &a, const Matrix<SX> &b);

/** \brief  Fill the matrix with the value val, make empty sparse if zero */
void fill(Matrix<SX>& mat, const SX& val);

/** 
* \brief univariate taylor series expansion
*
* Calculate the taylor expansion of expression 'ex' up to order 'order' with repsect to variable 'x' around the point 'a'
*
* \f$(x)=f(a)+f'(a)(x-a)+f''(a)\frac{(x-a)^2}{2!}+f'''(a)\frac{(x-a)^3}{3!}+\ldots\f$
* 
* Example usage:
* \code
* taylor(sin(x),x)
* \endcode
* \verbatim >>   x \endverbatim
*/
Matrix<SX> taylor(const Matrix<SX>& ex,const SX& x, const SX& a=casadi_limits<SX>::zero,int order=1);

/**
* \brief multivariate taylor series expansion
*
* Do taylor expansions until the aggregated order of a term is equal to 'order'.
* The aggregated order of \f$x^n y^m\f$ equals \f$n+m\f$.
*
*/
Matrix<SX> mtaylor(const Matrix<SX>& ex,const Matrix<SX>& x, const Matrix<SX>& a,int order=1);
/** 
* \brief multivariate taylor series expansion
*
* Do taylor expansions until the aggregated order of a term is equal to 'order'.
* The aggregated order of \f$x^n y^m\f$ equals \f$n+m\f$.
* 
* The argument order_contributions can denote how match each variable contributes to the aggregated order.
* If x=[x,y] and order_contributions=[1,2], then the aggregated order of \f$x^n y^m\f$ equals \f$1n+2m\f$.
*
* Example usage
*
* \code
* taylor(sin(x+y),[x,y],[a,b],1)
* \endcode
* \f$ \sin(b+a)+\cos(b+a)(x-a)+\cos(b+a)(y-b) \f$
* \code
* taylor(sin(x+y),[x,y],[0,0],4)
* \endcode
* \f$  y+x-(x^3+3y x^2+3 y^2 x+y^3)/6  \f$
* \code
* taylor(sin(x+y),[x,y],[0,0],4,[1,2])
* \endcode
* \f$  (-3 x^2 y-x^3)/6+y+x \f$
*
*/
Matrix<SX> mtaylor(const Matrix<SX>& ex,const Matrix<SX>& x, const Matrix<SX>& a,int order,const std::vector<int>&order_contributions);

/** \brief Count number of nodes */
int countNodes(const Matrix<SX>& A);

/** \brief Get a string representation for a binary SX, using custom arguments */
std::string getOperatorRepresentation(const SX& x, const std::vector<std::string>& args);

/** \brief Get all the free variables in an expression */
SXMatrix getFree(const SXMatrix& ex);

} // namespace CasADi

#endif // SX_TOOLS_HPP
