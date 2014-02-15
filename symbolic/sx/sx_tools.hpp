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

  /** \brief Create an n-by-m matrix with symbolic variables */
  SXMatrix ssym(const std::string& name, int n=1, int m=1);

  /** \brief Create an n-by-m matrix with symbolic variables */
  SXMatrix ssym(const std::string& name, const std::pair<int,int> & nm); 

  /** \brief Create a vector of length p with with matrices with symbolic variables of given sparsity */
  std::vector<SXMatrix> ssym(const std::string& name, const CRSSparsity& sp, int p);

  /** \brief Create a vector of length p with n-by-m matrices with symbolic variables */
  std::vector<SXMatrix> ssym(const std::string& name, int n, int m, int p);

  /** \brief Create a vector of length r of vectors of length p with matrices with symbolic variables with given sparsity */
  std::vector<std::vector<SXMatrix> > ssym(const std::string& name, const CRSSparsity& sp, int p, int r);
  
  /** \brief Create a vector of length r of vectors of length p with n-by-m matrices with symbolic variables */
  std::vector<std::vector<SXMatrix> > ssym(const std::string& name, int n, int m, int p, int r);

  /** \brief Create an matrix with symbolic variables, given a sparsity pattern */
  SXMatrix ssym(const std::string& name, const CRSSparsity& sp);

  /** \brief Create a symbolic matrix out of a numeric one */
  SXMatrix ssym(const Matrix<double>& x);

//@}

/** \brief  Expand the expression as a weighted sum (with constant weights)  */
void expand(const SXMatrix& ex, SXMatrix &weights, SXMatrix& terms);

/** \brief  Simplify the expression: formulates the expression as and eliminates terms */
void simplify(SX& ex);

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
SXMatrix heaviside(const SXMatrix &x);

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
SXMatrix rectangle(const SXMatrix &x);

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
SXMatrix triangle(const SXMatrix &x);

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
SXMatrix ramp(const SXMatrix &x);

/** \brief  Integrate f from a to b using Gaussian quadrature with n points */
SXMatrix gauss_quadrature(SXMatrix f, const SXMatrix &x, const SXMatrix &a, const SXMatrix &b, int order=5, const SXMatrix& w=SXMatrix());

/** \brief  Simplify an expression */
void simplify(SXMatrix &ex);

/** \brief  Remove identical calculations */
void compress(SXMatrix &ex, int level=5); 

/** \brief  Substitute variable v with expression vdef in an expression ex */
SXMatrix substitute(const SXMatrix& ex, const SXMatrix& v, const SXMatrix& vdef);

/** \brief  Substitute variable var with expression expr in multiple expressions */
std::vector<SXMatrix> substitute(const std::vector<SXMatrix>& ex, const std::vector<SXMatrix>& v, const std::vector<SXMatrix>& vdef);

/** \brief Substitute variable var out of or into an expression expr */
void substituteInPlace(const SXMatrix& v, SXMatrix &vdef, bool reverse=false);

/** \brief Substitute variable var out of or into an expression expr, with an arbitrary number of other expressions piggyback */
void substituteInPlace(const SXMatrix& v, SXMatrix &vdef, std::vector<SXMatrix>& ex, bool reverse=false);

/** \brief Substitute variable var out of or into an expression expr, with an arbitrary number of other expressions piggyback (vector version) */
void substituteInPlace(const std::vector<SXMatrix>& v, std::vector<SXMatrix>& vdef, std::vector<SXMatrix>& ex, bool reverse=false);

/** \brief Evaluate an SX graph numerically
* Note: this is not efficient. For critical parts (loops) of your code, always use SXFunction.
*/
Matrix<double> evalf(const SXMatrix &ex);

/** \brief Substitute variable v with value vdef in an expression ex, and evaluate numerically
* Note: this is not efficient. For critical parts (loops) of your code, always use SXFunction.
* @see numSample1D
*/
Matrix<double> evalf(const SXMatrix &ex, const SXMatrix &v, const Matrix<double> &vdef);

// /** \brief  Make the expression smooth by replacing non-smooth nodes with binary variables */
//void makeSmooth(SXMatrix &ex, SXMatrix &bvar, SXMatrix &bexpr);

/** \brief  Substitute derivatives with variables */
/** \brief void replaceDerivatives(SXMatrix &ex, const SXMatrix &var, const SXMatrix &dvar); */

//{@
/// Checks if expression does not contain NaN or Inf
bool isRegular(const SX& ex);
bool isRegular(const SXMatrix& ex);
//@}

#ifndef SWIG
// "operator?:" can not be overloaded
template<typename T>
T if_else(const SX& cond, const T& if_true, const T &if_false){
  return if_false + (if_true-if_false)*cond;
}
#endif

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
      row.append(array[i][j]);
    }
    
/** \brief  append row to matrix */
    ret.append(trans(row));
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
    ret.append(array[i]);
  }

  return ret;
}

template<> inline
SXMatrix GenericMatrix<SXMatrix>::sym(const std::string& name, const CRSSparsity& sp){ return ssym(name,sp);}

#endif

/// Check dependency: very inefficient algorithm
bool dependsOn(const SXMatrix& f, const SXMatrix &arg);


/** \brief Get all symbols contained in the supplied expression
* Get all symbols on which the supplied expression depends
* \see SXFunction::getFree()
*/
std::vector<SX> getSymbols(const SXMatrix& e);

/** \brief  check if smooth */
bool isSmooth(const SXMatrix& ex);

/** \brief  check if symbolic (Dense)

 Sparse matrices invariable return false 
  */
bool isSymbolic(const SXMatrix& ex);

/** \brief  check if symbolic

 Sparse matrices can return true if all non-zero elements are symbolic
 */
bool isSymbolicSparse(const SXMatrix& ex);

//@{
/** \brief Calculate jacobian via source code transformation

Uses CasADi::SXFunction::jac
 */
SXMatrix jacobian(const SXMatrix &ex, const SXMatrix &arg);
SXMatrix gradient(const SXMatrix &ex, const SXMatrix &arg);
SXMatrix tangent(const SXMatrix &ex, const SXMatrix &arg);
SXMatrix hessian(const SXMatrix &ex, const SXMatrix &arg);
void hessian(const SXMatrix &ex, const SXMatrix &arg, SXMatrix &H, SXMatrix &g); // hessian and gradient
//@}

/** \brief Calculate the Jacobian and multiply by a vector from the left
    This is equivalent to mul(jacobian(ex,arg),v) or mul(jacobian(ex,arg).T,v) for transpose_jacobian set to false and
    true respectively. If contrast to these expressions, it will use directional derivatives which is typically (but
    not necessarily) more efficient if the complete Jacobian is not needed and v has few columns.
 */
SXMatrix jacobianTimesVector(const SXMatrix &ex, const SXMatrix &arg, const SXMatrix &v, bool transpose_jacobian=false);

/** \brief  Obtain the values of a constant expression */
double getValue(const SXMatrix &ex, int i=0, int j=0);          // for constant expressions only
int getIntValue(const SXMatrix &ex, int i=0, int j=0);          // integer version
void getValue(const SXMatrix &ex, double *res); // for all constant expressions
void getIntValue(const SXMatrix &ex, int *res); // integer version
const std::string& getName(const SXMatrix &ex); // get the name (only for scalar variables)

/** \brief  Fill the matrix with the value val, make empty sparse if zero */
void fill(SXMatrix& mat, const SX& val);

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
SXMatrix taylor(const SXMatrix& ex,const SXMatrix& x, const SXMatrix& a=casadi_limits<SX>::zero,int order=1);

/**
* \brief multivariate taylor series expansion
*
* Do taylor expansions until the aggregated order of a term is equal to 'order'.
* The aggregated order of \f$x^n y^m\f$ equals \f$n+m\f$.
*
*/
SXMatrix mtaylor(const SXMatrix& ex,const SXMatrix& x, const SXMatrix& a,int order=1);
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
SXMatrix mtaylor(const SXMatrix& ex,const SXMatrix& x, const SXMatrix& a,int order,const std::vector<int>&order_contributions);

/** \brief Count number of nodes */
int countNodes(const SXMatrix& A);

/** \brief Get a string representation for a binary SX, using custom arguments */
std::string getOperatorRepresentation(const SX& x, const std::vector<std::string>& args);

  /** \brief Get all the free variables in an expression */
  SXMatrix getFree(const SXMatrix& ex);

  /** \brief Extract shared subexpressions from an set of expressions */
  void extractShared(std::vector<SX>& ex, 
                     std::vector<SX>& v, std::vector<SX>& vdef, 
                     const std::string& v_prefix="v_", const std::string& v_suffix="");
  
  /** \brief Print compact, introducing new variables for shared subexpressions */
  void printCompact(const SXMatrix& ex, std::ostream &stream=std::cout);
  
/** \brief extracts polynomial coefficients from an expression
*
* \parameter ex Scalar expression that represents a polynomial
* \paramater x  Scalar symbol that th epolynomial is build up with
*/  
SXMatrix poly_coeff(const SXMatrix& ex, const SXMatrix&x);

/** \brief Attempts to find the roots of a polynomial
*
*  This will only work for polynomials up to order 3
*  It is assumed that the roots are real.
*  
*/  
SXMatrix poly_roots(const SXMatrix& p);

/** \brief Attempts to find the eigenvalues of a symbolic matrix
*  This will only work for up to 3x3 matrices
*/  
SXMatrix eig_symbolic(const SXMatrix& m);

} // namespace CasADi

#endif // SX_TOOLS_HPP
