/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010-2014 Joel Andersson, Joris Gillis, Moritz Diehl,
 *                            K.U. Leuven. All rights reserved.
 *    Copyright (C) 2011-2014 Greg Horn
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


#ifndef CASADI_SX_TOOLS_HPP
#define CASADI_SX_TOOLS_HPP

#include "sx_element.hpp"
#include "../matrix/matrix_tools.hpp"
#include "../matrix/generic_matrix_tools.hpp"
#include "../casadi_options.hpp"

/** \defgroup expression_tools Expression tools
* Functions for manipulating DMatrix, SX, MX or Sparsity
*
*/


namespace casadi {

/**
\ingroup expression_tools
@{
*/

  /** \brief  Expand the expression as a weighted sum (with constant weights)
  */
  CASADI_CORE_EXPORT void expand(const SX& ex, SX &weights, SX& terms);

  /// \cond INTERNAL
  /** \brief  Simplify the expression: formulates the expression as and eliminates terms */
  CASADI_CORE_EXPORT void simplify(SXElement& ex);
  /// \endcond

  /** \brief Create a piecewise constant function
      Create a piecewise constant function with n=val.size() intervals

      Inputs:
      \param t a scalar variable (e.g. time)
      \param tval vector with the discrete values of t at the interval transitions (length n-1)
      \param val vector with the value of the function for each interval (length n)
  */
  CASADI_CORE_EXPORT SX pw_const(const SX &t, const SX &tval, const SX &val);

  /** Create a piecewise linear function
      Create a piecewise linear function:

      Inputs:
      \brief t a scalar variable (e.g. time)
      \brief tval vector with the the discrete values of t (monotonically increasing)
      \brief val vector with the corresponding function values (same length as tval)
  */
  CASADI_CORE_EXPORT SX pw_lin(const SXElement &t, const SX &tval, const SX &val);

  CASADI_CORE_EXPORT SX if_else(const SX &cond, const SX &if_true, const SX &if_false);
  /**  \brief Heaviside function
   *
   * \f[
   * \begin {cases}
   * H(x) = 0 & x<0 \\
   * H(x) = 1/2 & x=0 \\
   * H(x) = 1 & x>0 \\
   * \end {cases}
   * \f]
   */
  CASADI_CORE_EXPORT SX heaviside(const SX &x);

  /**
   * \brief rectangle function
   *
   * \f[
   * \begin {cases}
   * \Pi(x) = 1     & |x| < 1/2 \\
   * \Pi(x) = 1/2   & |x| = 1/2  \\
   * \Pi(x) = 0     & |x| > 1/2  \\
   * \end {cases}
   * \f]
   *
   * Also called: gate function, block function, band function, pulse function, window function
   */
  CASADI_CORE_EXPORT SX rectangle(const SX &x);

  /**
   * \brief triangle function
   *
   * \f[
   * \begin {cases}
   * \Lambda(x) = 0 &    |x| >= 1  \\
   * \Lambda(x) = 1-|x| &  |x| < 1
   * \end {cases}
   * \f]
   *
   */
  CASADI_CORE_EXPORT SX triangle(const SX &x);

  /**
   * \brief ramp function
   *
   *
   * \f[
   * \begin {cases}
   *  R(x) = 0   & x <= 1 \\
   *  R(x) = x   & x > 1 \\
   * \end {cases}
   * \f]
   *
   * Also called: slope function
   */
  CASADI_CORE_EXPORT SX ramp(const SX &x);

  /** \brief  Integrate f from a to b using Gaussian quadrature with n points */
  CASADI_CORE_EXPORT SX gauss_quadrature(SX f, const SX &x, const SX &a, const SX &b,
                                             int order=5, const SX& w=SX());

  /** \brief  Simplify an expression */
  CASADI_CORE_EXPORT void simplify(SX &ex);

  /** \brief  Remove identical calculations */
  CASADI_CORE_EXPORT void compress(SX &ex, int level=5);

  /** \brief  Substitute variable v with expression vdef in an expression ex */
  CASADI_CORE_EXPORT SX substitute(const SX& ex, const SX& v, const SX& vdef);

  /** \brief  Substitute variable var with expression expr in multiple expressions */
  CASADI_CORE_EXPORT std::vector<SX> substitute(const std::vector<SX>& ex,
                                                    const std::vector<SX>& v,
                                                    const std::vector<SX>& vdef);

  /** \brief Substitute variable var out of or into an expression expr */
  CASADI_CORE_EXPORT void substituteInPlace(const SX& v, SX &vdef, bool reverse=false);

  /** \brief Substitute variable var out of or into an expression expr,
   *  with an arbitrary number of other expressions piggyback */
  CASADI_CORE_EXPORT void substituteInPlace(const SX& v, SX &vdef,
                                                std::vector<SX>& ex, bool reverse=false);

  /** \brief Substitute variable var out of or into an expression expr,
   *  with an arbitrary number of other expressions piggyback (vector version) */
  CASADI_CORE_EXPORT void substituteInPlace(const std::vector<SX>& v,
                                                std::vector<SX>& vdef,
                                                std::vector<SX>& ex,
                                                bool reverse=false);

  /** \brief Evaluate an SX graph numerically
   * Note: this is not efficient. For critical parts (loops) of your code, always use SXFunction.
   */
  CASADI_CORE_EXPORT Matrix<double> evalf(const SX &ex);

  /** \brief Substitute variable v with value vdef in an expression ex, and evaluate numerically
   * Note: this is not efficient. For critical parts (loops) of your code, always use SXFunction.
   */
  CASADI_CORE_EXPORT Matrix<double> evalf(const SX &ex, const SX &v,
                                              const Matrix<double> &vdef);

  /** \brief  Get the sparsity pattern of a matrix */
  CASADI_CORE_EXPORT SX spy(const SX& A);

#ifndef SWIG
  // "operator?:" can not be overloaded
  template<typename T>
  T if_else(const SXElement& cond, const T& if_true, const T &if_false);

  /** \brief  Create a block matrix */
  template<int n, int m>
  SX blockmatrix(SX array[n][m]);

  /** \brief  Create a block matrix (vector) */
  template<int n>
  SX blockmatrix(SX array[n]);

#endif

// Implementation
#ifndef SWIG

  template<typename T>
  SX if_else(const SXElement& cond, const T& if_true, const T &if_false) {
    return if_false + (if_true-if_false)*cond;
  }

  template<int n, int m>
  SX blockmatrix(SX array[n][m]) {
    /** \brief  Return matrix */
    SX ret;

    /** \brief  loop over cols */
    for (int i=0; i<n; ++i) {
      /** \brief  Create a col */
      SX col;

      /** \brief  append components to the col */
      for (int j=0; j<m; ++j) {
        col.appendColumns(array[i][j]);
      }

      /** \brief  append col to matrix */
      ret.appendColumns(col.T());
    }

    return ret;
  }

  /** \brief  Create a block matrix (vector) */
  template<int n>
  SX blockmatrix(SX array[n]) {
    /** \brief  Return matrix */
    SX ret;

    /** \brief  loop over cols */
    for (int i=0; i<n; ++i) {
      /** \brief  append components */
      ret.appendColumns(array[i]);
    }

    return ret;
  }

#endif

  /** \brief Check if expression depends on the argument
    The argument must be symbolic
  */
  CASADI_CORE_EXPORT bool dependsOn(const SX& f, const SX &arg);


  /** \brief Get all symbols contained in the supplied expression
   * Get all symbols on which the supplied expression depends
   * \see SXFunction::getFree()
   */
  CASADI_CORE_EXPORT std::vector<SXElement> getSymbols(const SX& e);

  ///@{
  /** \brief Calculate jacobian via source code transformation

      Uses casadi::SXFunction::jac
  */
  CASADI_CORE_EXPORT SX jacobian(const SX &ex, const SX &arg);
  CASADI_CORE_EXPORT SX gradient(const SX &ex, const SX &arg);
  CASADI_CORE_EXPORT SX tangent(const SX &ex, const SX &arg);
  CASADI_CORE_EXPORT SX hessian(const SX &ex, const SX &arg);
  // hessian and gradient:
  CASADI_CORE_EXPORT void hessian(const SX &ex, const SX &arg, SX &H, SX &g);
  ///@}

  /** \brief Calculate the Jacobian and multiply by a vector from the right
      This is equivalent to <tt>mul(jacobian(ex, arg), v)</tt> or
      <tt>mul(jacobian(ex, arg).T, v)</tt> for
      transpose_jacobian set to false and true respectively. If contrast to these
      expressions, it will use directional derivatives which is typically (but
      not necessarily) more efficient if the complete Jacobian is not needed and v has few rows.
  */
  CASADI_CORE_EXPORT SX jacobianTimesVector(const SX &ex, const SX &arg, const SX &v,
                                                bool transpose_jacobian=false);

  /**
   * \brief univariate Taylor series expansion
   *
   * Calculate the Taylor expansion of expression 'ex' up to order 'order' with
   * respect to variable 'x' around the point 'a'
   *
   * \f$(x)=f(a)+f'(a)(x-a)+f''(a)\frac {(x-a)^2}{2!}+f'''(a)\frac{(x-a)^3}{3!}+\ldots\f$
   *
   * Example usage:
   * \code
   * taylor(sin(x), x)
   * \endcode
   * \verbatim >>   x \endverbatim
   */
  CASADI_CORE_EXPORT SX taylor(const SX& ex, const SX& x,
                                   const SX& a=casadi_limits<SXElement>::zero, int order=1);

  /**
   * \brief multivariate Taylor series expansion
   *
   * Do Taylor expansions until the aggregated order of a term is equal to 'order'.
   * The aggregated order of \f$x^n y^m\f$ equals \f$n+m\f$.
   *
   */
  CASADI_CORE_EXPORT SX mtaylor(const SX& ex, const SX& x, const SX& a, int order=1);
  /**
   * \brief multivariate Taylor series expansion
   *
   * Do Taylor expansions until the aggregated order of a term is equal to 'order'.
   * The aggregated order of \f$x^n y^m\f$ equals \f$n+m\f$.
   *
   * The argument order_contributions can denote how match each variable contributes
   * to the aggregated order. If x=[x, y] and order_contributions=[1, 2], then the
   * aggregated order of \f$x^n y^m\f$ equals \f$1n+2m\f$.
   *
   * Example usage
   *
   * \code
   * taylor(sin(x+y),[x, y],[a, b], 1)
   * \endcode
   * \f$ \sin(b+a)+\cos(b+a)(x-a)+\cos(b+a)(y-b) \f$
   * \code
   * taylor(sin(x+y),[x, y],[0, 0], 4)
   * \endcode
   * \f$  y+x-(x^3+3y x^2+3 y^2 x+y^3)/6  \f$
   * \code
   * taylor(sin(x+y),[x, y],[0, 0], 4,[1, 2])
   * \endcode
   * \f$  (-3 x^2 y-x^3)/6+y+x \f$
   *
   */
  CASADI_CORE_EXPORT SX mtaylor(const SX& ex, const SX& x, const SX& a, int order,
                                    const std::vector<int>&order_contributions);

  /** \brief Count number of nodes */
  CASADI_CORE_EXPORT int countNodes(const SX& A);

  /** \brief Get a string representation for a binary SX, using custom arguments */
  CASADI_CORE_EXPORT std::string getOperatorRepresentation(
                                       const SXElement& x,
                                       const std::vector<std::string>& args);

  /** \brief Get all the free variables in an expression */
  CASADI_CORE_EXPORT SX getFree(const SX& ex);

  /** \brief Extract shared subexpressions from an set of expressions */
  CASADI_CORE_EXPORT void extractShared(std::vector<SXElement>& ex,
                     std::vector<SXElement>& v, std::vector<SXElement>& vdef,
                     const std::string& v_prefix="v_", const std::string& v_suffix="");

  /** \brief Print compact, introducing new variables for shared subexpressions */
  CASADI_CORE_EXPORT void printCompact(const SX& ex, std::ostream &stream=std::cout);

  /** \brief extracts polynomial coefficients from an expression
   *
   * \parameter ex Scalar expression that represents a polynomial
   * \paramater x  Scalar symbol that the polynomial is build up with
   */
  CASADI_CORE_EXPORT SX poly_coeff(const SX& ex, const SX&x);

  /** \brief Attempts to find the roots of a polynomial
   *
   *  This will only work for polynomials up to order 3
   *  It is assumed that the roots are real.
   *
   */
  CASADI_CORE_EXPORT SX poly_roots(const SX& p);

  /** \brief Attempts to find the eigenvalues of a symbolic matrix
   *  This will only work for up to 3x3 matrices
   */
  CASADI_CORE_EXPORT SX eig_symbolic(const SX& m);

/*
@}
*/

} // namespace casadi



#endif // CASADI_SX_TOOLS_HPP
