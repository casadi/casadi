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
  /// \cond INTERNAL
  /** \brief  Simplify the expression: formulates the expression as and eliminates terms */
  inline void simplify(SXElement& ex) { ex = ex.zz_simplify();}
  /// \endcond

  /** \brief Evaluate an SX graph numerically
   * Note: this is not efficient. For critical parts (loops) of your code, always use SXFunction.
   */
  CASADI_EXPORT Matrix<double> evalf(const SX &ex);

  /** \brief Substitute variable v with value vdef in an expression ex, and evaluate numerically
   * Note: this is not efficient. For critical parts (loops) of your code, always use SXFunction.
   */
  CASADI_EXPORT Matrix<double> evalf(const SX &ex, const SX &v,
                                              const Matrix<double> &vdef);

#if !defined(SWIG) || !defined(SWIGMATLAB)

  /** \brief  Expand the expression as a weighted sum (with constant weights)
  */
  inline void expand(const SX& ex, SX &weights, SX& terms) { ex.zz_expand(weights, terms);}

  /** \brief Create a piecewise constant function
      Create a piecewise constant function with n=val.size() intervals

      Inputs:
      \param t a scalar variable (e.g. time)
      \param tval vector with the discrete values of t at the interval transitions (length n-1)
      \param val vector with the value of the function for each interval (length n)
  */
  inline SX pw_const(const SX &t, const SX &tval, const SX &val) { return t.zz_pw_const(tval, val);}

  /** Create a piecewise linear function
      Create a piecewise linear function:

      Inputs:
      \brief t a scalar variable (e.g. time)
      \brief tval vector with the the discrete values of t (monotonically increasing)
      \brief val vector with the corresponding function values (same length as tval)
  */
  inline SX pw_lin(const SX &t, const SX &tval, const SX &val) { return t.zz_pw_lin(tval, val);}

  inline SX if_else(const SX &cond, const SX &if_true, const SX &if_false) {
    return cond.zz_if_else(if_true, if_false);
  }
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
  inline SX heaviside(const SX &x) { return x.zz_heaviside();}

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
  inline SX rectangle(const SX &x) { return x.zz_rectangle();}

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
  inline SX triangle(const SX &x) { return x.zz_triangle();}

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
  inline SX ramp(const SX &x) { return x.zz_ramp();}

  /** \brief  Integrate f from a to b using Gaussian quadrature with n points */
  inline SX gauss_quadrature(const SX &f, const SX &x, const SX &a, const SX &b,
                             int order=5, const SX& w=SX()) {
    return f.zz_gauss_quadrature(x, a, b, order, w);
  }

  /** \brief  Simplify an expression */
  inline void simplify(SX &ex) { ex = ex.zz_simplify();}

  /** \brief  Substitute variable v with expression vdef in an expression ex */
  inline SX substitute(const SX& ex, const SX& v, const SX& vdef) {
    return ex.zz_substitute(v, vdef);
  }

  /** \brief  Substitute variable var with expression expr in multiple expressions */
  inline std::vector<SX> substitute(const std::vector<SX>& ex,
                                    const std::vector<SX>& v,
                                    const std::vector<SX>& vdef) {
    return SX::zz_substitute(ex, v, vdef);
  }

  /** \brief Substitute variable var out of or into an expression expr,
   *  with an arbitrary number of other expressions piggyback (vector version) */
  inline void substituteInPlace(const std::vector<SX>& v,
                                std::vector<SX>& vdef,
                                std::vector<SX>& ex,
                                bool reverse=false) {
    return SX::zz_substituteInPlace(v, vdef, ex, reverse);
  }

  /** \brief Substitute variable var out of or into an expression expr,
   *  with an arbitrary number of other expressions piggyback */
  inline void substituteInPlace(const SX& v, SX &vdef,
                                std::vector<SX>& ex, bool reverse=false) {
    std::vector<SX> v2(1, v);
    std::vector<SX> vdef2(1, vdef);
    substituteInPlace(v2, vdef2, ex, reverse);
    vdef = vdef2.front();
  }

  /** \brief Substitute variable var out of or into an expression expr */
  inline void substituteInPlace(const SX& v, SX &vdef, bool reverse=false) {
    // Empty vector
    std::vector<SX> ex;
    substituteInPlace(v, vdef, ex, reverse);
  }

  /** \brief  Get the sparsity pattern of a matrix */
  inline SX spy(const SX& A) { return A.zz_spy();}

  /** \brief Check if expression depends on the argument
    The argument must be symbolic
  */
  inline bool dependsOn(const SX& f, const SX &arg) { return f.zz_dependsOn(arg); }

  /** \brief Get all symbols contained in the supplied expression
   * Get all symbols on which the supplied expression depends
   * \see SXFunction::getFree()
   */
  inline SX getSymbols(const SX& e) { return e.zz_getSymbols().front();}

  /** \brief Get all the free variables in an expression */
  inline SX getFree(const SX& ex) { return getSymbols(ex);}

  ///@{
  /** \brief Calculate jacobian via source code transformation

      Uses casadi::SXFunction::jac
  */
  inline SX jacobian(const SX &ex, const SX &arg) { return ex.zz_jacobian(arg);}
  inline SX gradient(const SX &ex, const SX &arg) { return ex.zz_gradient(arg);}
  inline SX tangent(const SX &ex, const SX &arg) { return ex.zz_tangent(arg);}
  inline SX hessian(const SX &ex, const SX &arg) { return ex.zz_hessian(arg);}

  // Hessian and gradient:
  inline void hessian(const SX &ex, const SX &arg, SX &H, SX &g) {
    return ex.zz_hessian(arg, H, g);
  }
  ///@}

  /** \brief Calculate the Jacobian and multiply by a vector from the right
      This is equivalent to <tt>mul(jacobian(ex, arg), v)</tt> or
      <tt>mul(jacobian(ex, arg).T, v)</tt> for
      transpose_jacobian set to false and true respectively. If contrast to these
      expressions, it will use directional derivatives which is typically (but
      not necessarily) more efficient if the complete Jacobian is not needed and v has few rows.
  */
  inline SX jacobianTimesVector(const SX &ex, const SX &arg, const SX &v,
                                bool transpose_jacobian=false) {
    return ex.zz_jacobianTimesVector(arg, v, transpose_jacobian);
  }

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
  inline SX taylor(const SX& ex, const SX& x,
                   const SX& a=0, int order=1) { return ex.zz_taylor(x, a, order);}

  /**
   * \brief multivariate Taylor series expansion
   *
   * Do Taylor expansions until the aggregated order of a term is equal to 'order'.
   * The aggregated order of \f$x^n y^m\f$ equals \f$n+m\f$.
   *
   */
  inline SX mtaylor(const SX& ex, const SX& x, const SX& a, int order=1) {
    return ex.zz_mtaylor(x, a, order);
  }

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
  inline SX mtaylor(const SX& ex, const SX& x, const SX& a, int order,
                    const std::vector<int>& order_contributions) {
    return ex.zz_mtaylor(x, a, order, order_contributions);
  }

  /** \brief Count number of nodes */
  inline int countNodes(const SX& A) { return A.zz_countNodes();}

  /** \brief Get a string representation for a binary SX, using custom arguments */
  inline std::string getOperatorRepresentation(const SX& x,
                                               const std::vector<std::string>& args) {
    return x.zz_getOperatorRepresentation(args);
  }

  /** \brief Extract shared subexpressions from an set of expressions */
  inline void extractShared(std::vector<SX>& ex,
                            std::vector<SX>& v, std::vector<SX>& vdef,
                            const std::string& v_prefix="v_",
                            const std::string& v_suffix="") {
    SX::zz_extractShared(ex, v, vdef, v_prefix, v_suffix);
  }

  /** \brief Print compact, introducing new variables for shared subexpressions */
  inline void printCompact(const SX& ex, std::ostream &stream=CASADI_COUT) {
    ex.zz_printCompact(stream);
  }

  /** \brief extracts polynomial coefficients from an expression
   *
   * \param ex Scalar expression that represents a polynomial
   * \param x  Scalar symbol that the polynomial is build up with
   */
  inline SX poly_coeff(const SX& ex, const SX&x) { return ex.zz_poly_coeff(x);}

  /** \brief Attempts to find the roots of a polynomial
   *
   *  This will only work for polynomials up to order 3
   *  It is assumed that the roots are real.
   *
   */
  inline SX poly_roots(const SX& p) { return p.zz_poly_roots();}

  /** \brief Attempts to find the eigenvalues of a symbolic matrix
   *  This will only work for up to 3x3 matrices
   */
  inline SX eig_symbolic(const SX& m) { return m.zz_eig_symbolic();}

#endif // !defined(SWIG) || !defined(SWIGMATLAB)

/*
@}
*/

} // namespace casadi

#endif // CASADI_SX_TOOLS_HPP
