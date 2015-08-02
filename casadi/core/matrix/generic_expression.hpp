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


#ifndef CASADI_GENERIC_EXPRESSION_HPP
#define CASADI_GENERIC_EXPRESSION_HPP

#include "../casadi_math.hpp"

namespace casadi {

  /** \brief Empty Base
      This class is extended in SWIG.
   */
  struct CASADI_EXPORT GenericExpressionCommon {};

  /** \brief Expression interface
  *
  This is a common base class for SX, MX and Matrix<>, introducing a uniform syntax and implementing
  common functionality using the curiously recurring template pattern (CRTP) idiom.\n

  \author Joel Andersson
  \date 2012
*/
template<typename ExType>
class CASADI_EXPORT GenericExpression : public GenericExpressionCommon {
#ifndef SWIG
  protected:
    // Helper functions
    inline const ExType& self() const { return static_cast<const ExType&>(*this); }
    inline ExType& self() { return static_cast<ExType&>(*this); }
#endif // SWIG
  public:

#if !defined(SWIG) || defined(DOXYGEN)
/**
\ingroup expression_tools
@{
*/
    /// Addition
    friend inline ExType operator+(const ExType &x, const ExType &y) {
      return x.zz_plus(y);
    }

    /// Subtraction
    friend inline ExType operator-(const ExType &x, const ExType &y) {
      return x.zz_minus(y);
    }

    /// Elementwise multiplication
    friend inline ExType operator*(const ExType &x, const ExType &y) {
      return x.zz_times(y);
    }

    /// Elementwise division
    friend inline ExType operator/(const ExType &x, const ExType &y) {
      return x.zz_rdivide(y);
    }

    /// Logic less than
    friend inline ExType operator<(const ExType &x, const ExType &y) {
      return x.zz_lt(y);
    }

    /// Logic less or equal to
    friend inline ExType operator<=(const ExType &x, const ExType &y) {
      return x.zz_le(y);
    }

    /// Logic greater than
    friend inline ExType operator>(const ExType &x, const ExType &y) {
      return y.zz_lt(x);
    }

    /// Logic greater or equal to
    friend inline ExType operator>=(const ExType &x, const ExType &y) {
      return y.zz_le(x);
    }

    /// Logic equal to
    friend inline ExType operator==(const ExType &x, const ExType &y) {
      return x.zz_eq(y);
    }

    /// Logic not equal to
    friend inline ExType operator!=(const ExType &x, const ExType &y) {
      return x.zz_ne(y);
    }

    /** \brief Logical `and`
     * Returns (an expression evaluating to) 1 if both
     * expressions are nonzero and 0 otherwise
     */
    friend inline ExType operator&&(const ExType &x, const ExType &y) {
      return x.zz_and(y);
    }

    /** \brief  Logical `or`
     * returns (an expression evaluating to) 1 if at
     * least one expression is nonzero and 0 otherwise
     */
    friend inline ExType operator||(const ExType &x, const ExType &y) {
      return x.zz_or(y);
    }

    /// Absolute value
    friend inline ExType fabs(const ExType& x) {
      return x.zz_abs();
    }

    /// Absolute value
    friend inline ExType abs(const ExType& x) {
      return x.zz_abs();
    }

    /// Square root
    friend inline ExType sqrt(const ExType& x) {
      return x.zz_sqrt();
    }

    /// Sine
    friend inline ExType sin(const ExType& x) {
      return x.zz_sin();
    }

    /// Cosine
    friend inline ExType cos(const ExType& x) {
      return x.zz_cos();
    }

    /// Tangent
    friend inline ExType tan(const ExType& x) {
      return x.zz_tan();
    }

    /// Arc tangent
    friend inline ExType atan(const ExType& x) {
      return x.zz_atan();
    }

    /// Arc sine
    friend inline ExType asin(const ExType& x) {
      return x.zz_asin();
    }

    /// Arc cosine
    friend inline ExType acos(const ExType& x) {
      return x.zz_acos();
    }

    /// Hyperbolic tangent
    friend inline ExType tanh(const ExType& x) {
      return x.zz_tanh();
    }

    /// Hyperbolic sine
    friend inline ExType sinh(const ExType& x) {
      return x.zz_sinh();
    }

    /// Hyperbolic cosine
    friend inline ExType cosh(const ExType& x) {
      return x.zz_cosh();
    }

    /// Inverse hyperbolic tangent
    friend inline ExType atanh(const ExType& x) {
      return x.zz_atanh();
    }

    /// Inverse hyperbolic sine
    friend inline ExType asinh(const ExType& x) {
      return x.zz_asinh();
    }

    /// Inverse hyperbolic cosine
    friend inline ExType acosh(const ExType& x) {
      return x.zz_acosh();
    }

    /// Exponential function
    friend inline ExType exp(const ExType& x) {
      return x.zz_exp();
    }

    /// Natural logarithm
    friend inline ExType log(const ExType& x) {
      return x.zz_log();
    }

    /// Base-10 logarithm
    friend inline ExType log10(const ExType& x) {
      return x.zz_log10();
    }

    /// Round down to nearest integer
    friend inline ExType floor(const ExType& x) {
      return x.zz_floor();
    }

    /// Round up to nearest integer
    friend inline ExType ceil(const ExType& x) {
      return x.zz_ceil();
    }

    /// Error function
    friend inline ExType erf(const ExType& x) {
      return x.zz_erf();
    }

    /// Invers error function
    friend inline ExType erfinv(const ExType& x) {
      return x.zz_erfinv();
    }

    /** Sine function
        sign(x)   := -1 for x<0
        sign(x)   :=  1 for x>0,
        sign(0)   :=  0
        sign(NaN) :=  NaN
     */
    friend inline ExType sign(const ExType& x) {
      return x.zz_sign();
    }

    /// Elementwise power
    friend inline ExType pow(const ExType& x, const ExType& n) {
      return x.zz_power(n);
    }

    /// Remainder after division
    friend inline ExType fmod(const ExType& x, const ExType& y) {
      return x.zz_mod(y);
    }

    /// Two argument arc tangent
    friend inline ExType atan2(const ExType& x, const ExType& y) {
      return x.zz_atan2(y);
    }

    /// Smallest of two values
    friend inline ExType fmin(const ExType& x, const ExType& y) {
      return x.zz_min(y);
    }

    /// Largest of two values
    friend inline ExType fmax(const ExType& x, const ExType& y) {
      return x.zz_max(y);
    }

    /// Simplify an expression
    friend inline ExType simplify(const ExType &x) {
      return x.zz_simplify();
    }

    /** \brief Check if two nodes are equivalent up to a given depth.
     * Depth=0 checks if the expressions are identical, i.e. points to the same node.
     *
     * a = x*x
     * b = x*x
     *
     *  a.isEqual(b, 0)  will return false, but a.isEqual(b, 1) will return true
     */
    friend inline bool isEqual(const ExType& x, const ExType& y, int depth=0) {
      return x.zz_isEqual(y, depth);
    }

    friend inline bool iszero(const ExType& x) {
      return x.isZero();
    }

    /// Copy sign
    friend inline ExType copysign(const ExType& x, const ExType& n) {
      return x.zz_copysign(n);
    }

    /// Elementwise power with const power
    friend inline ExType constpow(const ExType& x, const ExType& n) {
      return x.zz_constpow(n);
    }

    /// In-place addition
    inline ExType& operator+=(const ExType &y) { return self() = self().zz_plus(y); }

    /// In-place subtraction
    inline ExType& operator-=(const ExType &y) { return self() = self().zz_minus(y); }

    /// In-place elementwise multiplication
    inline ExType& operator*=(const ExType &y) {return self() = self().zz_times(y);}

    /// In-place elementwise division
    inline ExType& operator/=(const ExType &y) {return self() = self().zz_rdivide(y);}

    /** \brief  Logical `not`
     * Returns (an expression evaluating to) 1 if
     * expression is zero and 0 otherwise
     */
    inline ExType operator!() const { return self().zz_not(); }

    /// Logical not, alternative syntax
    friend inline ExType logic_not(const ExType& x) {
      return !x;
    }

    /// Logical and, alternative syntax
    friend inline ExType logic_and(const ExType& x, const ExType& y) {
      return x && y;
    }

    /// Logical or, alterntive syntax
    friend inline ExType logic_or(const ExType& x, const ExType& y) {
      return x || y;
    }
/** @} */
#endif // SWIG

};

} // namespace casadi

#endif // CASADI_GENERIC_EXPRESSION_HPP
