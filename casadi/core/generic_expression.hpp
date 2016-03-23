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

#include "calculus.hpp"

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
      return ExType::binary(OP_ADD, x, y);
    }

    /// Subtraction
    friend inline ExType operator-(const ExType &x, const ExType &y) {
      return ExType::binary(OP_SUB, x, y);
    }

    /// Elementwise multiplication
    friend inline ExType operator*(const ExType &x, const ExType &y) {
      return ExType::binary(OP_MUL, x, y);
    }

    /// Elementwise division
    friend inline ExType operator/(const ExType &x, const ExType &y) {
      return ExType::binary(OP_DIV, x, y);
    }

    /// Logic less than
    friend inline ExType operator<(const ExType &x, const ExType &y) {
      return ExType::binary(OP_LT, x, y);
    }

    /// Logic less or equal to
    friend inline ExType operator<=(const ExType &x, const ExType &y) {
      return ExType::binary(OP_LE, x, y);
    }

    /// Logic greater than
    friend inline ExType operator>(const ExType &x, const ExType &y) {
      return ExType::binary(OP_LT, y, x);
    }

    /// Logic greater or equal to
    friend inline ExType operator>=(const ExType &x, const ExType &y) {
      return ExType::binary(OP_LE, y, x);
    }

    /// Logic equal to
    friend inline ExType operator==(const ExType &x, const ExType &y) {
      return ExType::binary(OP_EQ, x, y);
    }

    /// Logic not equal to
    friend inline ExType operator!=(const ExType &x, const ExType &y) {
      return ExType::binary(OP_NE, x, y);
    }

    /** \brief Logical `and`
     * Returns (an expression evaluating to) 1 if both
     * expressions are nonzero and 0 otherwise
     */
    friend inline ExType operator&&(const ExType &x, const ExType &y) {
      return ExType::binary(OP_AND, x, y);
    }

    /** \brief  Logical `or`
     * returns (an expression evaluating to) 1 if at
     * least one expression is nonzero and 0 otherwise
     */
    friend inline ExType operator||(const ExType &x, const ExType &y) {
      return ExType::binary(OP_OR, x, y);
    }

    /// Absolute value
    friend inline ExType fabs(const ExType& x) {
      return ExType::unary(OP_FABS, x);
    }

    /// Absolute value
    friend inline ExType abs(const ExType& x) {
      return ExType::unary(OP_FABS, x);
    }

    /// Square root
    friend inline ExType sqrt(const ExType& x) {
      return ExType::unary(OP_SQRT, x);
    }

    /// Square
    friend inline ExType sq(const ExType& x) {
      return ExType::unary(OP_SQ, x);
    }

    /// Sine
    friend inline ExType sin(const ExType& x) {
      return ExType::unary(OP_SIN, x);
    }

    /// Cosine
    friend inline ExType cos(const ExType& x) {
      return ExType::unary(OP_COS, x);
    }

    /// Tangent
    friend inline ExType tan(const ExType& x) {
      return ExType::unary(OP_TAN, x);
    }

    /// Arc tangent
    friend inline ExType atan(const ExType& x) {
      return ExType::unary(OP_ATAN, x);
    }

    /// Arc sine
    friend inline ExType asin(const ExType& x) {
      return ExType::unary(OP_ASIN, x);
    }

    /// Arc cosine
    friend inline ExType acos(const ExType& x) {
      return ExType::unary(OP_ACOS, x);
    }

    /// Hyperbolic tangent
    friend inline ExType tanh(const ExType& x) {
      return ExType::unary(OP_TANH, x);
    }

    /// Hyperbolic sine
    friend inline ExType sinh(const ExType& x) {
      return ExType::unary(OP_SINH, x);
    }

    /// Hyperbolic cosine
    friend inline ExType cosh(const ExType& x) {
      return ExType::unary(OP_COSH, x);
    }

    /// Inverse hyperbolic tangent
    friend inline ExType atanh(const ExType& x) {
      return ExType::unary(OP_ATANH, x);
    }

    /// Inverse hyperbolic sine
    friend inline ExType asinh(const ExType& x) {
      return ExType::unary(OP_ASINH, x);
    }

    /// Inverse hyperbolic cosine
    friend inline ExType acosh(const ExType& x) {
      return ExType::unary(OP_ACOSH, x);
    }

    /// Exponential function
    friend inline ExType exp(const ExType& x) {
      return ExType::unary(OP_EXP, x);
    }

    /// Natural logarithm
    friend inline ExType log(const ExType& x) {
      return ExType::unary(OP_LOG, x);
    }

    /// Base-10 logarithm
    friend inline ExType log10(const ExType& x) {
      return log(x)*(1/std::log(10.));
    }

    /// Round down to nearest integer
    friend inline ExType floor(const ExType& x) {
      return ExType::unary(OP_FLOOR, x);
    }

    /// Round up to nearest integer
    friend inline ExType ceil(const ExType& x) {
      return ExType::unary(OP_CEIL, x);
    }

    /// Error function
    friend inline ExType erf(const ExType& x) {
      return ExType::unary(OP_ERF, x);
    }

    /// Invers error function
    friend inline ExType erfinv(const ExType& x) {
      return ExType::unary(OP_ERFINV, x);
    }

    /** Sine function
        sign(x)   := -1 for x<0
        sign(x)   :=  1 for x>0,
        sign(0)   :=  0
        sign(NaN) :=  NaN
     */
    friend inline ExType sign(const ExType& x) {
      return ExType::unary(OP_SIGN, x);
    }

    /// Elementwise power
    friend inline ExType pow(const ExType& x, const ExType& y) {
      return ExType::binary(OP_POW, x, y);
    }

    /// Remainder after division
    friend inline ExType fmod(const ExType& x, const ExType& y) {
      return ExType::binary(OP_FMOD, x, y);
    }

    /// Two argument arc tangent
    friend inline ExType atan2(const ExType& x, const ExType& y) {
      return ExType::binary(OP_ATAN2, x, y);
    }

    /// Conditional assignment
    friend inline ExType if_else_zero(const ExType& x, const ExType& y) {
      return ExType::binary(OP_IF_ELSE_ZERO, x, y);
    }

    /// Smallest of two values
    friend inline ExType fmin(const ExType& x, const ExType& y) {
      return ExType::binary(OP_FMIN, x, y);
    }

    /// Largest of two values
    friend inline ExType fmax(const ExType& x, const ExType& y) {
      return ExType::binary(OP_FMAX, x, y);
    }

    /** \brief Check if two nodes are equivalent up to a given depth.
     * Depth=0 checks if the expressions are identical, i.e. points to the same node.
     *
     * a = x*x
     * b = x*x
     *
     *  a.is_equal(b, 0)  will return false, but a.is_equal(b, 1) will return true
     */
    friend inline bool is_equal(const ExType& x, const ExType& y, int depth=0) {
      return ExType::is_equal(x, y, depth);
    }

    /// Copy sign
    friend inline ExType copysign(const ExType& x, const ExType& y) {
      return ExType::binary(OP_COPYSIGN, x, y);
    }

    /// Elementwise power with const power
    friend inline ExType constpow(const ExType& x, const ExType& y) {
      return ExType::binary(OP_CONSTPOW, x, y);
    }

    /// Debug printing
    friend inline ExType printme(const ExType& x, const ExType& y) {
      return ExType::binary(OP_PRINTME, x, y);
    }

    /// In-place addition
    inline ExType& operator+=(const ExType &y) { return self() = self() + y; }

    /// In-place subtraction
    inline ExType& operator-=(const ExType &y) { return self() = self() - y; }

    /// In-place elementwise multiplication
    inline ExType& operator*=(const ExType &y) {return self() = self() * y;}

    /// In-place elementwise division
    inline ExType& operator/=(const ExType &y) {return self() = self() / y;}

    /** \brief  Logical `not`
     * Returns (an expression evaluating to) 1 if
     * expression is zero and 0 otherwise
     */
    inline ExType operator!() const {
      return ExType::unary(OP_NOT, self());
    }

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
