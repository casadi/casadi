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

  /** \brief Expression interface
  *
  This is a common base class for SX, MX and Matrix<>, introducing a uniform syntax and implementing
  common functionality using the curiously recurring template pattern (CRTP) idiom.\n

  \author Joel Andersson
  \date 2012
*/
template<typename ExType>
class CASADI_EXPORT GenericExpression {
#ifndef SWIG
  protected:
    // Helper functions
    inline const ExType& self() const { return static_cast<const ExType&>(*this); }
    inline ExType& self() { return static_cast<ExType&>(*this); }
#endif // SWIG
  public:

#ifndef SWIG
    /// Addition
    inline friend ExType operator+(const ExType &x, const ExType &y) { return x.zz_plus(y); }

    /// Subtraction
    inline friend ExType operator-(const ExType &x, const ExType &y) { return x.zz_minus(y); }

    /// Elementwise multiplication
    inline friend ExType operator*(const ExType &x, const ExType &y) { return x.zz_times(y); }

    /// Elementwise division
    inline friend ExType operator/(const ExType &x, const ExType &y) { return x.zz_rdivide(y); }

    /// In-place addition
    inline ExType& operator+=(const ExType &y) { return self() = self().zz_plus(y); }

    /// In-place subtraction
    inline ExType& operator-=(const ExType &y) { return self() = self().zz_minus(y); }

    /// In-place elementwise multiplication
    inline ExType& operator*=(const ExType &y) {return self() = self().zz_times(y);}

    /// In-place elementwise division
    inline ExType& operator/=(const ExType &y) {return self() = self().zz_rdivide(y);}

    /// Logic less than
    inline friend ExType operator<(const ExType &x, const ExType &y) { return x.zz_lt(y); }

    /// Logic less or equal to
    inline friend ExType operator<=(const ExType &x, const ExType &y) { return x.zz_le(y); }

    /// Logic greater than
    inline friend ExType operator>(const ExType &x, const ExType &y) { return x.zz_gt(y); }

    /// Logic greater or equal to
    inline friend ExType operator>=(const ExType &x, const ExType &y) { return x.zz_ge(y); }

    /// Logic equal to
    inline friend ExType operator==(const ExType &x, const ExType &y) { return x.zz_eq(y); }

    /// Logic not equal to
    inline friend ExType operator!=(const ExType &x, const ExType &y) { return x.zz_ne(y); }

    /// Logic not
    inline ExType operator!() const { return self().zz_not(); }

    /// Logic and
    inline friend ExType operator&&(const ExType &x, const ExType &y) { return x.zz_and(y); }

    /// Logic or
    inline friend ExType operator||(const ExType &x, const ExType &y) { return x.zz_or(y); }

    /** \brief  Simplify an expression */
    inline friend ExType simplify(const ExType &x) { return x.zz_simplify();}
    #endif // SWIG
    /**
    \ingroup expression_tools
    @{
    */
    #if !defined(SWIG) || defined(DOXYGEN)
    /** \brief Check if two nodes are equivalent up to a given depth.
     *  Depth=0 checks if the expressions are identical, i.e. points to the same node.
     *
     *  a = x*x
     *  b = x*x
     *
     *  a.isEqual(b, 0)  will return false, but a.isEqual(b, 1) will return true
     */
    inline friend bool isEqual(const ExType& x, const ExType& y, int depth=0) {
      return x.zz_isEqual(y, depth);
    }

    /** \brief  check if the matrix is 0 (note that false negative answers are possible) */
    inline friend bool iszero(const ExType& x) { return x.isZero();}

    /** \brief Absolute value, C++ syntax */
    inline friend ExType abs(const ExType& x) { return x.zz_abs();}

    /** \brief Absolute value, C syntax */
    inline friend ExType fabs(const ExType& x) { return x.zz_abs();}

    /** \brief Square root */
    inline friend ExType sqrt(const ExType& x) { return x.zz_sqrt();}

    /** \brief Sine */
    inline friend ExType sin(const ExType& x) { return x.zz_sin();}

    /** \brief Cosine */
    inline friend ExType cos(const ExType& x) { return x.zz_cos();}

    /** \brief Tangent */
    inline friend ExType tan(const ExType& x) { return x.zz_tan();}

    /** \brief Arc tangent */
    inline friend ExType atan(const ExType& x) { return x.zz_atan();}

    /** \brief Arc sine */
    inline friend ExType asin(const ExType& x) { return x.zz_asin();}

    /** \brief Arc cosine */
    inline friend ExType acos(const ExType& x) { return x.zz_acos();}

    /** \brief Hyperbolic tangent */
    inline friend ExType tanh(const ExType& x) { return x.zz_tanh();}

    /** \brief Hyperbolic sine */
    inline friend ExType sinh(const ExType& x) { return x.zz_sinh();}

    /** \brief Hyperbolic cosine */
    inline friend ExType cosh(const ExType& x) { return x.zz_cosh();}

    /** \brief Arc hyperbolic tangent */
    inline friend ExType atanh(const ExType& x) { return x.zz_atanh();}

    /** \brief Arc hyperbolic sine */
    inline friend ExType asinh(const ExType& x) { return x.zz_asinh();}

    /** \brief Arc hyperbolic cosine */
    inline friend ExType acosh(const ExType& x) { return x.zz_acosh();}

    /** \brief Natural exponential function (elementwise for matrix types) */
    inline friend ExType exp(const ExType& x) { return x.zz_exp();}

    /** \brief Natural logarithm */
    inline friend ExType log(const ExType& x) { return x.zz_log();}

    /** \brief 10-base logarithm */
    inline friend ExType log10(const ExType& x) { return x.zz_log10();}

    /** \brief Round down to nearest integer */
    inline friend ExType floor(const ExType& x) { return x.zz_floor();}

    /** \brief Round up to nearest integer */
    inline friend ExType ceil(const ExType& x) { return x.zz_ceil();}

    /** \brief Error function */
    inline friend ExType erf(const ExType& x) { return x.zz_erf();}

    /** \brief Sign function (note sign(nan) == nan, sign(0) == 0) */
    inline friend ExType sign(const ExType& x) { return x.zz_sign();}

    /** \brief Power (elementwise for matrix types) */
    inline friend ExType pow(const ExType& x, const ExType& n) { return x.zz_power(n);}

    /** \brief Modulo */
    inline friend ExType fmod(const ExType& x, const ExType& y) { return x.zz_mod(y);}

    /** \brief Arctan2 */
    inline friend ExType atan2(const ExType& x, const ExType& y) { return x.zz_atan2(y);}

    /** \brief Minimum of two values */
    inline friend ExType fmin(const ExType& x, const ExType& y) { return x.zz_min(y);}

    /** \brief Maximum of two values */
    inline friend ExType fmax(const ExType& x, const ExType& y) { return x.zz_max(y);}
    #endif // !SWIG || DOXYGEN
    /** @} */

    // \cond CLUTTER

    /// Matrix division from left
    inline ExType __mldivide__(const ExType& y) const
    { return y.__mrdivide__(self());}

    /// No need to have both < and >
    inline ExType zz_gt(const ExType& y) const
    { return y.zz_lt(self());}

    /// No need to have both <= and >=
    inline ExType zz_ge(const ExType& y) const
    { return y.zz_le(self());}

    /// Division (with <tt>__future__.division</tt> in effect)
    inline ExType __truediv__(const ExType& y) const {return self()/y;}
    /// \endcond
};


} // namespace casadi

#endif // CASADI_GENERIC_EXPRESSION_HPP
