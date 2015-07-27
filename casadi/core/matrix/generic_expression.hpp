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

  /** \brief Expression interface - non-tempalated base class 
      SWIG wrapping, to allow treating GenericExpression as a non-templated class in SWIG.
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


/*! inline friend bool isEqual(const ExType& x, const ExType& y, int depth=0)
 * \brief Check if two nodes are equivalent up to a given depth.
 *  Depth=0 checks if the expressions are identical, i.e. points to the same node.
 *
 *  a = x*x
 *  b = x*x
 *
 *  a.isEqual(b, 0)  will return false, but a.isEqual(b, 1) will return true
 */

#define GENERIC_EXPRESSION_FRIENDS_UNWRAPPED(DECL, M)   \
    DECL M operator+(const M &x, const M &y) {          \
      return x.zz_plus(y);                              \
    }                                                   \
    DECL M operator-(const M &x, const M &y) {          \
      return x.zz_minus(y);                             \
    }                                                   \
    DECL M operator*(const M &x, const M &y) {          \
      return x.zz_times(y);                             \
    }                                                   \
    DECL M operator/(const M &x, const M &y) {          \
      return x.zz_rdivide(y);                           \
    }                                                   \
    DECL M operator<(const M &x, const M &y) {          \
      return x.zz_lt(y);                                \
    }                                                   \
    DECL M operator<=(const M &x, const M &y) {         \
      return x.zz_le(y);                                \
    }                                                   \
    DECL M operator>(const M &x, const M &y) {          \
      return x.zz_gt(y);                                \
    }                                                   \
    DECL M operator>=(const M &x, const M &y) {         \
      return x.zz_ge(y);                                \
    }                                                   \
    DECL M operator==(const M &x, const M &y) {         \
      return x.zz_eq(y);                                \
    }                                                   \
    DECL M operator!=(const M &x, const M &y) {         \
      return x.zz_ne(y);                                \
    }                                                   \
    DECL M operator&&(const M &x, const M &y) {         \
      return x.zz_and(y);                               \
    }                                                   \
    DECL M operator||(const M &x, const M &y) {         \
      return x.zz_or(y);                                \
    }                                                   \
    DECL M abs(const M& x) {                            \
      return x.zz_abs();                                \
    }                                                   \
    DECL M fabs(const M& x) {                           \
      return x.zz_abs();                                \
    }                                                   \
    DECL M sqrt(const M& x) {                           \
      return x.zz_sqrt();                               \
    }                                                   \
    DECL M sin(const M& x) {                            \
      return x.zz_sin();                                \
    }                                                   \
    DECL M cos(const M& x) {                            \
      return x.zz_cos();                                \
    }                                                   \
    DECL M tan(const M& x) {                            \
      return x.zz_tan();                                \
    }                                                   \
    DECL M atan(const M& x) {                           \
      return x.zz_atan();                               \
    }                                                   \
    DECL M asin(const M& x) {                           \
      return x.zz_asin();                               \
    }                                                   \
    DECL M acos(const M& x) {                           \
      return x.zz_acos();                               \
    }                                                   \
    DECL M tanh(const M& x) {                           \
      return x.zz_tanh();                               \
    }                                                   \
    DECL M sinh(const M& x) {                           \
      return x.zz_sinh();                               \
    }                                                   \
    DECL M cosh(const M& x) {                           \
      return x.zz_cosh();                               \
    }                                                   \
    DECL M atanh(const M& x) {                          \
      return x.zz_atanh();                              \
    }                                                   \
    DECL M asinh(const M& x) {                          \
      return x.zz_asinh();                              \
    }                                                   \
    DECL M acosh(const M& x) {                          \
      return x.zz_acosh();                              \
    }                                                   \
    DECL M exp(const M& x) {                            \
      return x.zz_exp();                                \
    }                                                   \
    DECL M log(const M& x) {                            \
      return x.zz_log();                                \
    }                                                   \
    DECL M log10(const M& x) {                          \
      return x.zz_log10();                              \
    }                                                   \
    DECL M floor(const M& x) {                          \
      return x.zz_floor();                              \
    }                                                   \
    DECL M ceil(const M& x) {                           \
      return x.zz_ceil();                               \
    }                                                   \
    DECL M erf(const M& x) {                            \
      return x.zz_erf();                                \
    }                                                   \
    DECL M sign(const M& x) {                           \
      return x.zz_sign();                               \
    }                                                   \
    DECL M pow(const M& x, const M& n) {                \
      return x.zz_power(n);                             \
    }                                                   \
    DECL M fmod(const M& x, const M& y) {               \
      return x.zz_mod(y);                               \
    }                                                   \
    DECL M atan2(const M& x, const M& y) {              \
      return x.zz_atan2(y);                             \
    }                                                   \
    DECL M fmin(const M& x, const M& y) {               \
      return x.zz_min(y);                               \
    }                                                   \
    DECL M fmax(const M& x, const M& y) {               \
      return x.zz_max(y);                               \
    }                                                   \

#define GENERIC_EXPRESSION_FRIENDS(DECL, M)                     \
    DECL M simplify(const M &x) {                               \
      return x.zz_simplify();                                   \
    }                                                           \
    DECL bool isEqual(const M& x, const M& y, int depth=0) {    \
      return x.zz_isEqual(y, depth);                            \
    }                                                           \
    DECL bool iszero(const M& x) {                              \
      return x.isZero();                                        \
    }                                                           \

#ifndef SWIG
  GENERIC_EXPRESSION_FRIENDS_UNWRAPPED(inline friend, ExType)
  GENERIC_EXPRESSION_FRIENDS(inline friend, ExType)

    /// In-place addition
    inline ExType& operator+=(const ExType &y) { return self() = self().zz_plus(y); }

    /// In-place subtraction
    inline ExType& operator-=(const ExType &y) { return self() = self().zz_minus(y); }

    /// In-place elementwise multiplication
    inline ExType& operator*=(const ExType &y) {return self() = self().zz_times(y);}

    /// In-place elementwise division
    inline ExType& operator/=(const ExType &y) {return self() = self().zz_rdivide(y);}

    /// Logic not
    inline ExType operator!() const { return self().zz_not(); }

#endif // SWIG

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
