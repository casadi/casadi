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


/*! inline friend bool isEqual(const ExType& x, const ExType& y, int depth=0)
 * \brief Check if two nodes are equivalent up to a given depth.
 *  Depth=0 checks if the expressions are identical, i.e. points to the same node.
 *
 *  a = x*x
 *  b = x*x
 *
 *  a.isEqual(b, 0)  will return false, but a.isEqual(b, 1) will return true
 */

#define GENERIC_EXPRESSION_FRIENDS_UNWRAPPED(M)                         \
    inline SWIG_FRIEND M operator+(const M &x, const M &y) {            \
      return x.zz_plus(y);                                              \
    }                                                                   \
    inline SWIG_FRIEND M operator-(const M &x, const M &y) {            \
      return x.zz_minus(y);                                             \
    }                                                                   \
    inline SWIG_FRIEND M operator*(const M &x, const M &y) {            \
      return x.zz_times(y);                                             \
    }                                                                   \
    inline SWIG_FRIEND M operator/(const M &x, const M &y) {            \
      return x.zz_rdivide(y);                                           \
    }                                                                   \
    inline SWIG_FRIEND M operator<(const M &x, const M &y) {            \
      return x.zz_lt(y);                                                \
    }                                                                   \
    inline SWIG_FRIEND M operator<=(const M &x, const M &y) {           \
      return x.zz_le(y);                                                \
    }                                                                   \
    inline SWIG_FRIEND M operator>(const M &x, const M &y) {            \
      return x.zz_gt(y);                                                \
    }                                                                   \
    inline SWIG_FRIEND M operator>=(const M &x, const M &y) {           \
      return x.zz_ge(y);                                                \
    }                                                                   \
    inline SWIG_FRIEND M operator==(const M &x, const M &y) {           \
      return x.zz_eq(y);                                                \
    }                                                                   \
    inline SWIG_FRIEND M operator!=(const M &x, const M &y) {           \
      return x.zz_ne(y);                                                \
    }                                                                   \
    inline SWIG_FRIEND M operator&&(const M &x, const M &y) {           \
      return x.zz_and(y);                                               \
    }                                                                   \
    inline SWIG_FRIEND M operator||(const M &x, const M &y) {           \
      return x.zz_or(y);                                                \
    }                                                                   \
    inline SWIG_FRIEND M abs(const M& x) {                              \
      return x.zz_abs();                                                \
    }                                                                   \
    inline SWIG_FRIEND M fabs(const M& x) {                             \
      return x.zz_abs();                                                \
    }                                                                   \
    inline SWIG_FRIEND M sqrt(const M& x) {                             \
      return x.zz_sqrt();                                               \
    }                                                                   \
    inline SWIG_FRIEND M sin(const M& x) {                              \
      return x.zz_sin();                                                \
    }                                                                   \
    inline SWIG_FRIEND M cos(const M& x) {                              \
      return x.zz_cos();                                                \
    }                                                                   \
    inline SWIG_FRIEND M tan(const M& x) {                              \
      return x.zz_tan();                                                \
    }                                                                   \
    inline SWIG_FRIEND M atan(const M& x) {                             \
      return x.zz_atan();                                               \
    }                                                                   \
    inline SWIG_FRIEND M asin(const M& x) {                             \
      return x.zz_asin();                                               \
    }                                                                   \
    inline SWIG_FRIEND M acos(const M& x) {                             \
      return x.zz_acos();                                               \
    }                                                                   \
    inline SWIG_FRIEND M tanh(const M& x) {                             \
      return x.zz_tanh();                                               \
    }                                                                   \
    inline SWIG_FRIEND M sinh(const M& x) {                             \
      return x.zz_sinh();                                               \
    }                                                                   \
    inline SWIG_FRIEND M cosh(const M& x) {                             \
      return x.zz_cosh();                                               \
    }                                                                   \
    inline SWIG_FRIEND M atanh(const M& x) {                            \
      return x.zz_atanh();                                              \
    }                                                                   \
    inline SWIG_FRIEND M asinh(const M& x) {                            \
      return x.zz_asinh();                                              \
    }                                                                   \
    inline SWIG_FRIEND M acosh(const M& x) {                            \
      return x.zz_acosh();                                              \
    }                                                                   \
    inline SWIG_FRIEND M exp(const M& x) {                              \
      return x.zz_exp();                                                \
    }                                                                   \
    inline SWIG_FRIEND M log(const M& x) {                              \
      return x.zz_log();                                                \
    }                                                                   \
    inline SWIG_FRIEND M log10(const M& x) {                            \
      return x.zz_log10();                                              \
    }                                                                   \
    inline SWIG_FRIEND M floor(const M& x) {                            \
      return x.zz_floor();                                              \
    }                                                                   \
    inline SWIG_FRIEND M ceil(const M& x) {                             \
      return x.zz_ceil();                                               \
    }                                                                   \
    inline SWIG_FRIEND M erf(const M& x) {                              \
      return x.zz_erf();                                                \
    }                                                                   \
    inline SWIG_FRIEND M sign(const M& x) {                             \
      return x.zz_sign();                                               \
    }                                                                   \
    inline SWIG_FRIEND M pow(const M& x, const M& n) {                  \
      return x.zz_power(n);                                             \
    }                                                                   \
    inline SWIG_FRIEND M fmod(const M& x, const M& y) {                 \
      return x.zz_mod(y);                                               \
    }                                                                   \
    inline SWIG_FRIEND M atan2(const M& x, const M& y) {                \
      return x.zz_atan2(y);                                             \
    }                                                                   \
    inline SWIG_FRIEND M fmin(const M& x, const M& y) {                 \
      return x.zz_min(y);                                               \
    }                                                                   \
    inline SWIG_FRIEND M fmax(const M& x, const M& y) {                 \
      return x.zz_max(y);                                               \
    }                                                                   \

#define GENERIC_EXPRESSION_FRIENDS(M)                                   \
    inline SWIG_FRIEND M simplify(const M &x) {                         \
      return x.zz_simplify();                                           \
    }                                                                   \
    inline SWIG_FRIEND bool isEqual(const M& x, const M& y, int depth=0) { \
      return x.zz_isEqual(y, depth);                                    \
    }                                                                   \
    inline SWIG_FRIEND bool iszero(const M& x) {                        \
      return x.isZero();                                                \
    }                                                                   \

#ifndef SWIG
  GENERIC_EXPRESSION_FRIENDS_UNWRAPPED(ExType)
  GENERIC_EXPRESSION_FRIENDS(ExType)

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
