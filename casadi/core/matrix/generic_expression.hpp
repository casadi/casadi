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

#ifndef SWIG
    friend inline ExType operator+(const ExType &x, const ExType &y) {
      return x.zz_plus(y);
    }

    friend inline ExType operator-(const ExType &x, const ExType &y) {
      return x.zz_minus(y);
    }

    friend inline ExType operator*(const ExType &x, const ExType &y) {
      return x.zz_times(y);
    }

    friend inline ExType operator/(const ExType &x, const ExType &y) {
      return x.zz_rdivide(y);
    }

    friend inline ExType operator<(const ExType &x, const ExType &y) {
      return x.zz_lt(y);
    }

    friend inline ExType operator<=(const ExType &x, const ExType &y) {
      return x.zz_le(y);
    }

    friend inline ExType operator>(const ExType &x, const ExType &y) {
      return x.zz_gt(y);
    }

    friend inline ExType operator>=(const ExType &x, const ExType &y) {
      return x.zz_ge(y);
    }

    friend inline ExType operator==(const ExType &x, const ExType &y) {
      return x.zz_eq(y);
    }

    friend inline ExType operator!=(const ExType &x, const ExType &y) {
      return x.zz_ne(y);
    }

    friend inline ExType operator&&(const ExType &x, const ExType &y) {
      return x.zz_and(y);
    }

    friend inline ExType operator||(const ExType &x, const ExType &y) {
      return x.zz_or(y);
    }

    friend inline ExType fabs(const ExType& x) {
      return x.zz_abs();
    }

    friend inline ExType abs(const ExType& x) {
      return x.zz_abs();
    }

    friend inline ExType sqrt(const ExType& x) {
      return x.zz_sqrt();
    }

    friend inline ExType sin(const ExType& x) {
      return x.zz_sin();
    }

    friend inline ExType cos(const ExType& x) {
      return x.zz_cos();
    }

    friend inline ExType tan(const ExType& x) {
      return x.zz_tan();
    }

    friend inline ExType atan(const ExType& x) {
      return x.zz_atan();
    }

    friend inline ExType asin(const ExType& x) {
      return x.zz_asin();
    }

    friend inline ExType acos(const ExType& x) {
      return x.zz_acos();
    }

    friend inline ExType tanh(const ExType& x) {
      return x.zz_tanh();
    }

    friend inline ExType sinh(const ExType& x) {
      return x.zz_sinh();
    }

    friend inline ExType cosh(const ExType& x) {
      return x.zz_cosh();
    }

    friend inline ExType atanh(const ExType& x) {
      return x.zz_atanh();
    }

    friend inline ExType asinh(const ExType& x) {
      return x.zz_asinh();
    }

    friend inline ExType acosh(const ExType& x) {
      return x.zz_acosh();
    }

    friend inline ExType exp(const ExType& x) {
      return x.zz_exp();
    }

    friend inline ExType log(const ExType& x) {
      return x.zz_log();
    }

    friend inline ExType log10(const ExType& x) {
      return x.zz_log10();
    }

    friend inline ExType floor(const ExType& x) {
      return x.zz_floor();
    }

    friend inline ExType ceil(const ExType& x) {
      return x.zz_ceil();
    }

    friend inline ExType erf(const ExType& x) {
      return x.zz_erf();
    }

    friend inline ExType sign(const ExType& x) {
      return x.zz_sign();
    }

    friend inline ExType pow(const ExType& x, const ExType& n) {
      return x.zz_power(n);
    }

    friend inline ExType fmod(const ExType& x, const ExType& y) {
      return x.zz_mod(y);
    }

    friend inline ExType atan2(const ExType& x, const ExType& y) {
      return x.zz_atan2(y);
    }

    friend inline ExType fmin(const ExType& x, const ExType& y) {
      return x.zz_min(y);
    }

    friend inline ExType fmax(const ExType& x, const ExType& y) {
      return x.zz_max(y);
    }

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
