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
    inline ExType& operator+=(const ExType &y) {
      return static_cast<ExType&>(*this) = static_cast<ExType*>(this)->zz_plus(y);
    }

    /// In-place subtraction
    inline ExType& operator-=(const ExType &y) {
      return static_cast<ExType&>(*this) = static_cast<ExType*>(this)->zz_minus(y);
    }

    /// In-place elementwise multiplication
    inline ExType& operator*=(const ExType &y) {return static_cast<ExType&>(*this) =
            static_cast<ExType*>(this)->zz_times(y);}

    /// In-place elementwise division
    inline ExType& operator/=(const ExType &y) {return static_cast<ExType&>(*this) =
            static_cast<ExType*>(this)->zz_rdivide(y);}

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
    inline ExType operator!() const { return static_cast<const ExType &>(*this).zz_not(); }

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
    #endif // !SWIG || DOXYGEN
    /** @} */

    // \cond CLUTTER

    /// Matrix division from left
    inline ExType __mldivide__(const ExType& y) const
    { return y.__mrdivide__(static_cast<const ExType&>(*this));}

    /// No need to have both < and >
    inline ExType zz_gt(const ExType& y) const
    { return y.zz_lt(static_cast<const ExType&>(*this));}

    /// No need to have both <= and >=
    inline ExType zz_ge(const ExType& y) const
    { return y.zz_le(static_cast<const ExType&>(*this));}

    /// Division (with <tt>__future__.division</tt> in effect)
    inline ExType __truediv__(const ExType& y) const {return static_cast<const ExType&>(*this)/y;}
    /// \endcond
};


} // namespace casadi

#endif // CASADI_GENERIC_EXPRESSION_HPP
