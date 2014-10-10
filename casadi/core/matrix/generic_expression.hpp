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
class CASADI_CORE_EXPORT GenericExpression {
  public:

#ifndef SWIG
    /// Addition
    inline friend ExType operator+(const ExType &x, const ExType &y) { return x.__add__(y); }

    /// Subtraction
    inline friend ExType operator-(const ExType &x, const ExType &y) { return x.__sub__(y); }

    /// Elementwise multiplication
    inline friend ExType operator*(const ExType &x, const ExType &y) { return x.__mul__(y); }

    /// Elementwise division
    inline friend ExType operator/(const ExType &x, const ExType &y) { return x.__div__(y); }

    /// In-place addition
    inline ExType& operator+=(const ExType &y) {
      return static_cast<ExType&>(*this) = static_cast<ExType*>(this)->__add__(y);
    }

    /// In-place subtraction
    inline ExType& operator-=(const ExType &y) {
      return static_cast<ExType&>(*this) = static_cast<ExType*>(this)->__sub__(y);
    }

    /// In-place elementwise multiplication
    inline ExType& operator*=(const ExType &y) {return static_cast<ExType&>(*this) =
            static_cast<ExType*>(this)->__mul__(y);}

    /// In-place elementwise division
    inline ExType& operator/=(const ExType &y) {return static_cast<ExType&>(*this) =
            static_cast<ExType*>(this)->__div__(y);}

    /// Logic less than
    inline friend ExType operator<(const ExType &x, const ExType &y) { return x.__lt__(y); }

    /// Logic less or equal to
    inline friend ExType operator<=(const ExType &x, const ExType &y) { return x.__le__(y); }

    /// Logic greater than
    inline friend ExType operator>(const ExType &x, const ExType &y) { return x.__gt__(y); }

    /// Logic greater or equal to
    inline friend ExType operator>=(const ExType &x, const ExType &y) { return x.__ge__(y); }

    /// Logic equal to
    inline friend ExType operator==(const ExType &x, const ExType &y) { return x.__eq__(y); }

    /// Logic not equal to
    inline friend ExType operator!=(const ExType &x, const ExType &y) { return x.__ne__(y); }

    /// Logic not
    inline ExType operator!() const { return static_cast<const ExType &>(*this).logic_not(); }

    /// Logic and
    inline friend ExType operator&&(const ExType &x, const ExType &y) { return x.logic_and(y); }

    /// Logic or
    inline friend ExType operator||(const ExType &x, const ExType &y) { return x.logic_or(y); }

    #endif // SWIG

    // \cond SWIGINTERNAL

    /// Matrix division from left
    inline ExType __mldivide__(const ExType& y) const
    { return y.__mrdivide__(static_cast<const ExType&>(*this));}

    /// No need to have both < and >
    inline ExType __gt__(const ExType& y) const
    { return y.__lt__(static_cast<const ExType&>(*this));}

    /// No need to have both <= and >=
    inline ExType __ge__(const ExType& y) const
    { return y.__le__(static_cast<const ExType&>(*this));}

    /// Division (with <tt>__future__.division</tt> in effect)
    inline ExType __truediv__(const ExType& y) const {return static_cast<const ExType&>(*this)/y;}

    /** @name Operations from the left
     *  For Python
     */
    ///@{
    inline ExType __radd__(const ExType& y) const
    { return y.__add__(static_cast<const ExType&>(*this));}
    inline ExType __rsub__(const ExType& y) const
    { return y.__sub__(static_cast<const ExType&>(*this));}
    inline ExType __rmul__(const ExType& y) const
    { return y.__mul__(static_cast<const ExType&>(*this));}
    inline ExType __rdiv__(const ExType& y) const
    { return y.__div__(static_cast<const ExType&>(*this));}
    inline ExType __rlt__(const ExType& y) const
    { return y.__lt__(static_cast<const ExType&>(*this));}
    inline ExType __rle__(const ExType& y) const
    { return y.__le__(static_cast<const ExType&>(*this));}
    inline ExType __rgt__(const ExType& y) const
    { return y.__gt__(static_cast<const ExType&>(*this));}
    inline ExType __rge__(const ExType& y) const
    { return y.__ge__(static_cast<const ExType&>(*this));}
    inline ExType __req__(const ExType& y) const
    { return y.__eq__(static_cast<const ExType&>(*this));}
    inline ExType __rne__(const ExType& y) const
    { return y.__ne__(static_cast<const ExType&>(*this));}
    inline ExType __rtruediv__(const ExType& y) const
    { return y.__truediv__(static_cast<const ExType&>(*this));}
    ///@}
    /// \endcond

};


} // namespace casadi

#endif // CASADI_GENERIC_EXPRESSION_HPP
