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
#include "generic_expression_friends.hpp"

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
