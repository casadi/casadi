/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010 by Joel Andersson, Moritz Diehl, K.U.Leuven. All rights reserved.
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

#ifndef CASADI_LIMITS_HPP
#define CASADI_LIMITS_HPP

#include <cmath>

/** \brief casadi_limits class
The following class, which acts as a complements to the standard numeric_limits class, allows
to specify certain properties of scalar objects. The template can be specialized for 
e.g. symbolic scalars
\author Joel Andersson
\date 2011
*/


namespace CasADi{
  
template<class T>
class casadi_limits{
  public:
    static bool isZero(const T& val){ return val==0; }
    static bool isOne(const T& val){ return val==1;}
    static bool isConstant(const T& val){ return true;}
    static bool isInteger(const T& val){ return val==int(val);}
    static const T zero = 0;
    static const T one = 1;
    static const T two = 2;
    static const T minus_one = -1;
};

template<class T>
class casadi_operators{
  public:
    static T add(const T&x, const T&y){ return x+y;}
    static T sub(const T&x, const T&y){ return x-y;}
    static T mul(const T&x, const T&y){ return x*y;}
    static T div(const T&x, const T&y){ return x/y;}
    static T neg(const T&x){ return -x;}
    static T exp(const T&x){ return std::exp(x);}
    static T log(const T&x){ return std::log(x);}
    static T pow(const T&x, const T&y){ return std::pow(x,y);}
    static T sqrt(const T&x){ return std::sqrt(x);}
    static T sin(const T&x){ return std::sin(x);}
    static T cos(const T&x){ return std::cos(x);}
    static T tan(const T&x){ return std::tan(x);}
    static T asin(const T&x){ return std::asin(x);}
    static T acos(const T&x){ return std::acos(x);}
    static T atan(const T&x){ return std::atan(x);}
    static T floor(const T&x){ return std::floor(x);}
    static T ceil(const T&x){ return std::ceil(x);}
    static T equality(const T&x, const T&y){ return x==y;}
    static T fmin(const T&x, const T&y){ return std::min(x,y);}
    static T fmax(const T&x, const T&y){ return std::max(x,y);}
    static T fabs(const T&x){ return std::abs(x);}
};


} // namespace CasADi
#endif // CASADI_LIMITS_HPP

