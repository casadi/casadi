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
#include <limits>

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
    static bool isAlmostZero(const T& val, double tol){ return val<=tol && val>=-tol; }
    static bool isOne(const T& val){ return val==1;}
    static bool isMinusOne(const T& val){ return val==-1;}
    static bool isConstant(const T& val){ return true;}
    static bool isInteger(const T& val){ return val==int(val);}
    static bool isInf(const T& val){ return std::numeric_limits<T>::has_infinity ? val==std::numeric_limits<T>::infinity() : false;}
    static bool isMinusInf(const T& val){ return std::numeric_limits<T>::has_infinity ? val==-std::numeric_limits<T>::infinity() : false;}
    static bool isNaN(const T& val){ return std::numeric_limits<T>::has_quiet_NaN ? val!=val : false;}
    static const T zero;
    static const T one;
    static const T two;
    static const T minus_one;
};

template<class T>
const T casadi_limits<T>::zero = T(0);

template<class T>
const T casadi_limits<T>::one = 1;

template<class T>
const T casadi_limits<T>::two = 2;

template<class T>
const T casadi_limits<T>::minus_one = -1;

} // namespace CasADi
#endif // CASADI_LIMITS_HPP

