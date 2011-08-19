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

#ifndef CASADI_OPERATORS_HPP
#define CASADI_OPERATORS_HPP

#include <cmath>

/** \brief casadi_operators class
\author Joel Andersson
\date 2011
*/


namespace CasADi{

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
    static T constpow(const T&x, const T&y){ return std::pow(x,y);}
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
    static T erf(const T&x){ return ::erf(x);}
    static T sinh(const T&x){ return std::sinh(x);}
    static T cosh(const T&x){ return std::cosh(x);}
    static T tanh(const T&x){ return std::tanh(x);}
};

//@{
/** \brief  Forward declarations */
  class MX;
  //@}
} // namespace CasADi

#define MX CasADi::MX
namespace std{
//@{
/** \brief  Pre-C99 elementary functions from the math.h (cmath) header */
MX sqrt(const MX &x);
MX sin(const MX &x);
MX cos(const MX &x);
MX tan(const MX &x);
MX atan(const MX &x);
MX asin(const MX &x);
MX acos(const MX &x);
MX exp(const MX &x);
MX log(const MX &x);
MX log10(const MX &x);
MX constpow(const MX &x, const MX &n);
MX pow(const MX &x, const MX &n);
MX abs(const MX &x);
MX fabs(const MX &x); // same as abs
MX floor(const MX &x);
MX ceil(const MX &x);
MX sinh(const MX &x);
MX cosh(const MX &x);
MX tanh(const MX &x);
//@}
} // namespace std

/** \brief  C99 elementary functions from the math.h header */
MX erf(const MX &x);
MX fmin(const MX &a, const MX &b);
MX fmax(const MX &a, const MX &b);
#undef MX

#endif // CASADI_OPERATORS_HPP

