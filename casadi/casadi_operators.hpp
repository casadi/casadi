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

#include "pre_c99_support.hpp"
#include "casadi_types.hpp"

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
    static T pow(const T&x, const T&y){ return ::pow(x,y);}
    static T equality(const T&x, const T&y){ return x==y;}
    static T fmin(const T&x, const T&y){ return std::min(x,y);}
    static T fmax(const T&x, const T&y){ return std::max(x,y);}
    static T constpow(const T&x, const T&y){ return ::pow(x,y);}
    static T printme(const T&x, const T&y){ return printme(x,y);} // BUG? Infinite loop?
};
} // namespace CasADi

#endif // CASADI_OPERATORS_HPP

