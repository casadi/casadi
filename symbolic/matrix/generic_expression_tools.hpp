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

#ifndef GENERIC_EXPRESSION_TOOLS_HPP
#define GENERIC_EXPRESSION_TOOLS_HPP

#include "../casadi_math.hpp"

namespace CasADi{

//@ {
/** \brief  Conditional operators */
template<class T>
T logic_and(const T &a, const T &b){ return a && b; }

template<class T>
T logic_or(const T &a, const T &b){ return a || b; }

template<class T>
T logic_not(const T &a){ return !a; }

} // namespace CasADi

#ifdef SWIG

// map the template name to the instantiated name
#define GET_INST(T,function_name) \
%template(function_name) CasADi::function_name < T >;

// Define template instanciations
#define GENERIC_EXPRESSION_TOOLS_TEMPLATES(T) \
GET_INST(T,logic_and) \
GET_INST(T,logic_or) \
GET_INST(T,logic_not) \

#endif //SWIG

#endif // GENERIC_EXPRESSION_TOOLS_HPP

