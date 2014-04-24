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

namespace casadi {

  /** \brief  Logical `and`, returns (an expression evaluating to) 1 if both
   *          expressions are nonzero and 0 otherwise */
  template<typename DataType>
  DataType logic_and(const DataType& x, const DataType& y) { return x && y; }

  /** \brief  Logical `or`, returns (an expression evaluating to) 1 if at
   *          least one expression is nonzero and 0 otherwise */
  template<typename DataType>
  DataType logic_or(const DataType& x, const DataType& y) { return x || y; }

  /** \brief  Logical `not`, returns (an expression evaluating to) 1 if
   *          expression is zero and 0 otherwise */
  template<typename DataType>
  DataType logic_not(const DataType &x) { return !x; }

} // namespace casadi

#ifdef SWIG

// map the template name to the instantiated name
#define GET_INST(DataType, function_name) \
%template(function_name) casadi::function_name< DataType >;

// Define template instantiations
#define GENERIC_EXPRESSION_TOOLS_TEMPLATES(DataType) \
GET_INST(DataType, logic_and) \
GET_INST(DataType, logic_or) \
GET_INST(DataType, logic_not) \

#endif //SWIG

#endif // GENERIC_EXPRESSION_TOOLS_HPP
