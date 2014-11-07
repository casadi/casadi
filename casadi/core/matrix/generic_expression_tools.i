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
#ifndef CASADI_GENERIC_EXPRESSION_TOOLS_I
#define CASADI_GENERIC_EXPRESSION_TOOLS_I

%include <casadi/core/matrix/generic_expression_tools.hpp>

// map the template name to the instantiated name
#define GET_INST(DataType, function_name) \
%template(function_name) casadi::function_name< DataType >;

// Define template instantiations
#define GENERIC_EXPRESSION_TOOLS_TEMPLATES(DataType) \
GET_INST(DataType, logic_and) \
GET_INST(DataType, logic_or) \
GET_INST(DataType, logic_not) \

GENERIC_EXPRESSION_TOOLS_TEMPLATES(int)
GENERIC_EXPRESSION_TOOLS_TEMPLATES(double)
GENERIC_EXPRESSION_TOOLS_TEMPLATES(casadi::Matrix<int>)
GENERIC_EXPRESSION_TOOLS_TEMPLATES(casadi::Matrix<double>)
GENERIC_EXPRESSION_TOOLS_TEMPLATES(casadi::Matrix<casadi::SXElement>)
GENERIC_EXPRESSION_TOOLS_TEMPLATES(casadi::MX)
GENERIC_EXPRESSION_TOOLS_TEMPLATES(casadi::SXElement)

#endif // CASADI_GENERIC_EXPRESSION_TOOLS_I
