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
#ifndef CASADI_GENERIC_MATRIX_TOOLS_I
#define CASADI_GENERIC_MATRIX_TOOLS_I

%include <casadi/core/matrix/slice.i>
%include <casadi/core/matrix/sparsity.i>

%include <casadi/core/matrix/generic_matrix_tools.hpp>

// map the template name to the instantiated name
#define GMTT_INST(MatType, function_name) \
%template(function_name) casadi::function_name< MatType >;

// Define template instantiations
#define GENERIC_MATRIX_TOOLS_TEMPLATES(MatType) \
GMTT_INST(MatType, cross) \
GMTT_INST(MatType, quad_form) \
GMTT_INST(MatType, sum_square) \
GMTT_INST(MatType, tril2symm) \
GMTT_INST(MatType, triu2symm) \
GMTT_INST(MatType, triu) \
GMTT_INST(MatType, tril) \
GMTT_INST(MatType, isEqual) \
GMTT_INST(MatType, diagsplit)

// Define template instantiations
#define GENERIC_MATRIX_TOOLS_TEMPLATES_REAL_ONLY(MatType) \
GMTT_INST(MatType, linspace)

GENERIC_MATRIX_TOOLS_TEMPLATES(casadi::Matrix<int>)
GENERIC_MATRIX_TOOLS_TEMPLATES(casadi::Matrix<double>)
GENERIC_MATRIX_TOOLS_TEMPLATES(casadi::Matrix<casadi::SXElement>)
GENERIC_MATRIX_TOOLS_TEMPLATES(casadi::MX)

GENERIC_MATRIX_TOOLS_TEMPLATES_REAL_ONLY(casadi::Matrix<double>)
GENERIC_MATRIX_TOOLS_TEMPLATES_REAL_ONLY(casadi::Matrix<casadi::SXElement>)
GENERIC_MATRIX_TOOLS_TEMPLATES_REAL_ONLY(casadi::MX)

#endif // CASADI_GENERIC_MATRIX_TOOLS_I
