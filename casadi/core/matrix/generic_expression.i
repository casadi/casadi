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
#ifndef CASADI_GENERIC_EXPRESSION_I
#define CASADI_GENERIC_EXPRESSION_I

%include <casadi/core/printable_object.i>

%include <casadi/core/matrix/generic_expression.hpp>

%template(ExpIMatrix)        casadi::GenericExpression<casadi::Matrix<int> >;
%template(ExpDMatrix)        casadi::GenericExpression<casadi::Matrix<double> >;
%template(ExpSX)             casadi::GenericExpression<casadi::Matrix<casadi::SXElement> >;
%template(ExpMX)             casadi::GenericExpression<casadi::MX>;
%template(ExpSXElement)      casadi::GenericExpression<casadi::SXElement>;

#endif // CASADI_GENERIC_EXPRESSION_I

