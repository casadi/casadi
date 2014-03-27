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

#ifndef CASADI_HPP
#define CASADI_HPP

// Scalar expressions (why do i need to put it up here?)
#include "sx/sx_element.hpp"

// Generic tools
#include "polynomial.hpp"
#include "matrix/generic_matrix_tools.hpp"
#include "matrix/generic_expression_tools.hpp"

// Matrices
#include "matrix/matrix.hpp"
#include "matrix/matrix_tools.hpp"
#include "matrix/sparsity_tools.hpp"

// Scalar expressions
#include "sx/sx_tools.hpp"

// Matrix expressions
#include "mx/mx.hpp"
#include "mx/mx_tools.hpp"

// Misc functions
#include "function/sx_function.hpp"
#include "function/mx_function.hpp"
#include "function/external_function.hpp"

#endif //CASADI_HPP
