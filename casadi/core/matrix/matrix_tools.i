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
#ifndef CASADI_MATRIX_TOOLS_I
#define CASADI_MATRIX_TOOLS_I

%include <casadi/core/matrix/matrix.i>
%include <casadi/core/options_functionality.i>
%include <casadi/core/matrix/sparsity_tools.i>

%include <casadi/core/matrix/matrix_tools.hpp>

// map the template name to the instantiated name
#define MTT_INST(DataType, function_name)                       \
  %template(function_name) casadi::function_name <DataType >;

// Define template instantiations
#define MATRIX_TOOLS_TEMPLATES_COMMON(DataType)        \
  MTT_INST(DataType, transpose)                         \
  MTT_INST(DataType, mul)                               \
  MTT_INST(DataType, det)                               \
  MTT_INST(DataType, getMinor)                          \
  MTT_INST(DataType, cofactor)                          \
  MTT_INST(DataType, adj)                               \
  MTT_INST(DataType, inv)                               \
  MTT_INST(DataType, reshape)                           \
  MTT_INST(DataType, vec)                               \
  MTT_INST(DataType, vecNZ)                             \
  MTT_INST(DataType, blockcat)                          \
  MTT_INST(DataType, blocksplit)                        \
  MTT_INST(DataType, vertcat)                           \
  MTT_INST(DataType, vertsplit)                         \
  MTT_INST(DataType, horzcat)                           \
  MTT_INST(DataType, horzsplit)                         \
  MTT_INST(DataType, inner_prod)                        \
  MTT_INST(DataType, outer_prod)                        \
  MTT_INST(DataType, norm_1)                            \
  MTT_INST(DataType, norm_2)                            \
  MTT_INST(DataType, norm_inf)                          \
  MTT_INST(DataType, norm_F)                            \
  MTT_INST(DataType, norm_0_mul_nn)                     \
  MTT_INST(DataType, qr)                                \
  MTT_INST(DataType, nullspace)                         \
  MTT_INST(DataType, solve)                             \
  MTT_INST(DataType, pinv)                              \
  MTT_INST(DataType, repmat)                            \
  MTT_INST(DataType, unite)                             \
  MTT_INST(DataType, sumRows)                           \
  MTT_INST(DataType, sumCols)                           \
  MTT_INST(DataType, sumAll)                            \
  MTT_INST(DataType, trace)                             \
  MTT_INST(DataType, diag)                              \
  MTT_INST(DataType, blkdiag)                           \
  MTT_INST(DataType, polyval)                           \
  MTT_INST(DataType, addMultiple)                       \
  MTT_INST(DataType, veccat)                            \
  MTT_INST(DataType, vecNZcat)                          \
  MTT_INST(DataType, project)                           \
  MTT_INST(DataType, sprank)                            \
  MTT_INST(DataType, kron)

#ifdef SWIGOCTAVE
#define MATRIX_TOOLS_TEMPLATES(DataType) MATRIX_TOOLS_TEMPLATES_COMMON(DataType)
#else
#define MATRIX_TOOLS_TEMPLATES(DataType)               \
  MATRIX_TOOLS_TEMPLATES_COMMON(DataType)              \
  MTT_INST(DataType, sparse)                            \
  MTT_INST(DataType, dense)
#endif //SWIGOCTAVE

MATRIX_TOOLS_TEMPLATES(int)
MATRIX_TOOLS_TEMPLATES(double)
MATRIX_TOOLS_TEMPLATES(casadi::SXElement)

#endif // CASADI_MATRIX_TOOLS_I
