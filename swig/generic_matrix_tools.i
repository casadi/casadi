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

%include "slice.i"
%include "sparsity.i"

%include <casadi/core/matrix/generic_matrix_tools.hpp>

%define SPARSITY_INTERFACE_DECL(MatType...)
MatType horzcat(const std::vector< MatType > &v);
std::vector<MatType > horzsplit(const MatType &v, const std::vector<int>& offset);
std::vector<MatType > horzsplit(const MatType &v, int incr=1);
MatType vertcat(const std::vector< MatType > &v);
std::vector<MatType > vertsplit(const MatType &v, const std::vector<int>& offset);
std::vector<MatType > vertsplit(const MatType &v, int incr=1);
MatType blockcat(const std::vector< std::vector<MatType > > &v);
std::vector< std::vector< MatType > >
blocksplit(const MatType& x, const std::vector<int>& vert_offset, const std::vector<int>& horz_offset);
std::vector< std::vector< MatType > > blocksplit(const MatType& x, int vert_incr=1, int horz_incr=1);
MatType diagcat(const std::vector< MatType > &v);
std::vector< MatType > diagsplit(const MatType& x, const std::vector<int>& output_offset1,
                                 const std::vector<int>& output_offset2);
std::vector< MatType > diagsplit(const MatType& x, const std::vector<int>& output_offset);
std::vector< MatType > diagsplit(const MatType& x, int incr=1);
std::vector< MatType > diagsplit(const MatType& x, int incr1, int incr2);
MatType veccat(const std::vector< MatType >& x);
MatType vecNZcat(const std::vector< MatType >& x);
MatType mul(const MatType &x, const MatType &y);
MatType mul(const MatType &x, const MatType &y, const MatType &z);
MatType mul(const std::vector< MatType > &args);
MatType transpose(const MatType &X);
MatType vec(const MatType& a);
MatType vecNZ(const MatType& a);
MatType reshape(const MatType& a, const Sparsity& sp);
MatType reshape(const MatType& a, int nrow, int ncol);
MatType reshape(const MatType& a, std::pair<int, int> rc);
int sprank(const MatType& A);
MatType triu(const MatType& a, bool includeDiagonal=true);
MatType tril(const MatType& a, bool includeDiagonal=true);
%enddef

#ifndef SWIGMATLAB
SPARSITY_INTERFACE_DECL(casadi::Sparsity)
#endif // SWIGMATLAB

%define GENERIC_MATRIX_DECL(MatType...)
MatType quad_form(const MatType &X, const MatType &A);
MatType quad_form(const MatType &X);
MatType sum_square(const MatType &X);
MatType linspace(const MatType &a, const MatType &b, int nsteps);
MatType cross(const MatType &a, const MatType &b, int dim = -1);
MatType det(const MatType& A);
MatType inv(const MatType& A);
MatType trace(const MatType& a);
%enddef

%define MATRIX_DECL(MatType...)
MatType adj(const MatType& A);
MatType getMinor(const MatType &x, int i, int j);
MatType cofactor(const MatType &x, int i, int j);
%enddef

// map the template name to the instantiated name
%define GMTT_INST(MatType, function_name...)
%template(function_name) casadi::function_name< MatType >;
%enddef

// Define template instantiations
%define GENERIC_MATRIX_TOOLS_TEMPLATES(MatType...)
SPARSITY_INTERFACE_DECL(MatType)
GENERIC_MATRIX_DECL(MatType)
GMTT_INST(MatType, tril2symm)
GMTT_INST(MatType, triu2symm)
GMTT_INST(MatType, isEqual)
%enddef

%define GENERIC_MATRIX_TOOLS_TEMPLATES_MATRIX(DataType...)
GENERIC_MATRIX_TOOLS_TEMPLATES(casadi::Matrix<DataType>)
MATRIX_DECL(casadi::Matrix<DataType>)
%enddef

#ifndef SWIGMATLAB
GENERIC_MATRIX_TOOLS_TEMPLATES_MATRIX(int)
GENERIC_MATRIX_TOOLS_TEMPLATES_MATRIX(double)
GENERIC_MATRIX_TOOLS_TEMPLATES_MATRIX(casadi::SXElement)
GENERIC_MATRIX_TOOLS_TEMPLATES(casadi::MX)
#endif // SWIGMATLAB

#endif // CASADI_GENERIC_MATRIX_TOOLS_I
