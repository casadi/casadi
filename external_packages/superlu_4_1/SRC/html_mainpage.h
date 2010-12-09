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

/*! \mainpage SuperLU Documentation
 
SuperLU is a general purpose library for the direct solution of large,
sparse, nonsymmetric systems of linear equations on high performance
machines. The library is written in C and is callable from either C or
Fortran. The library routines perform an LU decomposition with
partial pivoting and triangular system solves through forward and back
substitution.  The library also provides threshold-based ILU factorization
preconditioners.

The factorization routines can handle non-square
matrices but the triangular solves are performed only for square
matrices. The matrix columns may be preordered (before factorization)
either through library or user supplied routines. This preordering for
sparsity is completely separate from the factorization. Working
precision iterative refinement subroutines are provided for improved
backward stability. Routines are also provided to equilibrate the
system, estimate the condition number, calculate the relative backward
error, and estimate error bounds for the refined solutions. 
 
 */
