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

#ifndef CASADI_RUNTIME_H
#define CASADI_RUNTIME_H

#include <math.h>
#ifdef __cplusplus
extern "C" {
#endif

  typedef double real_t;

  /// COPY: y <-x
  void casadi_copy(int n, const real_t* x, int inc_x, real_t* y, int inc_y);
  
  /// SWAP: x <-> y
  void casadi_swap(int n, real_t* x, int inc_x, real_t* y, int inc_y);
  
  /// COPY sparse: y <- x
  void casadi_copy_sparse(const real_t* x, const int* sp_x, real_t* y, const int* sp_y);
  
  /// SCAL: x <- alpha*x
  void casadi_scal(int n, real_t alpha, real_t* x, int inc_x);
  
  /// AXPY: y <- a*x + y
  void casadi_axpy(int n, real_t alpha, const real_t* x, int inc_x, real_t* y, int inc_y);
  
  /// DOT: inner_prod(x,y) -> return
  real_t casadi_dot(int n, const real_t* x, int inc_x, real_t* y, int inc_y);
  
  /// ASUM: ||x||_1 -> return
  real_t casadi_asum(int n, const real_t* x, int inc_x);
  
  /// IAMAX: index corresponding to the entry with the largest absolute value 
  int casadi_iamax(int n, const real_t* x, int inc_x);
  
  /// FILL: x <- alpha
  void casadi_fill(int n, real_t alpha, real_t* x, int inc_x);
  
  /// Sparse matrix-matrix multiplication, the second argument is transposed: z <- z + x*y'
  void casadi_mm_nt_sparse(const real_t* x, const int* sp_x, const real_t* trans_y, const int* sp_trans_y, real_t* z, const int* sp_z);
  
  /// NRM2: ||x||_2 -> return
  real_t casadi_nrm2(int n, const real_t* x, int inc_x);
  
  /// TRANS: y <- trans(x)
  void casadi_trans(const real_t* x, const int* sp_x, real_t* y, const int* sp_y, int *tmp);

#ifdef __cplusplus
} // extern "C"
#endif

#endif // CASADI_RUNTIME_H
