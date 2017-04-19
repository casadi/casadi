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


#ifndef CASADI_CASADI_RUNTIME_HPP
#define CASADI_CASADI_RUNTIME_HPP

#include <math.h>
#define CASADI_PREFIX(ID) casadi_##ID
#define CASADI_CAST(TYPE, ARG) static_cast<TYPE>(ARG)

/// \cond INTERNAL
namespace casadi {
  /// COPY: y <-x
  template<typename real_t>
  void CASADI_PREFIX(copy)(const real_t* x, int n, real_t* y);

  /// SWAP: x <-> y
  template<typename real_t>
  void CASADI_PREFIX(swap)(int n, real_t* x, int inc_x, real_t* y, int inc_y);

  /// Sparse copy: y <- x, w work vector (length >= number of rows)
  template<typename real_t>
  void CASADI_PREFIX(project)(const real_t* x, const int* sp_x, real_t* y, const int* sp_y, real_t* w);

  /// Convert sparse to dense
  template<typename real1_t, typename real2_t>
  void CASADI_PREFIX(densify)(const real1_t* x, const int* sp_x, real2_t* y, int tr);

  /// Convert dense to sparse
  template<typename real1_t, typename real2_t>
  void CASADI_PREFIX(sparsify)(const real1_t* x, real2_t* y, const int* sp_y, int tr);

  /// SCAL: x <- alpha*x
  template<typename real_t>
  void CASADI_PREFIX(scal)(int n, real_t alpha, real_t* x);

  /// AXPY: y <- a*x + y
  template<typename real_t>
  void CASADI_PREFIX(axpy)(int n, real_t alpha, const real_t* x, real_t* y);

  /// Inner product
  template<typename real_t>
  real_t CASADI_PREFIX(dot)(int n, const real_t* x, const real_t* y);

  /// Largest bound violation
  template<typename real_t>
  real_t CASADI_PREFIX(max_viol)(int n, const real_t* x, const real_t* lb, const real_t* ub);

  /// Sum of bound violations
  template<typename real_t>
  real_t CASADI_PREFIX(sum_viol)(int n, const real_t* x, const real_t* lb, const real_t* ub);

  /// IAMAX: index corresponding to the entry with the largest absolute value
  template<typename real_t>
  int CASADI_PREFIX(iamax)(int n, const real_t* x, int inc_x);

  /// FILL: x <- alpha
  template<typename real_t>
  void CASADI_PREFIX(fill)(real_t* x, int n, real_t alpha);

  /// Sparse matrix-matrix multiplication: z <- z + x*y
  template<typename real_t>
  void CASADI_PREFIX(mtimes)(const real_t* x, const int* sp_x, const real_t* y, const int* sp_y, real_t* z, const int* sp_z, real_t* w, int tr);

  /// Sparse matrix-vector multiplication: z <- z + x*y
  template<typename real_t>
  void CASADI_PREFIX(mv)(const real_t* x, const int* sp_x, const real_t* y, real_t* z, int tr);

  /// TRANS: y <- trans(x) , w work vector (length >= rows x)
  template<typename real_t>
  void CASADI_PREFIX(trans)(const real_t* x, const int* sp_x, real_t* y, const int* sp_y, int *tmp);

  /// NORM_1: ||x||_1 -> return
  template<typename real_t>
  real_t CASADI_PREFIX(norm_1)(int n, const real_t* x);

  /// NORM_2: ||x||_2 -> return
  template<typename real_t>
  real_t CASADI_PREFIX(norm_2)(int n, const real_t* x);

  /** Inf-norm of a vector *
      Returns the largest element in absolute value
   */
  template<typename real_t>
  real_t CASADI_PREFIX(norm_inf)(int n, const real_t* x);

  /** Inf-norm of a Matrix-matrix product,*
   * \param dwork  A real work vector that you must allocate
   *               Minimum size: y.size1()
   * \param iwork  A integer work vector that you must allocate
   *               Minimum size: y.size1()+x.size2()+1
   */
  template<typename real_t>
  real_t CASADI_PREFIX(norm_inf_mul)(const real_t* x, const int* sp_x, const real_t* y, const int* sp_y,
                             real_t *dwork, int *iwork);

  /** Calculates dot(x, mul(A, y)) */
  template<typename real_t>
  real_t CASADI_PREFIX(bilin)(const real_t* A, const int* sp_A, const real_t* x, const real_t* y);

  /// Adds a multiple alpha/2 of the outer product mul(x, trans(x)) to A
  template<typename real_t>
  void CASADI_PREFIX(rank1)(real_t* A, const int* sp_A, real_t alpha, const real_t* x);

  /// Get the nonzeros for the upper triangular half
  template<typename real_t>
  void CASADI_PREFIX(getu)(const real_t* x, const int* sp_x, real_t* v);

  /// Evaluate a polynomial
  template<typename real_t>
  real_t CASADI_PREFIX(polyval)(const real_t* p, int n, real_t x);

  // Loop over corners of a hypercube
  int CASADI_PREFIX(flip)(int* corner, int ndim);

  // Find the interval to which a value belongs
  template<typename real_t>
  int CASADI_PREFIX(low)(real_t x, const double* grid, int ng, int lookup_mode);

  // Get weights for the multilinear interpolant
  template<typename real_t>
  void CASADI_PREFIX(interpn_weights)(int ndim, const real_t* grid, const int* offset, const real_t* x, real_t* alpha, int* index);

  // Get coefficients for the multilinear interpolant
  template<typename real_t>
  real_t CASADI_PREFIX(interpn_interpolate)(int ndim, const int* offset, const real_t* values, const real_t* alpha, const int* index, const int* corner, real_t* coeff);

  // Multilinear interpolant
  template<typename real_t>
  real_t CASADI_PREFIX(interpn)(int ndim, const real_t* grid, const int* offset, const real_t* values, const real_t* x, int* iw, real_t* w);

  // Multilinear interpolant - calculate gradient
  template<typename real_t>
  void CASADI_PREFIX(interpn_grad)(real_t* grad, int ndim, const real_t* grid, const int* offset, const real_t* values, const real_t* x, int* iw, real_t* w);

  // De boor single basis evaluation
  template<typename real_t>
  void CASADI_PREFIX(de_boor)(real_t x, const real_t* knots, int n_knots, int degree, real_t* boor);

  // De boor nd evaluation
  template<typename real_t>
  void CASADI_PREFIX(nd_boor_eval)(real_t* ret, int n_dims, const real_t* knots, const int* offset, const int* degree, const int* strides, const real_t* c, int m, const real_t* x, const int* lookup_mode, int reverse, int* iw, real_t* w);

}

// Implementations

// Note: due to restrictions of cmake IO processing, make sure that
//      1)   semicolons (;) are never immediately preceded by whitespace
//      2)   line continuation slashes (\) are always immediately preceded by whitespace
#define CASADI_GEMM_NT(M, N, K, A, LDA, B, LDB, C, LDC) \
  for (i=0, rr=C; i<M; ++i) \
    for (j=0; j<N; ++j, ++rr) \
      for (k=0, ss=A+i*LDA, tt=B+j*LDB; k<K; ++k) \
        *rr += *ss++**tt++;


namespace casadi {
  #include "casadi_copy.hpp"
  #include "casadi_swap.hpp"
  #include "casadi_project.hpp"
  #include "casadi_densify.hpp"
  #include "casadi_sparsify.hpp"
  #include "casadi_scal.hpp"
  #include "casadi_iamax.hpp"
  #include "casadi_axpy.hpp"
  #include "casadi_dot.hpp"
  #include "casadi_fill.hpp"
  #include "casadi_max_viol.hpp"
  #include "casadi_sum_viol.hpp"
  #include "casadi_mtimes.hpp"
  #include "casadi_mv.hpp"
  #include "casadi_trans.hpp"
  #include "casadi_norm_1.hpp"
  #include "casadi_norm_2.hpp"
  #include "casadi_norm_inf.hpp"
  #include "casadi_norm_inf_mul.hpp"
  #include "casadi_bilin.hpp"
  #include "casadi_rank1.hpp"
  #include "casadi_low.hpp"
  #include "casadi_flip.hpp"
  #include "casadi_polyval.hpp"
  #include "casadi_de_boor.hpp"
  #include "casadi_nd_boor_eval.hpp"
  #include "casadi_interpn_weights.hpp"
  #include "casadi_interpn_interpolate.hpp"
  #include "casadi_interpn.hpp"
  #include "casadi_interpn_grad.hpp"

} // namespace casadi



/// \endcond

#endif // CASADI_CASADI_RUNTIME_HPP
