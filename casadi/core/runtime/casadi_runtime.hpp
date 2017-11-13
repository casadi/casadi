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

#include "../calculus.hpp"

#define CASADI_PREFIX(ID) casadi_##ID
#define CASADI_CAST(TYPE, ARG) static_cast<TYPE>(ARG)

/// \cond INTERNAL
namespace casadi {
  /// COPY: y <-x
  template<typename T1>
  void casadi_copy(const T1* x, int n, T1* y);

  /// SWAP: x <-> y
  template<typename T1>
  void casadi_swap(int n, T1* x, int inc_x, T1* y, int inc_y);

  /// Sparse copy: y <- x, w work vector (length >= number of rows)
  template<typename T1>
  void casadi_project(const T1* x, const int* sp_x, T1* y, const int* sp_y, T1* w);

  /// Convert sparse to dense
  template<typename T1, typename T2>
  void casadi_densify(const T1* x, const int* sp_x, T2* y, int tr);

  /// Convert dense to sparse
  template<typename T1, typename T2>
  void casadi_sparsify(const T1* x, T2* y, const int* sp_y, int tr);

  /// SCAL: x <- alpha*x
  template<typename T1>
  void casadi_scal(int n, T1 alpha, T1* x);

  /// AXPY: y <- a*x + y
  template<typename T1>
  void casadi_axpy(int n, T1 alpha, const T1* x, T1* y);

  /// Inner product
  template<typename T1>
  T1 casadi_dot(int n, const T1* x, const T1* y);

  /// Largest bound violation
  template<typename T1>
  T1 casadi_max_viol(int n, const T1* x, const T1* lb, const T1* ub);

  /// Sum of bound violations
  template<typename T1>
  T1 casadi_sum_viol(int n, const T1* x, const T1* lb, const T1* ub);

  /// IAMAX: index corresponding to the entry with the largest absolute value
  template<typename T1>
  int casadi_iamax(int n, const T1* x, int inc_x);

  /// FILL: x <- alpha
  template<typename T1>
  void casadi_fill(T1* x, int n, T1 alpha);

  /// Sparse matrix-matrix multiplication: z <- z + x*y
  template<typename T1>
  void casadi_mtimes(const T1* x, const int* sp_x, const T1* y, const int* sp_y,
                             T1* z, const int* sp_z, T1* w, int tr);

  /// Sparse matrix-vector multiplication: z <- z + x*y
  template<typename T1>
  void casadi_mv(const T1* x, const int* sp_x, const T1* y, T1* z, int tr);

  /// TRANS: y <- trans(x) , w work vector (length >= rows x)
  template<typename T1>
  void casadi_trans(const T1* x, const int* sp_x, T1* y, const int* sp_y, int* tmp);

  /// NORM_1: ||x||_1 -> return
  template<typename T1>
  T1 casadi_norm_1(int n, const T1* x);

  /// NORM_2: ||x||_2 -> return
  template<typename T1>
  T1 casadi_norm_2(int n, const T1* x);

  /** Inf-norm of a vector *
      Returns the largest element in absolute value
   */
  template<typename T1>
  T1 casadi_norm_inf(int n, const T1* x);

  /** Inf-norm of a Matrix-matrix product,*
   * \param dwork  A real work vector that you must allocate
   *               Minimum size: y.size1()
   * \param iwork  A integer work vector that you must allocate
   *               Minimum size: y.size1()+x.size2()+1
   */
  template<typename T1>
  T1 casadi_norm_inf_mul(const T1* x, const int* sp_x, const T1* y, const int* sp_y,
                             T1* dwork, int* iwork);

  /** Calculates dot(x, mul(A, y)) */
  template<typename T1>
  T1 casadi_bilin(const T1* A, const int* sp_A, const T1* x, const T1* y);

  /// Adds a multiple alpha/2 of the outer product mul(x, trans(x)) to A
  template<typename T1>
  void casadi_rank1(T1* A, const int* sp_A, T1 alpha, const T1* x);

  /// Get the nonzeros for the upper triangular half
  template<typename T1>
  void casadi_getu(const T1* x, const int* sp_x, T1* v);

  /// Evaluate a polynomial
  template<typename T1>
  T1 casadi_polyval(const T1* p, int n, T1 x);

  // Loop over corners of a hypercube
  int casadi_flip(int* corner, int ndim);

  // Find the interval to which a value belongs
  template<typename T1>
  int casadi_low(T1 x, const double* grid, int ng, int lookup_mode);

  // Get weights for the multilinear interpolant
  template<typename T1>
  void casadi_interpn_weights(int ndim, const T1* grid, const int* offset,
                                      const T1* x, T1* alpha, int* index);

  // Get coefficients for the multilinear interpolant
  template<typename T1>
  T1 casadi_interpn_interpolate(int ndim, const int* offset, const T1* values,
                                        const T1* alpha, const int* index,
                                        const int* corner, T1* coeff);

  // Multilinear interpolant
  template<typename T1>
  T1 casadi_interpn(int ndim, const T1* grid, const int* offset, const T1* values,
                            const T1* x, int* iw, T1* w);

  // Multilinear interpolant - calculate gradient
  template<typename T1>
  void casadi_interpn_grad(T1* grad, int ndim, const T1* grid, const int* offset,
                                   const T1* values, const T1* x, int* iw, T1* w);

  // De boor single basis evaluation
  template<typename T1>
  void casadi_de_boor(T1 x, const T1* knots, int n_knots, int degree, T1* boor);

  // De boor nd evaluation
  template<typename T1>
  void casadi_nd_boor_eval(T1* ret, int n_dims, const T1* knots, const int* offset,
                                   const int* degree, const int* strides, const T1* c, int m,
                                   const T1* x, const int* lookup_mode, int reverse, int* iw,
                                   T1* w);

  // Alias names
  inline void casadi_copy_int(const int* x, int n, int* y) {
    casadi_copy(x, n, y);
  }
  inline void casadi_fill_int(int* x, int n, int alpha) {
    casadi_fill(x, n, alpha);
  }

  // Dense matrix multiplication
  #define CASADI_GEMM_NT(M, N, K, A, LDA, B, LDB, C, LDC) \
    for (i=0, rr=C; i<M; ++i) \
      for (j=0; j<N; ++j, ++rr) \
        for (k=0, ss=A+i*LDA, tt=B+j*LDB; k<K; ++k) \
          *rr += *ss++**tt++;

  // Template function implementations
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
  #include "casadi_mv_dense.hpp"
  #include "casadi_finite_diff.hpp"
  #include "casadi_ldl.hpp"
  #include "casadi_qr.hpp"
} // namespace casadi

/// \endcond

#endif // CASADI_CASADI_RUNTIME_HPP
