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

  template<typename real_t>
  void CASADI_PREFIX(copy)(const real_t* x, int n, real_t* y) {
    int i;
    if (y) {
      if (x) {
        for (i=0; i<n; ++i) *y++ = *x++;
      } else {
        for (i=0; i<n; ++i) *y++ = 0.;
      }
    }
  }

  template<typename real_t>
  void CASADI_PREFIX(swap)(int n, real_t* x, int inc_x, real_t* y, int inc_y) {
    real_t t;
    int i;
    for (i=0; i<n; ++i) {
      t = *x;
      *x = *y;
      *y = t;
      x += inc_x;
      y += inc_y;
    }
  }

  template<typename real_t>
  void CASADI_PREFIX(project)(const real_t* x, const int* sp_x, real_t* y, const int* sp_y, real_t* w) {
    int ncol_x = sp_x[1];
    const int *colind_x = sp_x+2, *row_x = sp_x + 2 + ncol_x+1;
    int ncol_y = sp_y[1];
    const int *colind_y = sp_y+2, *row_y = sp_y + 2 + ncol_y+1;
    /* Loop over columns of x and y */
    int i, el;
    for (i=0; i<ncol_x; ++i) {
      /* Zero out requested entries in y */
      for (el=colind_y[i]; el<colind_y[i+1]; ++el) w[row_y[el]] = 0;
      /* Set x entries */
      for (el=colind_x[i]; el<colind_x[i+1]; ++el) w[row_x[el]] = x[el];
      /* Retrieve requested entries in y */
      for (el=colind_y[i]; el<colind_y[i+1]; ++el) y[el] = w[row_y[el]];
    }
  }

  template<typename real1_t, typename real2_t>
  void CASADI_PREFIX(densify)(const real1_t* x, const int* sp_x, real2_t* y, int tr) {
    /* Quick return - output ignored */
    if (!y) return;
    int nrow_x = sp_x[0], ncol_x = sp_x[1];
    const int *colind_x = sp_x+2, *row_x = sp_x+ncol_x+3;
    /* Zero out return value */
    CASADI_PREFIX(fill)(y, nrow_x*ncol_x, CASADI_CAST(real2_t, 0));
    /* Quick return - input is zero */
    if (!x) return;
    /* Copy nonzeros */
    int i, el;
    if (tr) {
      for (i=0; i<ncol_x; ++i) {
        for (el=colind_x[i]; el!=colind_x[i+1]; ++el) {
          y[i + row_x[el]*ncol_x] = CASADI_CAST(real2_t, *x++);
        }
      }
    } else {
      for (i=0; i<ncol_x; ++i) {
        for (el=colind_x[i]; el!=colind_x[i+1]; ++el) {
          y[row_x[el] + i*nrow_x] = CASADI_CAST(real2_t, *x++);
        }
      }
    }
  }

  template<typename real1_t, typename real2_t>
  void CASADI_PREFIX(sparsify)(const real1_t* x, real2_t* y, const int* sp_y, int tr) {
    int nrow_y = sp_y[0], ncol_y = sp_y[1];
    const int *colind_y = sp_y+2, *row_y = sp_y+ncol_y+3;
    int i, el;
    if (tr) {
      for (i=0; i<ncol_y; ++i) {
        for (el=colind_y[i]; el!=colind_y[i+1]; ++el) {
          *y++ = CASADI_CAST(real2_t, x[i + row_y[el]*ncol_y]);
        }
      }
    } else {
      for (i=0; i<ncol_y; ++i) {
        for (el=colind_y[i]; el!=colind_y[i+1]; ++el) {
          *y++ = CASADI_CAST(real2_t, x[row_y[el] + i*nrow_y]);
        }
      }
    }
  }

  template<typename real_t>
  void CASADI_PREFIX(scal)(int n, real_t alpha, real_t* x) {
    int i;
    for (i=0; i<n; ++i) *x++ *= alpha;
  }

  template<typename real_t>
  void CASADI_PREFIX(axpy)(int n, real_t alpha, const real_t* x, real_t* y) {
    int i;
    for (i=0; i<n; ++i) *y++ += alpha**x++;
  }

  template<typename real_t>
  real_t CASADI_PREFIX(dot)(int n, const real_t* x, const real_t* y) {
    real_t r = 0;
    int i;
    for (i=0; i<n; ++i) r += *x++ * *y++;
    return r;
  }

  template<typename real_t>
  real_t CASADI_PREFIX(max_viol)(int n, const real_t* x, const real_t* lb, const real_t* ub) {
    real_t r = 0;
    const real_t zero = 0;
    int i;
    for (i=0; i<n; ++i) {
      real_t x_i = x ? *x++ : zero;
      real_t lb_i = lb ? *lb++ : zero;
      real_t ub_i = ub ? *ub++ : zero;
      r = fmax(r, fmax(x_i-ub_i, zero));
      r = fmax(r, fmax(lb_i-x_i, zero));
    }
    return r;
  }

  template<typename real_t>
  real_t CASADI_PREFIX(sum_viol)(int n, const real_t* x, const real_t* lb, const real_t* ub) {
    real_t r = 0;
    const real_t zero = 0;
    int i;
    for (i=0; i<n; ++i) {
      real_t x_i = x ? *x++ : zero;
      real_t lb_i = lb ? *lb++ : zero;
      real_t ub_i = ub ? *ub++ : zero;
      r += fmax(x_i-ub_i, zero);
      r += fmax(lb_i-x_i, zero);
    }
    return r;
  }

  template<typename real_t>
  int CASADI_PREFIX(iamax)(int n, const real_t* x, int inc_x) {
    real_t t;
    real_t largest_value = -1.0;
    int largest_index = -1;
    int i;
    for (i=0; i<n; ++i) {
      t = fabs(*x);
      x += inc_x;
      if (t>largest_value) {
        largest_value = t;
        largest_index = i;
      }
    }
    return largest_index;
  }

  template<typename real_t>
  void CASADI_PREFIX(fill)(real_t* x, int n, real_t alpha) {
    int i;
    if (x) {
      for (i=0; i<n; ++i) *x++ = alpha;
    }
  }

  template<typename real_t>
  void CASADI_PREFIX(mtimes)(const real_t* x, const int* sp_x, const real_t* y, const int* sp_y, real_t* z, const int* sp_z, real_t* w, int tr) {
    /* Get sparsities */
    int ncol_x = sp_x[1];
    const int *colind_x = sp_x+2, *row_x = sp_x + 2 + ncol_x+1;
    int ncol_y = sp_y[1];
    const int *colind_y = sp_y+2, *row_y = sp_y + 2 + ncol_y+1;
    int ncol_z = sp_z[1];
    const int *colind_z = sp_z+2, *row_z = sp_z + 2 + ncol_z+1;

    int cc,kk, kk1;
    if (tr) {
      /* Loop over the columns of y and z */
      for (cc=0; cc<ncol_z; ++cc) {
        /* Get the dense column of y */
        for (kk=colind_y[cc]; kk<colind_y[cc+1]; ++kk) {
          w[row_y[kk]] = y[kk];
        }
        /* Loop over the nonzeros of z */
        for (kk=colind_z[cc]; kk<colind_z[cc+1]; ++kk) {
          int rr = row_z[kk];
          /* Loop over corresponding columns of x */
          for (kk1=colind_x[rr]; kk1<colind_x[rr+1]; ++kk1) {
            z[kk] += x[kk1] * w[row_x[kk1]];
          }
        }
      }
    } else {
      /* Loop over the columns of y and z */
      for (cc=0; cc<ncol_y; ++cc) {
        /* Get the dense column of z */
        for (kk=colind_z[cc]; kk<colind_z[cc+1]; ++kk) {
          w[row_z[kk]] = z[kk];
        }
        /* Loop over the nonzeros of y */
        for (kk=colind_y[cc]; kk<colind_y[cc+1]; ++kk) {
          int rr = row_y[kk];
          /* Loop over corresponding columns of x */
          for (kk1=colind_x[rr]; kk1<colind_x[rr+1]; ++kk1) {
            w[row_x[kk1]] += x[kk1]*y[kk];
          }
        }
        /* Get the sparse column of z */
        for (kk=colind_z[cc]; kk<colind_z[cc+1]; ++kk) {
          z[kk] = w[row_z[kk]];
        }
      }
    }
  }

  template<typename real_t>
  void CASADI_PREFIX(mv)(const real_t* x, const int* sp_x, const real_t* y, real_t* z, int tr) {
    /* Get sparsities */
    int ncol_x = sp_x[1];
    const int *colind_x = sp_x+2, *row_x = sp_x + 2 + ncol_x+1;
    int i, el;
    if (tr) {
      /* loop over the columns of x */
      for (i=0; i<ncol_x; ++i) {
        /* loop over the non-zeros of x */
        for (el=colind_x[i]; el<colind_x[i+1]; ++el) {
          z[i] += x[el] * y[row_x[el]];
        }
      }
    } else {
      /* loop over the columns of x */
      for (i=0; i<ncol_x; ++i) {
        /* loop over the non-zeros of x */
        for (el=colind_x[i]; el<colind_x[i+1]; ++el) {
          z[row_x[el]] += x[el] * y[i];
        }
      }
    }
  }

  template<typename real_t>
  void CASADI_PREFIX(trans)(const real_t* x, const int* sp_x, real_t* y, const int* sp_y, int *tmp) {
    int ncol_x = sp_x[1];
    int nnz_x = sp_x[2 + ncol_x];
    const int* row_x = sp_x + 2 + ncol_x+1;
    int ncol_y = sp_y[1];
    const int* colind_y = sp_y+2;
    int k;
    for (k=0; k<ncol_y; ++k) tmp[k] = colind_y[k];
    for (k=0; k<nnz_x; ++k) {
      y[tmp[row_x[k]]++] = x[k];
    }
  }

  template<typename real_t>
  real_t CASADI_PREFIX(norm_1)(int n, const real_t* x) {
    real_t ret = 0;
    int i;
    if (x) {
      for (i=0; i<n; ++i) ret += fabs(*x++);
    }
    return ret;
  }

  template<typename real_t>
  real_t CASADI_PREFIX(norm_2)(int n, const real_t* x) {
    return sqrt(CASADI_PREFIX(dot)(n, x, x));
  }

  template<typename real_t>
  real_t CASADI_PREFIX(norm_inf)(int n, const real_t* x) {
    real_t ret = 0;
    int i;
    for (i=0; i<n; ++i) ret = fmax(ret, fabs(*x++));
    return ret;
  }

  template<typename real_t>
  real_t CASADI_PREFIX(norm_inf_mul)(const real_t* x, const int* sp_x, const real_t* y, const int* sp_y, real_t* dwork, int* iwork) {
    real_t res = 0;
    /* Get sparsities */
    int nrow_x = sp_x[0], ncol_x = sp_x[1];
    const int *colind_x = sp_x+2, *row_x = sp_x + 2 + ncol_x+1;
    int ncol_y = sp_y[1];
    const int *colind_y = sp_y+2, *row_y = sp_y + 2 + ncol_y+1;

    /* Implementation borrowed from Scipy's sparsetools/csr.h */
    /* method that uses O(n) temp storage */
    int *mask = iwork + ncol_y+1;

    int i,jj,kk;
    // Pass 1
    for (i=0; i<nrow_x; ++i) mask[i] = -1;
    iwork[0] = 0;
    int nnz = 0;
    for (i=0; i<ncol_y; ++i) {
      int row_nnz = 0;
      for (jj=colind_y[i]; jj < colind_y[i+1]; jj++) {
        int j = row_y[jj];
        for (kk=colind_x[j]; kk < colind_x[j+1]; kk++) {
          int k = row_x[kk];
          if (mask[k] != i) {
            mask[k] = i;
            row_nnz++;
          }
        }
      }
      int next_nnz = nnz + row_nnz;
      nnz = next_nnz;
      iwork[i+1] = nnz;
    }

    // Pass 2
    int *next = iwork + ncol_y+1;
    for (i=0; i<nrow_x; ++i) next[i] = -1;
    real_t* sums = dwork;
    for (i=0; i<nrow_x; ++i) sums[i] = 0;
    nnz = 0;
    iwork[0] = 0;
    for (i=0; i<ncol_y; ++i) {
      int head   = -2;
      int length =  0;
      int jj_start = colind_y[i];
      int jj_end   = colind_y[i+1];
      for (jj=jj_start; jj<jj_end; ++jj) {
        int j = row_y[jj];
        real_t v = y[jj];
        int kk_start = colind_x[j];
        int kk_end   = colind_x[j+1];
        for (kk = kk_start; kk<kk_end; ++kk) {
          int k = row_x[kk];
          sums[k] += v*x[kk];
          if (next[k] == -1) {
            next[k] = head;
            head  = k;
            length++;
          }
        }
      }
      for (jj=0; jj<length; ++jj) {
        if (!is_zero(sums[head])) {
          real_t a = fabs(sums[head]);
          res = fmax(res, a);
          nnz++;
        }
        int temp = head;
        head = next[head];
        next[temp] = -1; //clear arrays
        sums[temp] =  0;
      }
      iwork[i+1] = nnz;
    }
    return res;
  }

  template<typename real_t>
  real_t CASADI_PREFIX(bilin)(const real_t* A, const int* sp_A, const real_t* x, const real_t* y) {
    /* Get sparsities */
    int ncol_A = sp_A[1];
    const int *colind_A = sp_A+2, *row_A = sp_A + 2 + ncol_A+1;

    /* Return value */
    real_t ret=0;

    /* Loop over the columns of A */
    int cc, rr, el;
    for (cc=0; cc<ncol_A; ++cc) {
      /* Loop over the nonzeros of A */
      for (el=colind_A[cc]; el<colind_A[cc+1]; ++el) {
        /* Get the row */
        rr=row_A[el];

        /* Add contribution */
        ret += x[rr]*A[el]*y[cc];
      }
    }
    return ret;
  }

  template<typename real_t>
  void CASADI_PREFIX(rank1)(real_t* A, const int* sp_A, real_t alpha, const real_t* x, const real_t* y) {
    /* Get sparsities */
    int ncol_A = sp_A[1];
    const int *colind_A = sp_A+2, *row_A = sp_A + 2 + ncol_A+1;

    /* Loop over the columns of A */
    int cc, rr, el;
    for (cc=0; cc<ncol_A; ++cc) {
      /* Loop over the nonzeros of A */
      for (el=colind_A[cc]; el<colind_A[cc+1]; ++el) {
        /* Get row */
        rr = row_A[el];

        /* Add the multiple */
        A[el] += alpha*x[rr]*y[cc];
      }
    }
  }

  template<typename real_t>
  void CASADI_PREFIX(getu)(const real_t* x, const int* sp_x, real_t* v) {
    /* Get sparsities */
    int ncol_x = sp_x[1];
    const int *colind_x = sp_x+2, *row_x = sp_x + 2 + ncol_x+1;

    /* Loop over the columns of x */
    int cc, el;
    for (cc=0; cc<ncol_x; ++cc) {
      /* Loop over the nonzeros of x */
      for (el=colind_x[cc]; el<colind_x[cc+1] && row_x[el]<=cc; ++el) {
        *v++ = x[el];
      }
    }
  }

  template<typename real_t>
  real_t CASADI_PREFIX(polyval)(const real_t* p, int n, real_t x) {
    real_t r=p[0];
    int i;
    for (i=1; i<=n; ++i) {
      r = r*x + p[i];
    }
    return r;
  }

  inline
  int CASADI_PREFIX(flip)(int* corner, int ndim) {
    int i;
    for (i=0; i<ndim; ++i) {
      if (corner[i]) {
        corner[i]=0;
      } else {
        corner[i]=1;
        return 1;
      }
    }
    return 0;
  }

  template<typename real_t>
  int CASADI_PREFIX(low)(real_t x, const double* grid, int ng, int lookup_mode) {
    if (lookup_mode) {
      double g0 = grid[0];
      int ret = (int) ((x-g0)*(ng-1)/(grid[ng-1]-g0));
      if (ret<0) ret=0;
      if (ret>ng-2) ret=ng-2;
      return ret;
    } else {
      int i;
      for (i=0; i<ng-2; ++i) {
        if (x < grid[i+1]) break;
      }
      return i;
    }
  }

  template<typename real_t>
  void CASADI_PREFIX(interpn_weights)(int ndim, const real_t* grid, const int* offset, const real_t* x, real_t* alpha, int* index, const int* lookup_mode) {
    /* Left index and fraction of interval */
    int i;
    for (i=0; i<ndim; ++i) {
      /* Grid point */
      real_t xi = x ? x[i] : 0;
      /* Grid */
      const real_t* g = grid + offset[i];
      int ng = offset[i+1]-offset[i];
      /* Find left index */
      int j = index[i] = CASADI_PREFIX(low)(xi, g, ng, lookup_mode[i]);
      /* Get interpolation/extrapolation alpha */
      alpha[i] = (xi-g[j])/(g[j+1]-g[j]);
    }
  }

  template<typename real_t>
  real_t CASADI_PREFIX(interpn_interpolate)(int ndim, const int* offset, const real_t* values, const real_t* alpha, const int* index, const int* corner, real_t* coeff) {
    /* Get weight and value for corner */
    real_t c=1;
    int ld=1; /* leading dimension */
    int i;
    for (i=0; i<ndim; ++i) {
      if (coeff) *coeff++ = c;
      if (corner[i]) {
        c *= alpha[i];
      } else {
        c *= 1-alpha[i];
      }
      values += (index[i]+corner[i])*ld;
      ld *= offset[i+1]-offset[i];
    }
    if (coeff) {
      return *values;
    } else {
      return c**values;
    }
  }

  template<typename real_t>
  real_t CASADI_PREFIX(interpn)(int ndim, const real_t* grid, const int* offset, const real_t* values, const real_t* x, const int* lookup_mode, int* iw, real_t* w) {
    /* Work vectors */
    real_t* alpha = w; w += ndim;
    int* index = iw; iw += ndim;
    int* corner = iw; iw += ndim;
    /* Left index and fraction of interval */
    CASADI_PREFIX(interpn_weights)(ndim, grid, offset, x, alpha, index, lookup_mode);
    /* Loop over all corners, add contribution to output */
    CASADI_PREFIX(fill)(corner, ndim, 0);
    real_t ret = 0;
    do {
      real_t* coeff = 0;
      ret += CASADI_PREFIX(interpn_interpolate)(ndim, offset, values,
        alpha, index, corner, coeff);
    } while (CASADI_PREFIX(flip)(corner, ndim));
    return ret;
  }

  template<typename real_t>
  void CASADI_PREFIX(interpn_grad)(real_t* grad, int ndim, const real_t* grid, const int* offset, const real_t* values, const real_t* x, const int* lookup_mode, int* iw, real_t* w) {
    /* Quick return */
    if (!grad) return;
    /* Work vectors */
    real_t* alpha = w; w += ndim;
    real_t* coeff = w; w += ndim;
    int* index = iw; iw += ndim;
    int* corner = iw; iw += ndim;
    /* Left index and fraction of interval */
    CASADI_PREFIX(interpn_weights)(ndim, grid, offset, x, alpha, index, lookup_mode);
    /* Loop over all corners, add contribution to output */
    CASADI_PREFIX(fill)(corner, ndim, 0);
    CASADI_PREFIX(fill)(grad, ndim, 0.);
    do {
      /* Get coefficients */
      real_t v = CASADI_PREFIX(interpn_interpolate)(ndim, offset, values,
        alpha, index, corner, coeff);
      /* Propagate to alpha */
      int i;
      for (i=ndim-1; i>=0; --i) {
        if (corner[i]) {
          grad[i] += v*coeff[i];
          v *= alpha[i];
        } else {
          grad[i] -= v*coeff[i];
          v *= 1-alpha[i];
        }
      }
    } while (CASADI_PREFIX(flip)(corner, ndim));
    /* Propagate to x */
    int i;
    for (i=0; i<ndim; ++i) {
      const real_t* g = grid + offset[i];
      int j = index[i];
      grad[i] /= g[j+1]-g[j];
    }
  }

  template<typename real_t>
  void CASADI_PREFIX(de_boor)(real_t x, const real_t* knots, int n_knots, int degree, real_t* boor) {
    // length boor: n_knots-1
    for (int d=1;d<degree+1;++d) {
      for (int i=0;i<n_knots-d-1;++i) {
        real_t b = 0;
        real_t bottom = knots[i + d] - knots[i];
        if (bottom) b = (x - knots[i]) * boor[i] / bottom;
        bottom = knots[i + d + 1] - knots[i + 1];
        if (bottom) b += (knots[i + d + 1] - x) * boor[i + 1] / bottom;
        boor[i] = b;
      }
    }
  }


  template<typename real_t>
  void CASADI_PREFIX(nd_boor_eval)(real_t* ret, int n_dims, const real_t* all_knots, const int* offset, const int* all_degree, const int* strides, const real_t* c, int m, const real_t* all_x, const int* lookup_mode, int reverse, int* iw, real_t* w) {
    int* boor_offset = iw; iw+=n_dims+1;
    int* starts = iw; iw+=n_dims;
    int* index = iw; iw+=n_dims;
    int* coeff_offset = iw;

    real_t* cumprod = w; w+= n_dims+1;
    real_t* all_boor = w;

    boor_offset[0] = 0;
    cumprod[n_dims] = 1;
    coeff_offset[n_dims] = 0;

    int n_iter = 1;
    for (int k=0;k<n_dims;++k) {
      real_t* boor = all_boor+boor_offset[k];

      int degree = all_degree[k];
      const real_t* knots = all_knots + offset[k];
      int n_knots = offset[k+1]-offset[k];
      int n_b = n_knots-degree-1;

      real_t x = all_x[k];
      int L = CASADI_PREFIX(low)(x, knots+degree, n_knots-2*degree, lookup_mode[k]);

      int start = L;
      if (start>n_b-degree-1) start = n_b-degree-1;

      starts[k] = start;

      CASADI_PREFIX(fill)(boor, 2*degree+1, 0.0);
      if (x>=knots[0] && x<=knots[n_knots-1]) {
        if (x==knots[1]) {
          CASADI_PREFIX(fill)(boor, degree+1, 1.0);
        } else if (x==knots[n_knots-1]) {
          boor[degree] = 1;
        } else if (knots[L+degree]==x) {
          boor[degree-1] = 1;
        } else {
          boor[degree] = 1;
        }
      }
      CASADI_PREFIX(de_boor)(x, knots+start, 2*degree+2, degree, boor);
      boor+= degree+1;
      n_iter*= degree+1;
      boor_offset[k+1] = boor_offset[k] + degree+1;
    }

    CASADI_PREFIX(fill)(index, n_dims, 0);

    // Prepare cumulative product
    for (int pivot=n_dims-1;pivot>=0;--pivot) {
      cumprod[pivot] = (*(all_boor+boor_offset[pivot]))*cumprod[pivot+1];
      coeff_offset[pivot] = starts[pivot]*strides[pivot]+coeff_offset[pivot+1];
    }

    for (int k=0;k<n_iter;++k) {

      // accumulate result
      for (int i=0;i<m;++i) {
        if (reverse) {
          ret[coeff_offset[0]+i] += c[i]*cumprod[0];
        } else {
          ret[i] += c[coeff_offset[0]+i]*cumprod[0];
        }
      }

      // Increment index
      index[0]++;
      int pivot = 0;

      // Handle index overflow
      {
        // increment next index (forward)
        while (index[pivot]==boor_offset[pivot+1]-boor_offset[pivot]) {
          index[pivot] = 0;
          if (pivot==n_dims-1) break;
          index[++pivot]++;
        }

        // update cumulative structures (reverse)
        while (pivot>0) {
          // Compute product
          cumprod[pivot] = (*(all_boor+boor_offset[pivot]+index[pivot]))*cumprod[pivot+1];
          // Compute offset
          coeff_offset[pivot] = (starts[pivot]+index[pivot])*strides[pivot]+coeff_offset[pivot+1];
          pivot--;
        }
      }

      // Compute product
      cumprod[0] = (*(all_boor+index[0]))*cumprod[1];

      // Compute offset
      coeff_offset[0] = (starts[0]+index[0])*m+coeff_offset[1];

    }
  }

} // namespace casadi



/// \endcond

#endif // CASADI_CASADI_RUNTIME_HPP
