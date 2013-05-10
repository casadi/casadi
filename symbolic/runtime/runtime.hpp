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

#ifndef CASADI_RUNTIME_HPP
#define CASADI_RUNTIME_HPP


#include <math.h>

namespace CasADi {

  /// COPY: y <-x
  template<typename real_t>
  void casadi_copy(int n, const real_t* x, int inc_x, real_t* y, int inc_y);
  
  /// SWAP: x <-> y
  template<typename real_t>
  void casadi_swap(int n, real_t* x, int inc_x, real_t* y, int inc_y);
  
  /// COPY sparse: y <- x
  template<typename real_t>
  void casadi_copy_sparse(const real_t* x, const int* sp_x, real_t* y, const int* sp_y);
  
  /// SCAL: x <- alpha*x
  template<typename real_t>
  void casadi_scal(int n, real_t alpha, real_t* x, int inc_x);
  
  /// AXPY: y <- a*x + y
  template<typename real_t>
  void casadi_axpy(int n, real_t alpha, const real_t* x, int inc_x, real_t* y, int inc_y);
  
  /// DOT: inner_prod(x,y) -> return
  template<typename real_t>
  real_t casadi_dot(int n, const real_t* x, int inc_x, const real_t* y, int inc_y);
  
  /// ASUM: ||x||_1 -> return
  template<typename real_t>
  real_t casadi_asum(int n, const real_t* x, int inc_x);
  
  /// IAMAX: index corresponding to the entry with the largest absolute value
  template<typename real_t>
  int casadi_iamax(int n, const real_t* x, int inc_x);
  
  /// FILL: x <- alpha
  template<typename real_t>
  void casadi_fill(int n, real_t alpha, real_t* x, int inc_x);
  
  /// Sparse matrix-matrix multiplication, the second argument is transposed: z <- z + x*y'
  template<typename real_t>
  void casadi_mm_nt_sparse(const real_t* x, const int* sp_x, const real_t* trans_y, const int* sp_trans_y, real_t* z, const int* sp_z);
  
  /// NRM2: ||x||_2 -> return
  template<typename real_t>
  real_t casadi_nrm2(int n, const real_t* x, int inc_x);
  
  /// TRANS: y <- trans(x)
  template<typename real_t>
  void casadi_trans(const real_t* x, const int* sp_x, real_t* y, const int* sp_y, int *tmp);

}

// Implementations


// Note: due to restrictions of cmake IO processing, make sure that
//      1)   semicolons (;) are never immediately preceded by whitespace
//      2)   line continuation slashes (\) are always immediately preceded by whitespace
#define CASADI_GEMM_NT(M,N,K,A,LDA,B,LDB,C,LDC) \
  for(i=0, rr=C; i<M; ++i) \
    for(j=0; j<N; ++j, ++rr) \
      for(k=0, ss=A+i*LDA, tt=B+j*LDB; k<K; ++k) \
        *rr += *ss++**tt++;


namespace CasADi {

  template<typename real_t>
  void casadi_copy(int n, const real_t* x, int inc_x, real_t* y, int inc_y){
    int i;
    for(i=0; i<n; ++i){
      *y = *x;
      x += inc_x;
      y += inc_y;
    }
  }
  
  template<typename real_t>
  void casadi_swap(int n, real_t* x, int inc_x, real_t* y, int inc_y){
    real_t t;
    int i;
    for(i=0; i<n; ++i){
      t = *x;
      *x = *y;
      *y = t;
      x += inc_x;
      y += inc_y;
    }
  }
  
  template<typename real_t>
  void casadi_copy_sparse(const real_t* x, const int* sp_x, real_t* y, const int* sp_y){
    int nrow_x = sp_x[0];
    int ncol_x = sp_x[1];
    const int* rowind_x = sp_x+2;
    const int* col_x = sp_x + 2 + nrow_x+1;
    int nnz_x = rowind_x[nrow_x];
    int nrow_y = sp_y[0];
    int ncol_y = sp_y[1];
    const int* rowind_y = sp_y+2;
    const int* col_y = sp_y + 2 + nrow_y+1;
    int nnz_y = rowind_y[nrow_y];
    if(sp_x==sp_y){
      casadi_copy(nnz_x,x,1,y,1);
    } else {
      int i;
      for(i=0; i<nrow_x; ++i){
        int el_x = rowind_x[i];
        int el_x_end = rowind_x[i+1];
        int j_x = el_x<el_x_end ? col_x[el_x] : ncol_x;
        int el_y;
        for(el_y=rowind_y[i]; el_y!=rowind_y[i+1]; ++el_y){
          int j=col_y[el_y];
          while(j_x<j){
            el_x++;
            j_x = el_x<el_x_end ? col_x[el_x] : ncol_x;
          }
          if(j_x==j){
            y[el_y] = x[el_x++];
            j_x = el_x<el_x_end ? col_x[el_x] : ncol_x;
          } else {
            y[el_y] = 0;
          }
        }
      }
    }
  }
  
  template<typename real_t>
  void casadi_scal(int n, real_t alpha, real_t* x, int inc_x){
    int i;
    for(i=0; i<n; ++i){
      *x *= alpha;
      x += inc_x;
    }
  }
  
  template<typename real_t>
  void casadi_axpy(int n, real_t alpha, const real_t* x, int inc_x, real_t* y, int inc_y){
    int i;
    for(i=0; i<n; ++i){
      *y += alpha**x;
      x += inc_x;
      y += inc_y;
    }
  }
  
  template<typename real_t>
  real_t casadi_dot(int n, const real_t* x, int inc_x, const real_t* y, int inc_y){
    real_t r = 0;
    int i;
    for(i=0; i<n; ++i){
      r += *x**y;
      x += inc_x;
      y += inc_y;
    }
    return r;
  }
  
  template<typename real_t>
  real_t casadi_asum(int n, const real_t* x, int inc_x){
    real_t r = 0;
    int i;
    for(i=0; i<n; ++i){
      r += fabs(*x);
      x += inc_x;
    }
    return r;
  }
  
  template<typename real_t>
  int casadi_iamax(int n, const real_t* x, int inc_x){
    real_t t;
    real_t largest_value = -1.0;
    int largest_index = -1;
    int i;
    for(i=0; i<n; ++i){
      t = fabs(*x);
      x += inc_x;
      if(t>largest_value){
        largest_value = t;
        largest_index = i;
      }
    }
    return largest_index;
  }
  
  template<typename real_t>
  void casadi_fill(int n, real_t alpha, real_t* x, int inc_x){
    int i;
    for(i=0; i<n; ++i){
      *x = alpha;
      x += inc_x;
    }
  }
  
  template<typename real_t>
  void casadi_mm_nt_sparse(const real_t* x, const int* sp_x, const real_t* trans_y, const int* sp_trans_y, real_t* z, const int* sp_z){

    int nrow_x = sp_x[0];
    int ncol_x = sp_x[1];
    const int* rowind_x = sp_x+2;
    const int* col_x = sp_x + 2 + nrow_x+1;
    int nnz_x = rowind_x[nrow_x];

    int ncol_y = sp_trans_y[0];
    int nrow_y = sp_trans_y[1];
    const int* colind_y = sp_trans_y+2;
    const int* row_y = sp_trans_y + 2 + ncol_y+1;
    int nnz_y = colind_y[ncol_y];

    int nrow_z = sp_z[0];
    int ncol_z = sp_z[1];
    const int* rowind_z = sp_z+2;
    const int* col_z = sp_z + 2 + nrow_z+1;
    int nnz_z = rowind_z[nrow_z];

    int i;
    for(i=0; i<nrow_z; ++i){
      int el;
      for(el=rowind_z[i]; el<rowind_z[i+1]; ++el){
        int j = col_z[el];
        int el1 = rowind_x[i];
        int el2 = colind_y[j];
        while(el1 < rowind_x[i+1] && el2 < colind_y[j+1]){ 
          int j1 = col_x[el1];
          int i2 = row_y[el2];
          if(j1==i2){
            z[el] += x[el1++] * trans_y[el2++];
          } else if(j1<i2) {
            el1++;
          } else {
            el2++;
          }
        }
      }
    }
  }
  
  template<typename real_t>
  real_t casadi_nrm2(int n, const real_t* x, int inc_x){
    real_t r = 0;
    int i;
    for(i=0; i<n; ++i){
      r += *x**x;
      x += inc_x;
    }
    return sqrt(r);
  }
  
  template<typename real_t>
  void casadi_trans(const real_t* x, const int* sp_x, real_t* y, const int* sp_y, int *tmp){
    int nrow_x = sp_x[0];
    int nnz_x = sp_x[2 + nrow_x];
    const int* col_x = sp_x + 2 + nrow_x+1;
    int nrow_y = sp_y[0];
    const int* rowind_y = sp_y+2;
    int k;
    for(k=0; k<nrow_y; ++k) tmp[k] = rowind_y[k];
    for(k=0; k<nnz_x; ++k){
      y[tmp[col_x[k]]++] = x[k];
    }
  }
  
}


#endif // CASADI_RUNTIME_HPP
