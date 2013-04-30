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

#include "runtime.h"
#ifdef __cplusplus
extern "C" {
#endif

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
#ifdef __cplusplus
} // extern "C"
#endif
