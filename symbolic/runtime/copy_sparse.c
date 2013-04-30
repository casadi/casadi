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

void casadi_copy_sparse(const d* x, const int* sp_x, d* y, const int* sp_y){
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
