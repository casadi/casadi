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

void casadi_trans(const d* x, const int* sp_x, d* y, const int* sp_y, int *tmp){
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
