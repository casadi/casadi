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

#include "sparsity_tools.hpp"

namespace CasADi{
  
CRSSparsity sp_tril(int n) {
  int c=0;
  int t=0;
  std::vector< int >  	col((n*(n+1))/2,0);
  for (int i=0;i<(n*(n+1))/2;i++) {
    col[i]=t++;
    if (t>c) {
      t=0;
      c++;
    }
  }

  std::vector< int >  	rowind(n+1,0);
  c=0;
  for (int i=1;i<n+1;i++)
    rowind[i]=rowind[i-1]+1+(c++);

  return CRSSparsity(n,n,col,rowind);
}

CRSSparsity sp_diag(int n) {
  int c=0;
  int t=0;
  std::vector< int >  	col(n);
  for (int i=0;i<n;i++) {
    col[i]=i;
  }

  std::vector< int >  	rowind(n+1);
  std::copy(col.begin(),col.end(),rowind.begin());
  rowind[n]=n;

  return CRSSparsity(n,n,col,rowind);
}

CRSSparsity sp_rowcol(std::vector<int> row, std::vector<int> col, int nrow, int ncol) {
  std::vector<int> rowind(nrow+1);
  std::vector<int> col_new(row.size()*col.size());
  for (int i=0;i<row.size();i++)
    std::copy(col.begin(),col.end(),col_new.begin()+col.size()*i);
  int cnt=0;
  int z=0;
  rowind[0]=0;
  for (int k=0; k < row.size(); k++) {
    while (z<row[k])
      rowind[++z]=cnt;
    rowind[row[k]+1]=(cnt+=col.size());
  }
  while (z<nrow)
    rowind[++z]=cnt;

  return CRSSparsity(nrow, ncol, col_new, rowind);
}

CRSSparsity sp_NZ(std::vector<int> row, std::vector<int> col, int nrow, int ncol, bool monotone) {
  std::vector<int> rowind(nrow+1);
  int cnt=0;
  int z=0;
  rowind[0]=0;
  for (int k=0; k < row.size(); k++) {
    while (z<row[k])
      rowind[++z]=cnt;
    rowind[row[k]+1]=cnt++;
  }
  while (z<nrow)
    rowind[++z]=cnt;
  return CRSSparsity(nrow, ncol, col, rowind);
}

std::vector<int> getNZDense(const CRSSparsity &sp) {
  std::vector<int> ret(sp.size());
  std::vector<int> row = sp.getRow();
  const std::vector<int> &col = sp.col();
  int s2 = sp.size2();
  for(int k=0;k<sp.size();k++) {
    ret[k] = col[k]+row[k]*s2;
  }
  return ret;
}

CRSSparsity reshape(const CRSSparsity& a, int n, int m){
  casadi_assert_message(a.numel() == n*m, "resize: number of elements must remain the same");
  
  std::vector<int> row = a.getRow();
  const std::vector<int> &col = a.col();

  std::vector<int> row_new(a.size());
  std::vector<int> col_new(a.size());
  
  int i=0;int j=0; int z =0;
  for(int k=0; k<a.size(); ++k){
    int i = row[k];
    int j = col[k];
    z = j+i*a.size2();
    row_new[k] = z/m;
    col_new[k] = z%m;
  }
  
  return  sp_NZ(row_new,col_new,n,m,true);
}

CRSSparsity vec(const CRSSparsity& a){
  return reshape(a,a.numel(),1);
}
  
} // namespace CasADi

