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
  
} // namespace CasADi

