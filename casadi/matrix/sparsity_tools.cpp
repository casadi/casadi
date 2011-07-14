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
#include "../casadi_exception.hpp"

using namespace std;

namespace CasADi{
  
CRSSparsity sp_dense(int n, int m) {
  return CRSSparsity(n,m,true);
}

CRSSparsity sp_sparse(int n, int m) {
  return CRSSparsity(n,m,false);
}

CRSSparsity sp_tril(int n) {
  if (n<0)
    throw CasadiException("sp_tril expects a positive integer as argument");

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
  casadi_assert_message(n>=0, "sp_diag expects a positive integer as argument");
/*  int c=0;
  int t=0;*/
  std::vector< int >  	col(n);
  for (int i=0;i<n;i++) {
    col[i]=i;
  }

  std::vector< int >  	rowind(n+1);
  std::copy(col.begin(),col.end(),rowind.begin());
  rowind[n]=n;

  return CRSSparsity(n,n,col,rowind);
}

CRSSparsity sp_band(int n, int p) {
  casadi_assert_message(n>=0, "sp_band expects a positive integer as argument");
  casadi_assert_message((p<0? -p : p)<n, "sp_band: position of band schould be smaller then size argument");
  
  int nc = n-(p<0? -p : p);
  
  std::vector< int >  	col(nc);
  
  int offset = max(p,0);
  for (int i=0;i<nc;i++) {
    col[i]=i+offset;
  }
  
  std::vector< int >  	rowind(n+1);
  
  offset = min(p,0);
  for (int i=0;i<n+1;i++) {
    rowind[i]=max(min(i+offset,nc),0);
  }
  
  return CRSSparsity(n,n,col,rowind);
  
}

CRSSparsity sp_banded(int n, int p) {
  throw CasadiException("sp_banded: Not implemented yet");
}


CRSSparsity sp_rowcol(std::vector<int> row, std::vector<int> col, int nrow, int ncol) {
  std::vector<int> rowind(nrow+1);
  std::vector<int> col_new(row.size()*col.size());
  
  // resulting col: the entries of col are repeated row.size() times
  for (int i=0;i<row.size();i++)
    std::copy(col.begin(),col.end(),col_new.begin()+col.size()*i);

  // resulting rowind: first entry is always zero
  int cnt=0;
  int z=0;
  rowind[0]=0;
  int k=0;
  try {
    for (k=0; k < row.size(); k++) {
      // resulting rowind: fill up rowind entries with copies
      while (z<row[k])
        rowind.at(++z)=cnt;
        
      // resulting rowind: add col.size() at each row[k]
      rowind.at(row[k]+1)=(cnt+=col.size());
    }
    while (z<nrow)
      rowind.at(++z)=cnt;                 
  }
  catch (out_of_range& oor) {
    stringstream ss;
    ss << "sp_rowcol: out-of-range error." << endl;
    ss << "The " << k << "th entry of row (" << row[k] << ") was bigger or equal to the specified total number of rows (" << nrow << ")" << endl;
    throw CasadiException(ss.str());
  }
  return CRSSparsity(nrow, ncol, col_new, rowind);
}

CRSSparsity sp_NZ(std::vector<int> row, std::vector<int> col, int nrow, int ncol, bool monotone) {
  if (row.size()!=col.size()) {
    stringstream ss;
    ss << "sp_NZ: row and col vectors must be of same length." << endl;
    ss << "row is length " << row.size() << " and " << " col has length " << col.size() << endl;
    throw CasadiException(ss.str());
  }
  if (monotone == false)
    throw CasadiException("sp_NZ: Not implemented for monotone false");
  // the given col is fine, we just need to calculate rowind.
  std::vector<int> rowind(nrow+1);
  int cnt=0;  // cumulative non-zero counter
  int z=0;
  rowind[0]=0;
  int k=0;
  try {
    for (k=0; k < row.size(); k++) {
      // fill up rowind entries copies of the increment counter
      while (z<row[k])
        rowind.at(++z)=cnt;
      // fill in the cumulative counter at the next row and increment counter 
      rowind.at(row[k]+1)=cnt++;
    }
    while (z<nrow)
      rowind.at(++z)=cnt;
  }
  catch (out_of_range& oor) {
    stringstream ss;
    ss << "sp_NZ: out-of-range error." << endl;
    ss << "The " << k << "th entry of row (" << row[k] << ") was bigger or equal to the specified total number of rows (" << nrow << ")." << endl;
    throw CasadiException(ss.str());
  }
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
  if (a.numel() != n*m) {
    stringstream ss;
    ss << "reshape: number of elements must remain the same." << endl;
    ss << "Input argument has shape " << a.size1() << " x " << a.size2() << " =  " << a.numel() << ", while you request a reshape to ";
    ss << n << " x " << m << " =  " << n*m << endl;
    throw CasadiException(ss.str());
  }
  // our strategy is: (col,rowind) -> (col,row) -> modulus calculus -> (col_new, row_new) -> sp_NZ
  std::vector<int> row = a.getRow();
  const std::vector<int> &col = a.col();

  std::vector<int> row_new(a.size());
  std::vector<int> col_new(a.size());
  
//  int i=0;int j=0; int z =0;
  for(int k=0; k<a.size(); ++k){
    int i = row[k];
    int j = col[k];
    int z = j+i*a.size2();
    row_new[k] = z/m;
    col_new[k] = z%m;
  }
  
  return  sp_NZ(row_new,col_new,n,m,true);
}

CRSSparsity vec(const CRSSparsity& a){
  return reshape(a,a.numel(),1);
}


CRSSparsity lowerSparsity(const CRSSparsity& a) {
  const std::vector<int> & col= a.col();
  std::vector<int> row = a.getRow();
  

  std::vector<int> new_col;   // new col
  std::vector<int> new_row;   // new row
  
  // Find size
  int n=0;
  for (int k=0;k<row.size();k++) n+= row[k] >= col[k];
  new_col.resize(n);
  new_row.resize(n);
  
  // populate return vector
  int cnt=0;
  for (int k=0;k<row.size();k++) {
    if (row[k] >= col[k]) {
      new_col[cnt]=col[k];
      new_row[cnt]=row[k];
      cnt++;
    }
  }
  return sp_NZ(new_row, new_col,a.size1(), a.size2(), true);
  
}

std::vector<int> lowerNZ(const CRSSparsity& a) {
  const std::vector<int> & col= a.col();
  std::vector<int> row = a.getRow();
  
  // Return vector
  std::vector<int> ret;
  
  // Find size of return vector
  int n=0;
  for (int k=0;k<row.size();k++) n+= (row[k] >= col[k]);
  ret.resize(n);
  
  // populate return vector
  int cnt=0;
  for (int k=0;k<row.size();k++) {
    if (row[k] >= col[k]) ret[cnt++]=k;
  }
  
  return ret;
}
  
} // namespace CasADi

