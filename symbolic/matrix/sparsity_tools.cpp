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
#include "../stl_vector_tools.hpp"
#include "../matrix/matrix.hpp"
#include "../matrix/matrix_tools.hpp"

using namespace std;

namespace CasADi{
  
  CRSSparsity sp_dense(int n, int m) {
    return CRSSparsity(n,m,true);
  }

  CRSSparsity sp_sparse(int n, int m) {
    return CRSSparsity(n,m,false);
  }

  CRSSparsity sp_dense(const std::pair<int,int> &nm) {
    return CRSSparsity(nm.first,nm.second,true);
  }

  CRSSparsity sp_sparse(const std::pair<int,int> &nm) {
    return CRSSparsity(nm.first,nm.second,false);
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

  CRSSparsity sp_diag(int n){
    casadi_assert_message(n>=0, "sp_diag expects a positive integer as argument");
  
    // Construct sparsity pattern
    std::vector<int> col(n);
    std::vector<int> rowind(n+1);
    int el = 0;
    for(int i=0; i<n; ++i){
      rowind[i] = el;
      col[el++] = i;
    }
    rowind.back() = el;

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


  CRSSparsity sp_rowcol(const std::vector<int>& row, const std::vector<int>& col, int nrow, int ncol) {
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
      casadi_error(
		   "sp_rowcol: out-of-range error." << endl <<
		   "The " << k << "th entry of row (" << row[k] << ") was bigger or equal to the specified total number of rows (" << nrow << ")"
		   );
    }
    return CRSSparsity(nrow, ncol, col_new, rowind);
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
    casadi_assert_message(a.numel() == n*m,
			  "reshape: number of elements must remain the same." << endl <<
			  "Input argument has shape " << a.size1() << " x " << a.size2() << " =  " << a.numel() << ", while you request a reshape to " <<
			  n << " x " << m << " =  " << n*m
			  );

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
  
    return  sp_triplet(n,m,row_new,col_new);
  }

  CRSSparsity vec(const CRSSparsity& a){
    return reshape(trans(a),a.numel(),1);
  }

  CRSSparsity trans(const CRSSparsity& a) {
    return a.transpose();
  }


  CRSSparsity lowerSparsity(const CRSSparsity& a, bool includeDiagonal) {
    const std::vector<int> & col= a.col();
    std::vector<int> row = a.getRow();
  

    std::vector<int> new_col;   // new col
    std::vector<int> new_row;   // new row
  
    int offset = (includeDiagonal ? 1 : 0);
  
    // Find size
    int n=0;
    for (int k=0;k<row.size();k++) n+= row[k] + offset > col[k];
    new_col.resize(n);
    new_row.resize(n);
  
    // populate return vector
    int cnt=0;
    for (int k=0;k<row.size();k++) {
      if (row[k] + offset > col[k]) {
	new_col[cnt]=col[k];
	new_row[cnt]=row[k];
	cnt++;
      }
    }
    return sp_triplet(a.size1(), a.size2(), new_row, new_col);
  
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

  CRSSparsity sp_triplet(int nrow, int ncol, const std::vector<int>& row, const std::vector<int>& col, std::vector<int>& mapping, bool invert_mapping){
    // Assert dimensions
    casadi_assert_message(row.size()==col.size(),"inconsistent lengths");

    // Create the return sparsity pattern and access vectors
    CRSSparsity ret(nrow,ncol);
    vector<int> &r_rowind = ret.rowindRef();
    vector<int> &r_col = ret.colRef();
    r_col.reserve(col.size());

    // Consistency check and check if elements are already perfectly ordered with no duplicates
    int last_row=-1, last_col=-1;
    bool perfectly_ordered=true;
    for(int k=0; k<row.size(); ++k){
      // Consistency check
      casadi_assert_message(row[k]>=0 && row[k]<nrow,"Row index out of bounds");
      casadi_assert_message(col[k]>=0 && col[k]<ncol,"Col index out of bounds");
    
      // Check if ordering is already perfect
      perfectly_ordered = perfectly_ordered && (row[k]<last_row || (row[k]==last_row && col[k]<=last_col));
      last_row = row[k];
      last_col = col[k];
    }
  
    // Quick return if perfectly ordered
    if(perfectly_ordered){
      // Save columns
      r_col.resize(col.size());
      copy(col.begin(),col.end(),r_col.begin());
    
      // Find offset index
      int el=0;
      for(int i=0; i<nrow; ++i){
	while(el<row.size() && row[el]==i) el++; 
	r_rowind[i+1] = el;
      }
    
      // Identity mapping
      mapping.resize(col.size());
      for(int k=0; k<col.size(); ++k)
	mapping[k] = k;
    
      // Quick return
      return ret;
    }
    
    // Reuse data
    vector<int>& mapping1 = invert_mapping ? r_col : mapping;
    vector<int>& mapping2 = invert_mapping ? mapping : r_col;
  
    // Make sure that enough memory is allocated to use as a work vector
    mapping1.reserve(std::max(ncol+1,int(row.size())));
  
    // Number of elements in each column
    vector<int>& colcount = mapping1; // reuse memory
    colcount.resize(ncol+1);
    fill(colcount.begin(),colcount.end(),0);
    for(vector<int>::const_iterator it=col.begin(); it!=col.end(); ++it){
      colcount[*it+1]++;
    }
  
    // Cumsum to get index offset for each column
    for(int i=0; i<ncol; ++i){
      colcount[i+1] += colcount[i];
    }
  
    // New column for each old column
    mapping2.resize(col.size());
    for(int k=0; k<col.size(); ++k){
      mapping2[colcount[col[k]]++] = k;
    }
  
    // Number of elements in each row
    vector<int>& rowcount = r_rowind; // reuse memory, r_rowind is already the right size and is filled with zeros
    for(vector<int>::const_iterator it=mapping2.begin(); it!=mapping2.end(); ++it){
      rowcount[row[*it]+1]++;
    }
  
    // Cumsum to get index offset for each row
    for(int i=0; i<nrow; ++i){
      rowcount[i+1] += rowcount[i];
    }

    // New row for each old row
    mapping1.resize(row.size());
    for(vector<int>::const_iterator it=mapping2.begin(); it!=mapping2.end(); ++it){
      mapping1[rowcount[row[*it]]++] = *it;
    }

    // Current element in the return matrix
    int r_el = 0;
    r_col.resize(row.size());

    // Current nonzero
    vector<int>::const_iterator it=mapping1.begin();

    // Loop over rows
    r_rowind[0] = 0;
    for(int i=0; i<nrow; ++i){
    
      // Previous column (to detect duplicates)
      int j_prev = -1;
    
      // Loop over nonzero elements of the row
      while(it!=mapping1.end() && row[*it]==i){

	// Get the element
	int el = *it;
	it++;

	// Get the column
	int j = col[el];
      
	// If not a duplicate, save to return matrix
	if(j!=j_prev)
	  r_col[r_el++] = j;
      
	if(invert_mapping){
	  // Save to the inverse mapping
	  mapping2[el] = r_el-1;        
	} else {
	  // If not a duplicate, save to the mapping vector
	  if(j!=j_prev)
	    mapping1[r_el-1] = el;
	}
      
	// Save column
	j_prev = j;
      }
    
      // Update row offset
      r_rowind[i+1] = r_el;
    }

    // Resize the column vector
    r_col.resize(r_el);

    // Resize mapping matrix
    if(!invert_mapping){
      mapping1.resize(r_el);
    }
  
    return ret;
  }


  CRSSparsity sp_triplet(int n, int m, const std::vector<int>& row, const std::vector<int>& col){
    std::vector<int> mapping;
    return sp_triplet(n,m,row,col,mapping,false);
  }


  CRSSparsity mul(const  CRSSparsity& a, const  CRSSparsity &b) {
    return (mul(DMatrix(a,1),DMatrix(b,1))).sparsity();
  }

  std::size_t hash_sparsity(int nrow, int ncol, const std::vector<int>& col, const std::vector<int>& rowind){
    // Condense the sparsity pattern to a single, deterministric number
    std::size_t ret=0;
    hash_combine(ret,nrow);
    hash_combine(ret,ncol);
    hash_combine(ret,rowind);
    hash_combine(ret,col);
    return ret;
  }
  
  std::vector<int> sp_compress(const CRSSparsity& a){
    // Get the sparsity pattern
    int nrow = a.size1();
    int ncol = a.size2();
    const vector<int>& rowind = a.rowind();
    const vector<int>& col = a.col();
    
    // Create compressed pattern
    vector<int> ret;
    ret.reserve(1 + 1 + rowind.size() + col.size());
    ret.push_back(nrow);
    ret.push_back(ncol);
    ret.insert(ret.end(),rowind.begin(),rowind.end());
    ret.insert(ret.end(),col.begin(),col.end());
    return ret;
  }
  
  CRSSparsity sp_compress(const std::vector<int>& v){
    // Check consistency
    casadi_assert(v.size() >= 2);
    int nrow = v[0];
    //int ncol = v[1];
    casadi_assert(v.size() >= 2 + nrow+1);
    int nnz = v[2 + nrow];
    casadi_assert(v.size() == 2 + nrow+1 + nnz);

    // Call array version
    return sp_compress(&v.front());
  }
  
  CRSSparsity sp_compress(const int* v){
    // Get sparsity pattern
    int nrow = v[0];
    int ncol = v[1];
    const int *rowind = v+2;
    int nnz = rowind[nrow];
    const int *col = v + 2 + nrow+1;
    
    // Construct sparsity pattern
    return CRSSparsity(nrow, ncol, vector<int>(col,col+nnz), vector<int>(rowind,rowind+nrow+1));
  }
  
  int rank(const CRSSparsity& a) {
    std::vector< int > rowperm;
    std::vector< int > colperm;
    std::vector< int > rowblock;
    std::vector< int > colblock;
    std::vector< int > coarse_rowblock;
    std::vector< int > coarse_colblock;
    a.dulmageMendelsohn(rowperm, colperm, rowblock, colblock, coarse_rowblock, coarse_colblock);
    return coarse_colblock.at(3);
  }

  bool isSingular(const CRSSparsity& a) {
    casadi_assert_message(a.size1()==a.size2(),"isSingular: only defined for square matrices, but got " << a.dimString());
    return rank(a)!=a.size1();
  }

  CRSSparsity sp_unit(int n, int el){
    CRSSparsity ret = sp_sparse(n,1);
    ret.getNZ(el,0);
    return ret;
  }


} // namespace CasADi

