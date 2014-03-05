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
  
  Sparsity sp_dense(int nrow, int ncol) {
    return Sparsity(nrow,ncol,true);
  }

  Sparsity sp_sparse(int nrow, int ncol) {
    return Sparsity(nrow,ncol,false);
  }

  Sparsity sp_dense(const std::pair<int,int> &nm) {
    return Sparsity(nm.first,nm.second,true);
  }

  Sparsity sp_sparse(const std::pair<int,int> &nm) {
    return Sparsity(nm.first,nm.second,false);
  }

  Sparsity sp_triu(int n) {
    casadi_assert_message(n>=0,"sp_triu expects a positive integer as argument");
    int nrow=n, ncol=n;
    std::vector<int> colind, row;
    colind.reserve(ncol+1);
    row.reserve((n*(n+1))/2);
    
    // Loop over columns
    colind.push_back(0);
    for(int cc=0; cc<ncol; ++cc){
      // Loop over rows for the upper triangular half
      for(int rr=0; rr<=cc; ++rr){
        row.push_back(rr);
      }
      colind.push_back(row.size());
    }

    // Return the pattern
    return Sparsity(nrow,ncol,colind,row);
  }

  Sparsity sp_tril(int n) {
    casadi_assert_message(n>=0,"sp_tril expects a positive integer as argument");
    int nrow=n, ncol=n;
    std::vector<int> colind, row;
    colind.reserve(ncol+1);
    row.reserve((n*(n+1))/2);
    
    // Loop over columns
    colind.push_back(0);
    for(int cc=0; cc<ncol; ++cc){
      // Loop over rows for the lower triangular half
      for(int rr=cc; rr<nrow; ++rr){
        row.push_back(rr);
      }
      colind.push_back(row.size());
    }

    // Return the pattern
    return Sparsity(nrow,ncol,colind,row);
  }

  Sparsity sp_diag(int n){
    casadi_assert_message(n>=0, "sp_diag expects a positive integer as argument");
  
    // Construct sparsity pattern
    std::vector<int> row(n);
    std::vector<int> colind(n+1);
    int el = 0;
    for(int i=0; i<n; ++i){
      colind[i] = el;
      row[el++] = i;
    }
    colind.back() = el;

    return Sparsity(n,n,colind,row);
  }

  Sparsity sp_band(int n, int p) {
    casadi_assert_message(n>=0, "sp_band expects a positive integer as argument");
    casadi_assert_message((p<0? -p : p)<n, "sp_band: position of band schould be smaller then size argument");
  
    int nc = n-(p<0? -p : p);
  
    std::vector< int >          row(nc);
  
    int offset = max(p,0);
    for (int i=0;i<nc;i++) {
      row[i]=i+offset;
    }
  
    std::vector< int >          colind(n+1);
  
    offset = min(p,0);
    for (int i=0;i<n+1;i++) {
      colind[i]=max(min(i+offset,nc),0);
    }
  
    return Sparsity(n,n,colind,row);
  
  }

  Sparsity sp_banded(int n, int p) {
    // This is not an efficient implementation
    Sparsity ret = Sparsity(n,n);
    for (int i=-p;i<=p;++i) {
      ret = ret + sp_band(n,i);
    }
    return ret;
  }

  Sparsity sp_rowcol(const std::vector<int>& row, const std::vector<int>& col, int nrow, int ncol) {
    std::vector<int> all_rows, all_cols;
    all_rows.reserve(row.size()*col.size());
    all_cols.reserve(row.size()*col.size());
    for(vector<int>::const_iterator c_it=col.begin(); c_it!=col.end(); ++c_it){
      casadi_assert_message(*c_it>=0 && *c_it<ncol, "sp_rowcol: Column index out of bounds");
      for(vector<int>::const_iterator r_it=row.begin(); r_it!=row.end(); ++r_it){
        casadi_assert_message(*r_it>=0 && *r_it<nrow, "sp_rowcol: Row index out of bounds");
        all_rows.push_back(*r_it);
        all_cols.push_back(*c_it);
      }
    }
    return sp_triplet(nrow,ncol,all_rows,all_cols);
  }

  std::vector<int> getNZDense(const Sparsity &sp) {
    std::vector<int> ret(sp.size());
    std::vector<int> col = sp.getCol();
    const std::vector<int> &row = sp.row();
    int s2 = sp.size1();
    for(int k=0;k<sp.size();k++) {
      ret[k] = row[k]+col[k]*s2;
    }
    return ret;
  }

  Sparsity reshape(const Sparsity& a, int nrow, int ncol){
    return a.reshape(nrow,ncol);
  }

  Sparsity vec(const Sparsity& a){ 
    return reshape(a,a.numel(),1);
  }

  Sparsity flatten(const Sparsity& a){
    return reshape(trans(a),a.numel(),1);
  }

  Sparsity trans(const Sparsity& a) {
    return a.transpose();
  }

  Sparsity upperSparsity(const Sparsity& a, bool includeDiagonal){
    return a.upper(includeDiagonal);
  }

  Sparsity lowerSparsity(const Sparsity& a, bool includeDiagonal){
    return a.lower(includeDiagonal);
  }

  std::vector<int> upperNZ(const Sparsity& a) {
    return a.upperNZ();
  }

  std::vector<int> lowerNZ(const Sparsity& a) {
    return a.lowerNZ();
  }

  Sparsity sp_triplet(int nrow, int ncol, const std::vector<int>& row, const std::vector<int>& col, std::vector<int>& mapping, bool invert_mapping){
    // Assert dimensions
    casadi_assert_message(col.size()==row.size(),"inconsistent lengths");

    // Create the return sparsity pattern and access vectors
    Sparsity ret = Sparsity(nrow,ncol);
    vector<int> &r_colind = ret.colindRef();
    vector<int> &r_row = ret.rowRef();
    r_row.reserve(row.size());

    // Consistency check and check if elements are already perfectly ordered with no duplicates
    int last_col=-1, last_row=-1;
    bool perfectly_ordered=true;
    for(int k=0; k<col.size(); ++k){
      // Consistency check
      casadi_assert_message(col[k]>=0 && col[k]<ncol,"Col index out of bounds");
      casadi_assert_message(row[k]>=0 && row[k]<nrow,"Row index out of bounds");
    
      // Check if ordering is already perfect
      perfectly_ordered = perfectly_ordered && (col[k]<last_col || (col[k]==last_col && row[k]<=last_row));
      last_col = col[k];
      last_row = row[k];
    }
  
    // Quick return if perfectly ordered
    if(perfectly_ordered){
      // Save rows
      r_row.resize(row.size());
      copy(row.begin(),row.end(),r_row.begin());
    
      // Find offset index
      int el=0;
      for(int i=0; i<ncol; ++i){
        while(el<col.size() && col[el]==i) el++; 
        r_colind[i+1] = el;
      }
    
      // Identity mapping
      mapping.resize(row.size());
      for(int k=0; k<row.size(); ++k)
        mapping[k] = k;
    
      // Quick return
      return ret;
    }
    
    // Reuse data
    vector<int>& mapping1 = invert_mapping ? r_row : mapping;
    vector<int>& mapping2 = invert_mapping ? mapping : r_row;
  
    // Make sure that enough memory is allocated to use as a work vector
    mapping1.reserve(std::max(nrow+1,int(col.size())));
  
    // Number of elements in each row
    vector<int>& rowcount = mapping1; // reuse memory
    rowcount.resize(nrow+1);
    fill(rowcount.begin(),rowcount.end(),0);
    for(vector<int>::const_iterator it=row.begin(); it!=row.end(); ++it){
      rowcount[*it+1]++;
    }
  
    // Cumsum to get index offset for each row
    for(int i=0; i<nrow; ++i){
      rowcount[i+1] += rowcount[i];
    }
  
    // New row for each old row
    mapping2.resize(row.size());
    for(int k=0; k<row.size(); ++k){
      mapping2[rowcount[row[k]]++] = k;
    }
  
    // Number of elements in each col
    vector<int>& colcount = r_colind; // reuse memory, r_colind is already the right size and is filled with zeros
    for(vector<int>::const_iterator it=mapping2.begin(); it!=mapping2.end(); ++it){
      colcount[col[*it]+1]++;
    }
  
    // Cumsum to get index offset for each col
    for(int i=0; i<ncol; ++i){
      colcount[i+1] += colcount[i];
    }

    // New col for each old col
    mapping1.resize(col.size());
    for(vector<int>::const_iterator it=mapping2.begin(); it!=mapping2.end(); ++it){
      mapping1[colcount[col[*it]]++] = *it;
    }

    // Current element in the return matrix
    int r_el = 0;
    r_row.resize(col.size());

    // Current nonzero
    vector<int>::const_iterator it=mapping1.begin();

    // Loop over cols
    r_colind[0] = 0;
    for(int i=0; i<ncol; ++i){
    
      // Previous row (to detect duplicates)
      int j_prev = -1;
    
      // Loop over nonzero elements of the col
      while(it!=mapping1.end() && col[*it]==i){

        // Get the element
        int el = *it;
        it++;

        // Get the row
        int j = row[el];
      
        // If not a duplicate, save to return matrix
        if(j!=j_prev)
          r_row[r_el++] = j;
      
        if(invert_mapping){
          // Save to the inverse mapping
          mapping2[el] = r_el-1;        
        } else {
          // If not a duplicate, save to the mapping vector
          if(j!=j_prev)
            mapping1[r_el-1] = el;
        }
      
        // Save row
        j_prev = j;
      }
    
      // Update col offset
      r_colind[i+1] = r_el;
    }

    // Resize the row vector
    r_row.resize(r_el);

    // Resize mapping matrix
    if(!invert_mapping){
      mapping1.resize(r_el);
    }
  
    return ret;
  }


  Sparsity sp_triplet(int nrow, int ncol, const std::vector<int>& row, const std::vector<int>& col){
    std::vector<int> mapping;
    return sp_triplet(nrow,ncol,row,col,mapping,false);
  }


  Sparsity mul(const  Sparsity& a, const  Sparsity &b) {
    return (mul(DMatrix(a,1),DMatrix(b,1))).sparsity();
  }

  std::size_t hash_sparsity(int nrow, int ncol, const std::vector<int>& colind, const std::vector<int>& row){
    // Condense the sparsity pattern to a single, deterministric number
    std::size_t ret=0;
    hash_combine(ret,nrow);
    hash_combine(ret,ncol);
    hash_combine(ret,colind);
    hash_combine(ret,row);
    return ret;
  }
  
  std::vector<int> sp_compress(const Sparsity& a){
    // Get the sparsity pattern
    int nrow = a.size1();
    int ncol = a.size2();
    const vector<int>& colind = a.colind();
    const vector<int>& row = a.row();
    
    // Create compressed pattern
    vector<int> ret;
    ret.reserve(1 + 1 + colind.size() + row.size());
    ret.push_back(nrow);
    ret.push_back(ncol);
    ret.insert(ret.end(),colind.begin(),colind.end());
    ret.insert(ret.end(),row.begin(),row.end());
    return ret;
  }
  
  Sparsity sp_compress(const std::vector<int>& v){
    // Check consistency
    casadi_assert(v.size() >= 2);
    //int nrow = v[0];
    int ncol = v[1];
    casadi_assert(v.size() >= 2 + ncol+1);
    int nnz = v[2 + ncol];
    casadi_assert(v.size() == 2 + ncol+1 + nnz);

    // Call array version
    return sp_compress(&v.front());
  }
  
  Sparsity sp_compress(const int* v){
    // Get sparsity pattern
    int nrow = v[0];
    int ncol = v[1];
    const int *colind = v+2;
    int nnz = colind[ncol];
    const int *row = v + 2 + ncol+1;
    
    // Construct sparsity pattern
    return Sparsity(nrow, ncol, vector<int>(colind,colind+ncol+1),vector<int>(row,row+nnz));
  }
  
  int rank(const Sparsity& a) {
    std::vector<int> rowperm, colperm, rowblock, colblock, coarse_rowblock, coarse_colblock;
    a.dulmageMendelsohn(rowperm, colperm, rowblock, colblock, coarse_rowblock, coarse_colblock);
    return coarse_colblock.at(3);
  }

  bool isSingular(const Sparsity& a) {
    casadi_assert_message(a.size2()==a.size1(),"isSingular: only defined for square matrices, but got " << a.dimString());
    return rank(a)!=a.size2();
  }

  Sparsity sp_unit(int n, int el){
    Sparsity ret = sp_sparse(1,n);
    ret.getNZ(0,el);
    return ret;
  }
  
  Sparsity horzcat(const std::vector<Sparsity> & sp) {
    if(sp.empty()){
      return Sparsity();
    } else {
      Sparsity ret = sp[0];
      for(int i=1; i<sp.size(); ++i) {
        ret.appendColumns(sp[i]);
      }
      return ret;
    }
  }
  
  Sparsity horzcat(const Sparsity & a, const Sparsity & b) {
    Sparsity ret = a;
    ret.appendColumns(b);
    return ret;
  }

  Sparsity vertcat(const std::vector<Sparsity> & sp) {
    if(sp.empty()){
      return Sparsity();
    } else if(sp[0].vector()){
      Sparsity ret = sp[0];
      for(int i=1; i<sp.size(); ++i) {
        ret.append(sp[i]);
      }
      return ret;
    } else {
      Sparsity ret = trans(sp[0]);
      for(int i=1; i<sp.size(); ++i) {
        ret.appendColumns(trans(sp[i]));
      }
      return trans(ret);
    }
  }
  
  Sparsity vertcat(const Sparsity & a, const Sparsity & b) {
    if(a.vector()){
      Sparsity ret = a;
      ret.append(b);
      return ret;
    } else {
      Sparsity ret = trans(a);
      ret.appendColumns(trans(b));
      return trans(ret);
    }
  }
  
  Sparsity blkdiag(const std::vector< Sparsity > &v) {
    int n = 0;
    int m = 0;
    
    std::vector<int> colind(1,0);
    std::vector<int> row;
    
    int nz = 0;
    for (int i=0;i<v.size();++i) {
      const std::vector<int> &colind_ = v[i].colind();
      const std::vector<int> &row_ = v[i].row();
      for (int k=1;k<colind_.size();++k) {
        colind.push_back(colind_[k]+nz);
      }
      for (int k=0;k<row_.size();++k) {
        row.push_back(row_[k]+m);
      }
      n+= v[i].size2();
      m+= v[i].size1();
      nz+= v[i].size();
    }
    
    return Sparsity(m,n,colind,row);
  }
  
  Sparsity blkdiag(const Sparsity &a, const Sparsity &b) {
    
    std::vector<Sparsity> v;
    v.push_back(a);
    v.push_back(b);
    
    return blkdiag(v);
  }

  std::vector<Sparsity> horzsplit(const Sparsity& sp, const std::vector<int>& output_offset){
    // Number of outputs
    int n = output_offset.size()-1;

    // Get the sparsity of the input
    const vector<int>& colind_x = sp.colind();
    const vector<int>& row_x = sp.row();
    
    // Allocate result
    std::vector<Sparsity> ret;
    ret.reserve(n);

    // Sparsity pattern as CCS vectors
    vector<int> colind, row;
    int ncol, nrow = sp.size1();

    // Get the sparsity patterns of the outputs
    for(int i=0; i<n; ++i){
      int first_col = output_offset[i];
      int last_col = output_offset[i+1];
      ncol = last_col - first_col;

      // Construct the sparsity pattern
      colind.resize(ncol+1);
      copy(colind_x.begin()+first_col, colind_x.begin()+last_col+1, colind.begin());
      for(vector<int>::iterator it=colind.begin()+1; it!=colind.end(); ++it) *it -= colind[0];
      colind[0] = 0;
      row.resize(colind.back());
      copy(row_x.begin()+colind_x[first_col],row_x.begin()+colind_x[last_col],row.begin());
      
      // Append to the list
      ret.push_back(Sparsity(nrow,ncol,colind,row));
    }

    // Return (RVO)
    return ret;
  }

  std::vector<Sparsity> vertsplit(const Sparsity& sp, const std::vector<int>& output_offset){
    std::vector<Sparsity> ret = horzsplit(sp.transpose(),output_offset);
    for(std::vector<Sparsity>::iterator it=ret.begin(); it!=ret.end(); ++it){
      *it = it->transpose();
    }
    return ret;
  }

} // namespace CasADi

