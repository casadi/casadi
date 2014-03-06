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

  Sparsity mul(const  Sparsity& a, const  Sparsity &b) {
    return (mul(DMatrix(a,1),DMatrix(b,1))).sparsity();
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

