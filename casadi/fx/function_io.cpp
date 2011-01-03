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

#include "function_io.hpp"
#include <sstream>
#include "../stl_vector_tools.hpp"
#include "../casadi_exception.hpp"

using namespace std;

namespace CasADi{

FunctionIO::FunctionIO(){
  setSize(0);
  nfdir_ = 0;
  nadir_ = 0;
}

void FunctionIO::init(){
  data_.resize(nadir_+1+nfdir_);
  
  // Dense matrix if no sparsity defined
  if(rowind_.empty()){
    rowind_.resize(nrow_+1);
    col_.resize(nrow_*ncol_);
    int el=0;
    for(int i=0; i<nrow_; ++i){
      rowind_[i] = el;
      for(int j=0; j<ncol_; ++j)
        col_[el++] = j;
    }
    rowind_[nrow_] = nrow_*ncol_;
  }

  // Non-zeros
  for(vector<vector<double> >::iterator it=data_.begin(); it!=data_.end(); ++it)
    it->resize(col_.size());
}

void FunctionIO::setSize(int nrow, int ncol){
  nrow_ = nrow;
  ncol_ = ncol;
}

void FunctionIO::setSparsityCRS(const vector<int>& rowind, const vector<int> &col){
  rowind_ = rowind;
  col_ = col;
}

void FunctionIO::setNumFwdDir(int nfdir){
  nfdir_ = nfdir;
}

void FunctionIO::setNumAdjDir(int nadir){
  nadir_ = nadir;
}

int FunctionIO::numFwdDir() const{
  return nfdir_;
}

int FunctionIO::numAdjDir() const{
  return nadir_;
}

vector<double>& FunctionIO::dataF(int dir){
  return data(1+dir);
}

const vector<double>& FunctionIO::dataF(int dir) const{
  return data(1+dir);
}

vector<double>& FunctionIO::dataA(int dir){
  return data(-1-dir);
}

const vector<double>& FunctionIO::dataA(int dir) const{
  return data(-1-dir);
}

int FunctionIO::numel() const{
  return nrow_*ncol_;
}

int FunctionIO::size1() const{
  return nrow_;
}

int FunctionIO::size2() const{
  return ncol_;
}


void FunctionIO::setSparsity(const vector<int>& rr, const vector<int>& cc, vector<int>& elind){
  // Number of non-zeros
  int nnz = cc.size();

  // Row indices
  rowind_.resize(size1()+1);
  
  // Columns
  col_ = cc;
  
  // Transformation
  elind.resize(nnz);
  
  rowind_[0]=0;
  int el=0;
  for(int r=0; r<size1(); ++r){
    while(el < nnz && rr.at(el)==r){
      // Transformation
      elind.at(el) = el;
      
      el++;
      if(el<nnz)
	assert(rr.at(el)>=r);
    }
    rowind_.at(r+1) = el;
  }
  assert(el==nnz);
}

void FunctionIO::getSparseSym(double *res, int dir) const{
  // copy to the result vector
  int nz = 0;
  for(int row=0; row<size1(); ++row)
  {
    // Loop over the elements in the row
    for(int el=rowind_[row]; el<rowind_[row+1]; ++el){ // loop over the non-zero elements
      if(col_[el] > row) break; // break inner loop (only lower triangular part is used)
      res[nz] = data(dir)[el];
      nz++;
    }
  }
}

void FunctionIO::getTimesVector(const double *v, double *res, int dir) const{
  // copy the result
  for(int i=0; i<size1(); ++i){ // loop over rows
    res[i] = 0;
    for(int el=rowind_[i]; el<rowind_[i+1]; ++el){ // loop over the non-zero elements
      int j=col_[el];  // column

      // Multiply with the vector
      res[i] += v[j]*data(dir)[el];
    }
  }
}

void FunctionIO::getBand(int kl, int ku, int ldres, double *res, int dir) const{
  // delete the content of the matrix
  for(int j=0; j<size2(); ++j) // loop over columns
    for(int s=0; s<kl+ku+1; ++s) // loop over the subdiagonals
      res[s + ldres*j] = 0;
  
  // loop over rows
  for(int i=0; i<size1(); ++i){ 
    
    // loop over the non-zero elements
    for(int el=rowind_[i]; el<rowind_[i+1]; ++el){ 
      int j=col_[el];  // column
      
      // Check if we have not yet inside the band
      if(j<i-kl) continue;

      // Check if we are already outside the band
      if(j>i+ku) break;

      // Get the subdiagonal
      int s = i - j + ku;

      // Store the element
      res[s + ldres*j] = data(dir)[el];
    }
  }
}


void FunctionIO::getSparsityCRS(vector<int>& rowind, vector<int> &col) const{
    rowind = rowind_;
    col = col_;
}

void FunctionIO::getSparsity(vector<int>& row, vector<int> &col) const{
  int nnz = col_.size();
  row.resize(nnz);
  col = col_;
  for(int r=0; r<size1(); ++r){
    for(int el = rowind_[r]; el < rowind_[r+1]; ++el){
        row[el] = r;
      }
  }
}

void FunctionIO::assertNNZ(int sz, Sparsity sp) const{
  int nnz_correct = -1;
  switch(sp){
    case SPARSE:
      nnz_correct = size();
      break;
    case DENSE:
      nnz_correct = numel();
      break;
    default:
      throw CasadiException("FunctionIO::assertNNZ: unknown sparsity");
  }

  if(nnz_correct!=sz){
    stringstream ss;
    ss << "FunctionIO::assertNNZ: wrong number of elements (" << sz << "), but should be " << nnz_correct << flush;
    throw CasadiException(ss.str());
  }
}

void FunctionIO::assertNumEl(int sz) const{
  if(numel()!=sz){
    stringstream ss;
    ss << "FunctionIO::assertNumEl: wrong number of elements (" << sz << " ), but should be " << numel();
    throw CasadiException(ss.str());
  }
}

int FunctionIO::size() const{
  return col_.size();
}

int FunctionIO::sizeU() const{
  int nnz = 0;
  for(int r=0; r<size1(); ++r){
    for(int el = rowind_[r]; el < rowind_[r+1]; ++el){
      nnz += col_[el]>=r;
    }
  }
  return nnz;
}

int FunctionIO::sizeL() const{
  int nnz = 0;
  for(int r=0; r<size1(); ++r){
    for(int el = rowind_[r]; el < rowind_[r+1]; ++el){
      nnz += col_[el]<=r;
    }
  }
  return nnz;
}

void FunctionIO::setv(double val, vector<double>& v, Sparsity sp) const{
  assertNumEl(1);
  assertNNZ(1,sp);
  v[0] = val;
}
    
void FunctionIO::getv(double& val, const vector<double>& v, Sparsity sp) const{
  assertNumEl(1);
  assertNNZ(1,sp);
  val = v[0];
}

void FunctionIO::setv(const vector<double>& val, vector<double>& v, Sparsity sp) const{
  assertNNZ(val.size(),sp);
  copy(val.begin(),val.end(),v.begin());  
}

void FunctionIO::getv(vector<double>& val, const vector<double>& v, Sparsity sp) const{
  if(sp==SPARSE || (sp==DENSE && numel()==size())){
    assertNNZ(val.size(),sp);
    copy(v.begin(),v.end(),val.begin());
  } else {
    assertNNZ(val.size(),sp);
    // general sparse
    int k=0; // index of the result
    for(int i=0; i<size1(); ++i) // loop over rows
      for(int el=rowind_[i]; el<rowind_[i+1]; ++el){ // loop over the non-zero elements
        int j=col_[el];  // column
        for(; k<i*size2()+j; ++k)
          val[k] = 0; // add zeros before the non-zero element

      // add the non-zero element
      val[k] = v[el];
      k++;
    }
    // add sparse zeros at the end of the matrix
    for(; k<numel(); ++k)
     val[k] = 0;
  }
}

void FunctionIO::setv(const double* val, vector<double>& v, Sparsity sp) const{
  copy(val,val+v.size(),v.begin());  
}

void FunctionIO::getv(double* val, const vector<double>& v, Sparsity sp) const{
  copy(v.begin(),v.end(),val);
}

vector<double>& FunctionIO::data(int dir){
  return data_.at(nadir_+dir);
}

const vector<double>& FunctionIO::data(int dir) const{
  return data_.at(nadir_+dir);
}

void FunctionIO::set(double val, int dir, Sparsity sp){
  setv(val,data(dir),sp);
}
    
void FunctionIO::get(double& val, int dir, Sparsity sp) const{
  getv(val,data(dir),sp);
}

void FunctionIO::set(const std::vector<double>& val, int dir, Sparsity sp){
  setv(val,data(dir),sp);
}

void FunctionIO::get(std::vector<double>& val, int dir, Sparsity sp) const{
  getv(val,data(dir),sp);
}

void FunctionIO::set(const double* val, int dir, Sparsity sp){
  setv(val,data(dir),sp);
}

void FunctionIO::get(double* val, int dir, Sparsity sp) const{
  getv(val,data(dir),sp);
}




} // namespace CasADi

