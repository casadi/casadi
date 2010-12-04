/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A minimalistic computer algebra system with automatic differentiation 
 *    and framework for dynamic optimization.
 *    Copyright (C) 2010 by Joel Andersson, Moritz Diehl et al., K.U.Leuven. All rights reserved.
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
}

void FunctionIO::init(){
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
  data_.resize(col_.size());

  // forward derivatives
  for(int i=0; i<dataF_.size(); ++i)
    dataF_[i].resize(col_.size());

  // adjoint derivatives
  for(int i=0; i<dataA_.size(); ++i)
    dataA_[i].resize(col_.size());
}

void FunctionIO::setSize(int nrow, int ncol){
  nrow_ = nrow;
  ncol_ = ncol;
}

void FunctionIO::setSparsityCRS(const vector<int>& rowind, const vector<int> &col){
  rowind_ = rowind;
  col_ = col;
}

void FunctionIO::setNumFwdDir(int ndir){
  dataF_.resize(ndir);
}

void FunctionIO::setNumAdjDir(int ndir){
  dataA_.resize(ndir);
}

int FunctionIO::numFwdDir() const{
  return dataF_.size();
}

int FunctionIO::numAdjDir() const{
  return dataA_.size();
}

vector<double>& FunctionIO::dataF(int dir){
  return dataF_.at(dir);
}

const vector<double>& FunctionIO::dataF(int dir) const{
  return dataF_.at(dir);
}

vector<double>& FunctionIO::dataA(int dir){
  return dataA_.at(dir);
}

const vector<double>& FunctionIO::dataA(int dir) const{
  return dataA_.at(dir);
}

vector<double>& FunctionIO::data(){
  return data_;
}

const vector<double>& FunctionIO::data() const{
  return data_;
}


vector<double>& FunctionIO::data(int ord){
  assert(ord==0 || ord==1);
  if(ord==0)
    return data_;
  else 
    return dataF_.at(0);
}

const vector<double>& FunctionIO::data(int ord) const{
  assert(ord==0 || ord==1);
  if(ord==0)
    return data_;
  else 
    return dataF_.at(0);
}







// MatrixSize FunctionIO::size() const{
//   return MatrixSize(nrow_,ncol_);
// }

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

void FunctionIO::getSparseSym(double *res, int ord) const{
  // copy to the result vector
  int nz = 0;
  for(int row=0; row<size1(); ++row)
  {
    // Loop over the elements in the row
    for(int el=rowind_[row]; el<rowind_[row+1]; ++el){ // loop over the non-zero elements
      if(col_[el] > row) break; // break inner loop (only lower triangular part is used)
      res[nz] = data(ord)[el];
      nz++;
    }
  }
}

void FunctionIO::getTimesVector(const double *v, double *res, int ord) const{
  // copy the result
  for(int i=0; i<size1(); ++i){ // loop over rows
    res[i] = 0;
    for(int el=rowind_[i]; el<rowind_[i+1]; ++el){ // loop over the non-zero elements
      int j=col_[el];  // column

      // Multiply with the vector
      res[i] += v[j]*data(ord)[el];
    }
  }
}

void FunctionIO::getBand(int kl, int ku, int ldres, double *res, int ord) const{
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
      res[s + ldres*j] = data(ord)[el];
    }
  }
}

void FunctionIO::getDense(double *x, int ord) const{
  const vector<double>& arg = data(ord);

  
  if(col_.size() == numel()){
    // if dense matrix
    for(int i=0; i<arg.size(); ++i)
      x[i] = arg[i];
  } else {
    // general sparse
    int k=0; // index of the result
    for(int i=0; i<size1(); ++i) // loop over rows
      for(int el=rowind_[i]; el<rowind_[i+1]; ++el){ // loop over the non-zero elements
        int j=col_[el];  // column
        for(; k<i*size2()+j; ++k) 
          x[k] = 0; // add zeros before the non-zero element

      // add the non-zero element
      x[k] = data(ord)[el];
      k++;
    }
    // add sparse zeros at the end of the matrix
    for(; k<numel(); ++k)
     x[k] = 0;
   
    assert(k==numel());
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





} // namespace CasADi

