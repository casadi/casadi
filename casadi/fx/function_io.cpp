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
  data_.resize(1);
  dense_ = true;
  nfdir_ = 0;
  nadir_ = 0;
}

void FunctionIO::init(){
  
  // Make dense if necessary
  if(dense_){
    data_[0] = Matrix<double>(size1(),size2(),0);
  }

  // Non-zeros
  data_.resize(nadir_+1+nfdir_);
  for(int i=0; i<nfdir_+nadir_; ++i)
    data_[i+1] = data_[i];
}

void FunctionIO::setSize(int nrow, int ncol){
  data_[0].resize(nrow,ncol);
}

void FunctionIO::setSparsityCRS(const vector<int>& rowind, const vector<int> &col){
  data_[0] = Matrix<double>(size1(),size2(),col,rowind);
  dense_ = false;
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
  return size1()*size2();
}

int FunctionIO::size1() const{
  return data_[0].size1();
}

int FunctionIO::size2() const{
  return data_[0].size2();
}


// void FunctionIO::setSparsity(const vector<int>& rr, const vector<int>& cc, vector<int>& elind){
//   // Number of non-zeros
//   int nnz = cc.size();
// 
//   // Row indices
//   rowind().resize(size1()+1);
//   
//   // Columns
//   col_ = cc;
//   
//   // Transformation
//   elind.resize(nnz);
//   
//   rowind()[0]=0;
//   int el=0;
//   for(int r=0; r<size1(); ++r){
//     while(el < nnz && rr.at(el)==r){
//       // Transformation
//       elind.at(el) = el;
//       
//       el++;
//       if(el<nnz)
// 	assert(rr.at(el)>=r);
//     }
//     rowind().at(r+1) = el;
//   }
//   assert(el==nnz);
//   
//   dense_ = false;
// 
// }

void FunctionIO::getSparseSym(double *res, int dir) const{
  // copy to the result vector
  int nz = 0;
  for(int row=0; row<size1(); ++row)
  {
    // Loop over the elements in the row
    for(int el=rowind(row); el<rowind(row+1); ++el){ // loop over the non-zero elements
      if(col(el) > row) break; // break inner loop (only lower triangular part is used)
      res[nz] = data(dir)[el];
      nz++;
    }
  }
}

void FunctionIO::getTimesVector(const double *v, double *res, int dir) const{
  // copy the result
  for(int i=0; i<size1(); ++i){ // loop over rows
    res[i] = 0;
    for(int el=rowind(i); el<rowind(i+1); ++el){ // loop over the non-zero elements
      int j=col(el);  // column

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
    for(int el=rowind(i); el<rowind(i+1); ++el){ 
      int j=col(el);  // column
      
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
  rowind = data_[0].rowind();
  col = data_[0].col();
}

void FunctionIO::getSparsity(vector<int>& row, vector<int> &col) const{
  int nnz = size();
  row.resize(nnz);
  col = this->col();
  for(int r=0; r<size1(); ++r){
    for(int el = rowind(r); el < rowind(r+1); ++el){
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
  return col().size();
}

int FunctionIO::sizeU() const{
  int nnz = 0;
  for(int r=0; r<size1(); ++r){
    for(int el = rowind(r); el < rowind(r+1); ++el){
      nnz += col(el)>=r;
    }
  }
  return nnz;
}

int FunctionIO::sizeL() const{
  int nnz = 0;
  for(int r=0; r<size1(); ++r){
    for(int el = rowind(r); el < rowind(r+1); ++el){
      nnz += col(el)<=r;
    }
  }
  return nnz;
}

vector<double>& FunctionIO::data(int dir){
  return data_.at(nadir_+dir);
}

const vector<double>& FunctionIO::data(int dir) const{
  return data_.at(nadir_+dir);
}

void FunctionIO::set(double val, int dir, Sparsity sp){
  assertNumEl(1);
  assertNNZ(1,sp);
  data(dir)[0] = val;
}
    
void FunctionIO::get(double& val, int dir, Sparsity sp) const{
  assertNumEl(1);
  assertNNZ(1,sp);
  val = data(dir)[0];
}

void FunctionIO::set(const std::vector<double>& val, int dir, Sparsity sp){
  assertNNZ(val.size(),sp);
  set(&val[0],dir,sp);
}

void FunctionIO::get(std::vector<double>& val, int dir, Sparsity sp) const{
  assertNNZ(val.size(),sp);
  get(&val[0],dir,sp);
}

void FunctionIO::set(const double* val, int dir, Sparsity sp){
  vector<double> &v = data(dir);
  if(sp==SPARSE || (sp==DENSE && numel()==size())){
    copy(val,val+v.size(),v.begin());
  } else if(sp==DENSE){
    for(int i=0; i<size1(); ++i) // loop over rows
      for(int el=rowind(i); el<rowind(i+1); ++el){ // loop over the non-zero elements
        // column
        int j=col(el);
        
        // Set the element
        v[el] = val[i*size2()+j];
    }
  } else {
    throw CasadiException("FunctionIO::set: not SPARSE or DENSE");
  }
}

void FunctionIO::get(double* val, int dir, Sparsity sp) const{
  const vector<double> &v = data(dir);
  if(sp==SPARSE || (sp==DENSE && numel()==size())){
    copy(v.begin(),v.end(),val);
  } else if(sp==DENSE){
    int k=0; // index of the result
    for(int i=0; i<size1(); ++i) // loop over rows
      for(int el=rowind(i); el<rowind(i+1); ++el){ // loop over the non-zero elements
        int j=col(el);  // column
        for(; k<i*size2()+j; ++k)
          val[k] = 0; // add zeros before the non-zero element

      // add the non-zero element
      val[k] = v[el];
      k++;
    }
    // add sparse zeros at the end of the matrix
    for(; k<numel(); ++k)
     val[k] = 0;
  } else {
    throw CasadiException("FunctionIO::get: not SPARSE or DENSE");
  }
}


const std::vector<int>& FunctionIO::rowind() const{
  return data_[0].rowind();
}

std::vector<int>& FunctionIO::rowind(){
  return data_[0].rowind();
}

const std::vector<int>& FunctionIO::col() const{
  return data_[0].col();
}

std::vector<int>& FunctionIO::col(){
  return data_[0].col();
}

int FunctionIO::rowind(int i) const{
  return data_[0].rowind(i);
}

int FunctionIO::col(int el) const{
  return data_[0].col(el);
}


} // namespace CasADi

