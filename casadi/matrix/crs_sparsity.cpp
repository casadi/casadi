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

#include "crs_sparsity.hpp"
#include "../stl_vector_tools.hpp"

using namespace std;

namespace CasADi{

CRSSparsity::CRSSparsity(){
}

CRSSparsity::CRSSparsity(int nrow, int ncol, bool dense){
  vector<int> col, rowind(nrow+1,0);
  if(dense){
    col.resize(nrow*ncol);
    rowind.resize(nrow+1);
    for(int i=0; i<nrow+1; ++i)
      rowind[i] = i*ncol;
    for(int i=0; i<nrow; ++i)
      for(int j=0; j<ncol; ++j)
        col[j+i*ncol] = j;
  }
  assignNode(new CRSSparsityNode(nrow, ncol, col, rowind));
}

CRSSparsity::CRSSparsity(int nrow, int ncol, std::vector<int> col, std::vector<int> rowind){
  assignNode(new CRSSparsityNode(nrow, ncol, col, rowind));
}
    
CRSSparsityNode* CRSSparsity::operator->(){
  makeUnique();
  return (CRSSparsityNode*)(SharedObject::operator->());
}

const CRSSparsityNode* CRSSparsity::operator->() const{
  return (const CRSSparsityNode*)(SharedObject::operator->());
}
  
bool CRSSparsity::checkNode() const{
  return dynamic_cast<const CRSSparsityNode*>(get());
}

int CRSSparsity::size1() const{
  return (*this)->nrow_;
}
    
int CRSSparsity::size2() const{
  return (*this)->ncol_;
}
    
int CRSSparsity::numel() const{
  return size1()*size2();
}
    
int CRSSparsity::size() const{
  return (*this)->col_.size();
}
    
const std::vector<int>& CRSSparsity::col() const{
  return (*this)->col_;
}
    
const std::vector<int>& CRSSparsity::rowind() const{
  return (*this)->rowind_;
}
    
std::vector<int>& CRSSparsity::col(){
  makeUnique();
  return (*this)->col_;
}
    
std::vector<int>& CRSSparsity::rowind(){
  makeUnique();
  return (*this)->rowind_;
}
    
int CRSSparsity::col(int el) const{
  return col()[el];
}
    
int CRSSparsity::rowind(int row) const{
  return rowind()[row];
}

void CRSSparsity::resize(int nrow, int ncol){
  if(nrow != size1() || ncol != size2()){
    if(nrow < size1() || ncol < size2()){
      // Row and column index of the new
      vector<int> col_new, rowind_new(nrow+1,0);

      // Loop over the rows which may contain nonzeros
      int i;
      for(i=0; i<size1() && i<nrow; ++i){
        // First nonzero element of the row
        rowind_new[i] = col_new.size();
        
        // Record columns of the nonzeros
        for(int el=rowind(i); el<rowind(i+1) && col(el)<ncol; ++el){
          col_new.push_back(col(el));
        }
      }
      
      // Save row-indices for the rest of the rows
      for(; i<nrow+1; ++i){
        rowind_new[i] = col_new.size();
      }
        
      // Save the sparsity
      *this = CRSSparsity(nrow,ncol,col_new,rowind_new);
      
    } else {
      // Make larger: Very cheap operation
      (*this)->nrow_ = nrow;
      (*this)->ncol_ = ncol;
      (*this)->rowind_.resize(size1()+1,size());
    }
  }
}

int CRSSparsity::getNZ(int i, int j){
  if(i >= size1() || j>=size2()) throw CasadiException("CRSSparsity::getNZ: out of bounds");

  // Quick return if matrix is dense
  if(numel()==size())
    return j+i*size2();
  
  // go to the place where the element should be
  int ind;
  for(ind=rowind(i); ind<rowind(i+1); ++ind){ // better: loop from the back to the front
    if(col(ind) == j){
      return ind; // element exists
    } else if(col(ind) > j)
      break;                // break at the place where the element should be added
  }
  
  // Make sure that there no other objects are affected
  makeUnique();
  
  // insert the element
  col().insert(col().begin()+ind,j);
  for(int row=i+1; row<size1()+1; ++row)
    rowind()[row]++;
  
  // Return the location of the new element
  return ind;
}

int CRSSparsity::getNZ(int i, int j) const{
  if(i >= size1() || j>=size2()) throw CasadiException("CRSSparsity::getNZ: out of bounds");
  
  // Quick return if matrix is dense
  if(numel()==size())
    return j+i*size2();

  // Find sparse element
  for(int ind=rowind(i); ind<rowind(i+1); ++ind){
    if(col(ind) == j)
      return ind;     // element exists
    else if(col(ind) > j)
      break;                // break at the place where the element should be added
  }
  return -1;
}

int CRSSparsity::sizeU() const{
  int nnz = 0;
  for(int r=0; r<size1(); ++r){
    for(int el = rowind(r); el < rowind(r+1); ++el){
      nnz += col(el)>=r;
    }
  }
  return nnz;
}

int CRSSparsity::sizeL() const{
  int nnz = 0;
  for(int r=0; r<size1(); ++r){
    for(int el = rowind(r); el < rowind(r+1); ++el){
      nnz += col(el)<=r;
    }
  }
  return nnz;
}


void CRSSparsityNode::repr(std::ostream &stream) const{
  stream << "Compressed Row Storage: " << nrow_ << "-by-" << ncol_ << " matrix, " << col_.size() << " structural non-zeros";
}

void CRSSparsityNode::print(std::ostream &stream) const{
  repr(stream);
  stream << endl;
  stream << "col:    " << col_ << endl;
  stream << "rowind: " << rowind_ << endl;
}

vector<int> CRSSparsity::getRow() const{
  vector<int> row(size());
  for(int r=0; r<size1(); ++r){
    for(int el = rowind(r); el < rowind(r+1); ++el){
        row[el] = r;
      }
  }
  return row;
}

void CRSSparsity::getSparsityCRS(std::vector<int>& rowind, std::vector<int> &col) const{
  rowind = this->rowind();
  col = this->col();
}

void CRSSparsity::getSparsity(std::vector<int>& row, std::vector<int> &col) const{
  row = this->getRow();
  col = this->col();
}

void CRSSparsity::bucketSort(std::vector<std::vector<int> >& buckets, std::vector<int>& row) const{
  // Assert dimensions
  buckets.resize(size2());

  // Create a vector with the rows for each non-zero element
  row.resize(size());

  // Empty the buckets
  for(std::vector<std::vector<int> >::iterator it=buckets.begin(); it!=buckets.end(); ++it)
    it->clear();
  
  // Loop over the rows of the original matrix
  for(int i=0; i<size1(); ++i)
  {
    // Loop over the elements in the row
    for(int el=rowind(i); el<rowind(i+1); ++el){ // loop over the non-zero elements
      int j=col(el);  // column
      
     // put the element into the right bucket
     buckets[j].push_back(el);
     
     // save the row index
     row[el] = i;
    }
  }
}

CRSSparsity CRSSparsity::transpose(std::vector<int>& mapping) const{
  // Non-zero entries on each column
  std::vector<std::vector<int> > buckets;

  // Create a vector with the rows for each non-zero element
  std::vector<int> row;

  // Do a bucket sorting
  bucketSort(buckets,row);

  // create the return object
  CRSSparsity ret(size2(),size1());

  // reserve space (to make the calculations quicker)
  ret.reserve(size(),size2());

  // Store the mapping of the nonzero entries
  mapping.resize(size());

  // loop over the columns of the object to be transposed
  for(int j=0; j<size2(); ++j){
    
    // Loop over the non-zero entries of the column
    for(int r=0; r<buckets[j].size(); ++r){
      // the index of the non-zero element
     int el =  buckets[j][r];
     
     // Get the row of the element
     int i = row[el]; 
     
     // add the element (can be done much more efficiently!)
     int el_ret = ret.getNZ(j,i);
     
     // Store the mapping
     mapping[el_ret] = el;
    }
  }
  
  // Return the sparsity
  return ret;

}

void CRSSparsity::reserve(int nnz, int nrow){
  col().reserve(nnz);
  rowind().reserve(nrow+1);
}


} // namespace CasADi


