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
    
std::vector<int>& CRSSparsity::colRef(){
  makeUnique();
  return (*this)->col_;
}
    
std::vector<int>& CRSSparsity::rowindRef(){
  makeUnique();
  return (*this)->rowind_;
}
    
int CRSSparsity::col(int el) const{
  return col().at(el);
}
    
int CRSSparsity::rowind(int row) const{
  return rowind().at(row);
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
  casadi_assert_message(i>=0 && i<size1() && j>=0 && j<size2(),"Indices out of bounds");

  // Quick return if matrix is dense
  if(numel()==size())
    return j+i*size2();
  
// TEST THIS  
#ifdef CRS_QUICK_APPEND
  // Quick return if we are adding an element to the end
  if(rowind(i)==size() || (rowind(i+1)==size() && col.back()<j)){
    vector<int>& colv = col();
    vector<int>& rowindv = rowind();
    colv.push_back(j);
    for(int ii=i; ii<size1(); ++ii){
      rowindv[ii+1]++;
    }
  }
#endif

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
  colRef().insert(colRef().begin()+ind,j);
  for(int row=i+1; row<size1()+1; ++row)
    rowindRef()[row]++;
  
  // Return the location of the new element
  return ind;
}

int CRSSparsity::getNZ(int i, int j) const{
  casadi_assert_message(i>=0 && i<size1() && j>=0 && j<size2(),"Indices out of bounds");
  
  // Quick return if matrix is dense
  if(numel()==size())
    return j+i*size2();

  // Find sparse element
  for(int ind=rowind(i); ind<rowind(i+1); ++ind){
    if(col(ind) == j){
      return ind;     // element exists
    }
    else if(col(ind) > j)
      break;                // break at the place where the element should be added
  }
  return -1;
}

std::vector<int> CRSSparsity::getNZNew(std::vector<int> i, std::vector<int> j){
  std::vector<int> ret;
  ret.reserve(i.size());

    // Quick return if matrix is dense
  if(numel()==size()){
    for(int k=0; k<i.size(); ++k)
      ret.push_back(j[k]+i[k]*size2());
    return ret;
  }

  // Very inefficient algorithm
  for(int k=0; k<i.size(); ++k){
    ret.push_back(getNZ(i[k],j[k]));
  }
  return ret;
}

std::vector<int> CRSSparsity::getNZNew(std::vector<int> i, std::vector<int> j) const{
  std::vector<int> ret;
  ret.reserve(i.size());

    // Quick return if matrix is dense
  if(numel()==size()){
    for(int k=0; k<i.size(); ++k)
      ret.push_back(j[k]+i[k]*size2());
    return ret;
  }

  // Very inefficient algorithm
  for(int k=0; k<i.size(); ++k){
    ret.push_back(getNZ(i[k],j[k]));
  }
  return ret;
}

std::vector<int> CRSSparsity::getNZ(std::vector<int> ii, std::vector<int> jj) const{
  std::vector<int> ret;
  for(vector<int>::const_iterator it=ii.begin(); it!=ii.end(); ++it){
    int el=rowind(*it);
    for(vector<int>::const_iterator jt=jj.begin(); jt!=jj.end(); ++jt){
      // Continue to the non-zero element
      for(; el<rowind(*it+1) && col(el)<*jt; ++el){}
      
      // Add the non-zero element, if there was an element in the location exists
      if(el<rowind(*it+1) && col(el)== *jt)
        ret.push_back(el);
      else
        ret.push_back(-1);
    }
  }
  return ret;
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
    for(int el = rowind(r); el < rowind(r+1) && col(el)<=r; ++el){
      nnz ++;
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

void CRSSparsity::bucketSort(std::vector<std::list<int> >& buckets, std::vector<int>& row) const{
  // Assert dimensions
  buckets.resize(size2());

  // Create a vector with the rows for each non-zero element
  row.resize(size());

  // Empty the buckets
  for(std::vector<std::list<int> >::iterator it=buckets.begin(); it!=buckets.end(); ++it)
    it->clear();
  
  // Quick access
  const vector<int>& c = col();
  
  // Loop over the rows of the original matrix
  for(int i=0; i<size1(); ++i)
  {
    // Loop over the elements in the row
    for(int el=rowind(i); el<rowind(i+1); ++el){ // loop over the non-zero elements
      int j=c[el];  // column
      
      // put the element into the right bucket
      buckets[j].push_back(el);
     
     // save the row index
     row[el] = i;
    }
  }
}

CRSSparsity CRSSparsity::transpose(std::vector<int>& mapping) const{
  // Non-zero entries on each column
  vector<list<int> > buckets;

  // Create a vector with the rows for each non-zero element
  vector<int> row;

  // Do a bucket sorting
  bucketSort(buckets,row);

  // create the return object
  CRSSparsity ret(size2(),size1());
  
  // Get references to the column vector and row indices
  vector<int> &rrowind = ret->rowind_;
  vector<int> &rcol = ret->col_;
  
  // reserve space (to make the calculations quicker)
  rcol.reserve(size());

  // Store the mapping of the nonzero entries
  mapping.clear();
  mapping.reserve(size());

  // loop over the rows of the resulting object
  for(int i=0; i<size2(); ++i){
    
    // Loop over the non-zero entries of the row
    for(list<int>::const_iterator it=buckets[i].begin(); it!=buckets[i].end(); ++it){
      // the index of the non-zero element
     int el =  *it;
     
     // Get the column of the element
     int j = row[el];
     
     // Store the mapping
     mapping.push_back(el);
     
     // Store the column index
     rcol.push_back(j);
    }
    // Log the row index
    rrowind[i+1] = rcol.size();
  }
  
  // Return the sparsity
  return ret;

}

void CRSSparsity::reserve(int nnz, int nrow){
  colRef().reserve(nnz);
  rowindRef().reserve(nrow+1);
}

bool CRSSparsity::operator==(const CRSSparsity& y) const{
  // Quick true if the objects are the same
  if(get() == y.get())
    return true;
  
  // Quick false if dimensions or number of non-zeros are different
  if(size1() != y.size1() || size2() != y.size2() || size() != y.size())
    return false;
  
  // Check rowind vectors
  const std::vector<int>& r1 = rowind();
  const std::vector<int>& r2 = y.rowind();
  for(int i=0; i<r1.size(); ++i)
    if(r1[i]!=r2[i])
      return false;
  
  // Check column vector
  const std::vector<int>& c1 = col();
  const std::vector<int>& c2 = y.col();
  for(int i=0; i<c1.size(); ++i)
    if(c1[i]!=c2[i])
      return false;

  // Patterns are identical if we reached this point
  return true;
}

void CRSSparsity::append(const CRSSparsity& sp){
  // Assert dimensions
  casadi_assert_message(size2()==sp.size2(),"Dimension mismatch");
  
  // Get current sparsity pattern
  vector<int>& col_ = colRef();
  vector<int>& rowind_ = rowindRef();

  // Get current number of non-zeros
  int sz = size();
  
  // Add column indices
  col_.insert(col_.end(),sp.col().begin(),sp.col().end());
  
  // Add row indices
  rowind_.pop_back();
  rowind_.insert(rowind_.end(),sp.rowind().begin(),sp.rowind().end());
  for(int i = size1(); i<rowind_.size(); ++i)
    rowind_[i] += sz;
  
  // Update dimensions
  (*this)->nrow_ += sp.size1();
}


} // namespace CasADi


