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

CRSSparsity::CRSSparsity(int nrow, int ncol, vector<int> col, vector<int> rowind){
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
    
const vector<int>& CRSSparsity::col() const{
  return (*this)->col_;
}
    
const vector<int>& CRSSparsity::rowind() const{
  return (*this)->rowind_;
}
    
vector<int>& CRSSparsity::colRef(){
  makeUnique();
  return (*this)->col_;
}
    
vector<int>& CRSSparsity::rowindRef(){
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
  
  // Quick return if we are adding an element to the end
  if(rowind(i)==size() || (rowind(i+1)==size() && col().back()<j)){
    vector<int>& colv = colRef();
    vector<int>& rowindv = rowindRef();
    colv.push_back(j);
    for(int ii=i; ii<size1(); ++ii){
      rowindv[ii+1]++;
    }
    return colv.size()-1;
  }

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
  
  // Quick return if past the end
  if(rowind(i)==size() || (rowind(i+1)==size() && col().back()<j)){
    return -1;
  }

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

vector<int> CRSSparsity::getNZNew(vector<int> i, vector<int> j){
  vector<int> ret;
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

vector<int> CRSSparsity::getNZNew(vector<int> i, vector<int> j) const{
  vector<int> ret;
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

vector<int> CRSSparsity::getNZ(vector<int> ii, vector<int> jj) const{
  vector<int> ret;
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

bool CRSSparsity::dense() const{
  return size() == numel();
}

CRSSparsity CRSSparsity::getSub(const std::vector<int>& ii, const std::vector<int>& jj, std::vector<int>& mapping) const{
  if(dense()){
    CRSSparsity ret(ii.size(),jj.size(),true);
    mapping.resize(ret.size());
    for(int i=0; i!=ii.size(); ++i){
      for(int j=0; j!=jj.size(); ++j){
        mapping[j+i*ret.size2()] = jj[j] + ii[i]*size2();
      }
    }
    return ret;
  } else {
    CRSSparsity ret(ii.size(),jj.size(),false);
    mapping.resize(0);
    for(int i=0; i!=ii.size(); ++i){
      int j=0;
      for(int el=rowind(ii[i]); el<rowind(ii[i]+1); ++el){
        while(jj[j]<col(el) && j<jj.size()){
          ++j;
        }
          
        if(j<jj.size() && jj[j]==col(el)){
          int ind = ret.getNZ(i,j);
          casadi_assert(ind==ret.size()-1); // make sure that we are adding to the end
          mapping.push_back(el);
        }
      }
    }
    return ret;
  }
}

CRSSparsity CRSSparsity::eraseSub(const std::vector<int>& ii, const std::vector<int>& jj, std::vector<int>& mapping) const{
  CRSSparsity ret(size1(),size2(),false);
  mapping.resize(0);
  for(int i=0; i!=ii.size(); ++i){
    int j=0;
    for(int el=rowind(ii[i]); el<rowind(ii[i]+1); ++el){
      while(j<jj.size() && jj[j]<col(el)){
        ++j;
      }
        
      if(!(j<jj.size() && jj[j]==col(el))){
        int ind = ret.getNZ(ii[i],col(el));
        casadi_assert(ind==ret.size()-1); // make sure that we are adding to the end
        mapping.push_back(el);
      }
    }
  }
  return ret;
}

void CRSSparsity::setSub(const CRSSparsity& sub, const std::vector<int>& ii, const std::vector<int>& jj, std::vector<int>& mapping_nz, std::vector<int>& mapping_ind){
  if(dense() && sub.dense()){
    mapping_nz.resize(size());
    for(int k=0; k<size(); ++k) mapping_nz[k] = k;
    mapping_ind.resize(sub.size());
    fill(mapping_ind.begin(),mapping_ind.end(),0);
    for(int i=0; i!=ii.size(); ++i){
      for(int j=0; j!=jj.size(); ++j){
        mapping_nz[jj[j] + ii[i]*size2()] = 1;
        mapping_ind[jj[j] + ii[i]*size2()] = j+i*sub.size2();
      }
    }
  } else {
    // Find out the non-zeros to be erased
    vector<int> to_be_erased;
    getSub(ii,jj,to_be_erased);
    
    // Get the sparsity vectors
    vector<int> &c = colRef();
    vector<int> &rind = rowindRef();
    
    casadi_assert_message(0, "Not implemented");
  }
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


void CRSSparsityNode::repr(ostream &stream) const{
  stream << "Compressed Row Storage: " << nrow_ << "-by-" << ncol_ << " matrix, " << col_.size() << " structural non-zeros";
}

void CRSSparsityNode::print(ostream &stream) const{
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

void CRSSparsity::getSparsityCRS(vector<int>& rowind, vector<int> &col) const{
  rowind = this->rowind();
  col = this->col();
}

void CRSSparsity::getSparsity(vector<int>& row, vector<int> &col) const{
  row = this->getRow();
  col = this->col();
}

void CRSSparsity::bucketSort(vector<list<int> >& buckets, vector<int>& row) const{
  // Assert dimensions
  buckets.resize(size2());

  // Create a vector with the rows for each non-zero element
  row.resize(size());

  // Empty the buckets
  for(vector<list<int> >::iterator it=buckets.begin(); it!=buckets.end(); ++it)
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

CRSSparsity CRSSparsity::transpose(vector<int>& mapping) const{
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

CRSSparsity CRSSparsity::patternUnion(const CRSSparsity& y, vector<int>& mapping) const{
  // Assert dimensions
  casadi_assert_message(size1()==y.size1(), "The number of rows does not match");
  casadi_assert_message(size2()==y.size2(), "The number of columns does not match");
  
  // Quick return if the patterns are equal
  if(*this == y){
    mapping.resize(size());
    fill(mapping.begin(),mapping.end(),0);
    return *this;
  }
  
  // Create return object
  CRSSparsity ret(size1(),size2());
  
  // Get refences to the sparsity vectors
  vector<int>& r = ret.rowindRef();
  vector<int>& c = ret.colRef();
  
  // Prepare the assembly of the rowind vector below
  r.clear();
  r.push_back(0);
  
  // Clear the mapping
  mapping.clear();
  
  // Loop over rows of both patterns
  for(int i=0; i<size1(); ++i){
    // Non-zero element of the two matrices
    int el1 = rowind(i);
    int el2 = y.rowind(i);
    
    // End of the non-zero elements of the row for the two matrices
    int el1_last = rowind(i+1);
    int el2_last = y.rowind(i+1);
    
    // Loop over the non-zeros of both matrices
    while(el1<el1_last || el2<el2_last){
      // Get the columns
      int col1 = col(el1);
      int col2 = col(el2);
      
      // Add to the return matrix
      if(col1==col2){
        c.push_back(col1);
        mapping.push_back(0);
        el1++; el2++;
      } else if(col1<col2){
        c.push_back(col1);
        mapping.push_back(-1);
        el1++;
      } else {
        c.push_back(col2);
        mapping.push_back(1);
        el2++;
      }
    }
    
    // Save the index of the last nonzero on the row
    r.push_back(c.size());
  }
  
  // Make sure that the object was correctly created
  casadi_assert(r.size()==size1()+1);
  casadi_assert(mapping.size()==c.size());
  casadi_assert(c.size()==r.back());
  
  // Return 
  return ret;
}

CRSSparsity CRSSparsity::patternIntersection(const CRSSparsity& y, vector<int>& mapping) const{
  throw CasadiException("CRSSparsity::patternIntersection not implemented");
}

CRSSparsity CRSSparsity::patternProduct(const CRSSparsity& y_trans) const{
  vector< vector< pair<int,int> > > dummy; // dummy argument
  return patternProduct(y_trans, dummy, false);
}

CRSSparsity CRSSparsity::patternProduct(const CRSSparsity& y_trans, vector< vector< pair<int,int> > >& mapping, bool with_mapping) const{
  // return object
  CRSSparsity ret(size1(),y_trans.size1());
  
  // Get the vectors for the return pattern
  vector<int>& c = ret.colRef();
  vector<int>& r = ret.rowindRef();
  
  // Direct access to the arrays
  const vector<int> &x_col = col();
  const vector<int> &y_row = y_trans.col();
  const vector<int> &x_rowind = rowind();
  const vector<int> &y_colind = y_trans.rowind();

  if(with_mapping){
  
    // Clear the mapping
    mapping.clear();

    // the entry of the matrix to be calculated
    vector< pair<int,int> > d; 

    // loop over the row of the resulting matrix)
    for(int i=0; i<size1(); ++i){
      for(int j=0; j<y_trans.size1(); ++j){ // loop over the column of the resulting matrix
        int el1 = x_rowind[i];
        int el2 = y_colind[j];
        d.clear();
        while(el1 < x_rowind[i+1] && el2 < y_colind[j+1]){ // loop over non-zero elements
          int j1 = x_col[el1];
          int i2 = y_row[el2];      
          if(j1==i2){
            d.push_back(pair<int,int>(el1++,el2++));
          } else if(j1<i2) {
            el1++;
          } else {
            el2++;
          }
        }
        if(!d.empty()){
          c.push_back(j);
          mapping.push_back(d);
        }
      }
      r[i+1] = c.size();
    }
  } else {
    // loop over the row of the resulting matrix)
    for(int i=0; i<size1(); ++i){
      for(int j=0; j<y_trans.size1(); ++j){ // loop over the column of the resulting matrix
        int el1 = x_rowind[i];
        int el2 = y_colind[j];
        while(el1 < x_rowind[i+1] && el2 < y_colind[j+1]){ // loop over non-zero elements
          int j1 = x_col[el1];
          int i2 = y_row[el2];
          if(j1==i2){
            c.push_back(j);
            break;
          } else if(j1<i2) {
            el1++;
          } else {
            el2++;
          }
        }
      }
      r[i+1] = c.size();
    }
  }
    
  return ret;
}

bool CRSSparsity::operator==(const CRSSparsity& y) const{
  // Quick true if the objects are the same
  if(get() == y.get())
    return true;
  
  // First check dimensions and number of non-zeros
  if(size()!=y.size() || size1()!=y.size1() || size2()!=y.size2())
    return false;

  // Check if dense
  if(size()==numel())
    return true;
  
  // Check the number of non-zeros per row
  if(!equal(rowind().begin(),rowind().end(),y.rowind().begin()))
    return false;
  
  // Finally check the column indices
  if(!equal(col().begin(),col().end(),y.col().begin()))
    return false;
  
  // Equal if reached this point
  return true;
}

void CRSSparsity::reserve(int nnz, int nrow){
  colRef().reserve(nnz);
  rowindRef().reserve(nrow+1);
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


