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
#include "sparsity_tools.hpp"
#include "../stl_vector_tools.hpp"
#include <climits>

using namespace std;

namespace CasADi{

CRSSparsity::CRSSparsity(int null){
  casadi_assert(null==0);
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
  return dynamic_cast<const CRSSparsityNode*>(get())!=0;
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

void CRSSparsity::sanityCheck(bool complete) const { 
  (*this)->sanityCheck(complete);
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
  casadi_assert_message(i<size1() && j<size2(),"Indices out of bounds");

  if (i<0) i += size1();
  if (j<0) j += size2();
  
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
  casadi_assert_message(i<size1(),"First index out of bounds");
  casadi_assert_message(j<size2(),"Second index out of bounds");
  
  if (i<0) i += size1();
  if (j<0) j += size2();
  
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

CRSSparsity CRSSparsity::reshape(int n, int m) const{
  casadi_assert_message(numel() == n*m, "reshape: number of elements must remain the same");
  CRSSparsity ret(n,m);
  ret.reserve(size(), n);
  for(int i=0; i<size1(); ++i){
    for(int el=rowind(i); el<rowind(i+1); ++el){
      int j = col(el);
      
      // Element number
      int k_ret = j+i*size2();
      
      // Row and column in the new matrix
      int i_ret = k_ret/m;
      int j_ret = k_ret%m;
      ret.getNZ(i_ret,j_ret);
    }
  }
  return ret;
}

// vector<int> CRSSparsity::getNZNew(vector<int> i, vector<int> j){
//   vector<int> ret;
//   ret.reserve(i.size());
// 
//     // Quick return if matrix is dense
//   if(numel()==size()){
//     for(int k=0; k<i.size(); ++k)
//       ret.push_back(j[k]+i[k]*size2());
//     return ret;
//   }
// 
//   // Very inefficient algorithm
//   for(int k=0; k<i.size(); ++k){
//     ret.push_back(getNZ(i[k],j[k]));
//   }
//   return ret;
// }
// 
// vector<int> CRSSparsity::getNZNew(vector<int> i, vector<int> j) const{
//   vector<int> ret;
//   ret.reserve(i.size());
// 
//     // Quick return if matrix is dense
//   if(numel()==size()){
//     for(int k=0; k<i.size(); ++k)
//       ret.push_back(j[k]+i[k]*size2());
//     return ret;
//   }
// 
//   // Very inefficient algorithm
//   for(int k=0; k<i.size(); ++k){
//     ret.push_back(getNZ(i[k],j[k]));
//   }
//   return ret;
// }

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

bool CRSSparsity::diagonal() const{
  // Check if matrix is square
  if(size1() != size2()) return false;
    
  // Check if correct number of non-zeros (one per row)
  if(size() != size1()) return false;

  // Check that the column indices are correct
  for(int i=0; i<size(); ++i){
    if(col(i)!=i)
      return false;
  }
   
  // Make sure that the row indices are correct
  for(int i=0; i<size1(); ++i){
    if(rowind(i)!=i)
      return false;
  }
  
  // Diagonal if reached this point
  return true;
}


CRSSparsity CRSSparsity::getSub(const vector<int>& ii, const vector<int>& jj, vector<int>& mapping) const{
  // Get non-zeros
  vector<int> kk = getNZ(ii,jj);

  // count the number of non-zeros
  int nnz = 0;
  for(int k=0; k<kk.size(); ++k)
    if(kk[k]!=-1)
      nnz++;
      
  // Allocate sparsity vectors
  vector<int> rowind(ii.size()+1,0);
  vector<int> col(nnz);
  mapping.resize(nnz);
  
  // Get sparsity
  int k=0;
  int el=0;
  for(int i=0; i<ii.size(); ++i){
    for(int j=0; j<jj.size(); ++j){
      if(k<kk.size()){
        if(kk[k]!=-1){
          col[el] = j;
          mapping[el] = kk[k];
          el++;
        }
        k++;
      }
    }
    rowind[i+1]=el;
  }
  
  // Create sparsity pattern
  CRSSparsity sp(ii.size(),jj.size(),col,rowind);
  return sp;
}

vector<int> CRSSparsity::erase(const vector<int>& ii, const vector<int>& jj){
  // Mapping
  vector<int> mapping;
  
  // Quick return if no elements
  if(numel()==0)
    return mapping;
  
  // Reserve memory
  mapping.reserve(size());
  
  // References to sparsity pattern
  vector<int>& c = colRef();
  vector<int>& r = rowindRef();
  
  // Number of non-zeros
  int nz=0;
  
  // Rows to be erased
  vector<int>::const_iterator ie = ii.begin();
  
  // First and last index for the row
  int el_first=0, el_last=0;
  
  // Loop over rows
  for(int i=0; i<size1(); ++i){
    // Update beginning and end of non-zero indices
    el_first = el_last;
    el_last = r[i+1];
    
    // Is it a row that can be deleted
    bool deletable_row = ie!=ii.end() && *ie==i;
    if(deletable_row){
      ie++;
      
      // Columns to be erased
      vector<int>::const_iterator je = jj.begin();

      // Loop over nonzero elements of the row
      for(int el=el_first; el<el_last; ++el){
        // Column
        int j=c[el];
        
        // Continue to the next column to skip
        for(; je!=jj.end() && *je<j; ++je);
        
        // Remove column if necessary
        if(je!=jj.end() && *je==j){
          je++;
          continue;
        }
        
        // Save old nonzero for each new nonzero
        mapping.push_back(el);
        
        // Update column and increase nonzero counter
        c[nz++] = j;
      }
    } else {
      // Loop over nonzero elements of the row
      for(int el=el_first; el<el_last; ++el){
        // Column
        int j=c[el];
      
        // Save old nonzero for each new nonzero
        mapping.push_back(el);
      
        // Update column and increase nonzero counter
        c[nz++] = j;
      }
    }
    
    // Register last nonzero of the row
    r[i+1]=nz;
  }
  
  // Truncate column matrix
  c.resize(nz);
  
  return mapping;
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

void CRSSparsityNode::sanityCheck(bool complete) const{
  if (rowind_.size() != nrow_+1) {
    std::stringstream s;
    s << "CRSSparsityNode:Compressed Row Storage is not sane. The following must hold:" << std::endl;
    s << "  rowind.size() = nrow + 1, but got   rowind.size() = " << rowind_.size() << "   and   nrow = "  << nrow_ << std::endl;
    s << "  Note that the signature is as follows: CRSSparsity (nrow, ncol, col, rowind)." << std::endl;
    throw CasadiException(s.str());
  }
  if (complete) {
  
    if (rowind_.size()>0) {
      for (int k=1;k<rowind_.size();k++) {
        if (rowind_[k]<rowind_[k-1]) {
          throw CasadiException("CRSSparsityNode:Compressed Row Storage is not sane. rowind must be monotone. Note that the signature is as follows: CRSSparsity (nrow, ncol, col, rowind).");
        }
      }
      
      if (rowind_[0]!=0) {
        throw CasadiException("CRSSparsityNode:Compressed Row Storage is not sane. First element of rowind must be zero. Note that the signature is as follows: CRSSparsity (nrow, ncol, col, rowind).");
      }
      if (rowind_[(rowind_.size()-1)]!=col_.size()) {
        std::stringstream s;
        s << "CRSSparsityNode:Compressed Row Storage is not sane. The following must hold:" << std::endl;
        s << "  rowind[lastElement] = col.size(), but got   rowind[lastElement] = " << rowind_[(rowind_.size()-1)] << "   and   col.size() = "  << col_.size() << std::endl;
        s << "  Note that the signature is as follows: CRSSparsity (nrow, ncol, col, rowind)." << std::endl;
        throw CasadiException(s.str());
      }
      if (col_.size()>nrow_*ncol_) {
        std::stringstream s;
        s << "CRSSparsityNode:Compressed Row Storage is not sane. The following must hold:" << std::endl;
        s << "  col.size() <= nrow * ncol, but got   col.size()  = " << col_.size() << "   and   nrow * ncol = "  << nrow_*ncol_ << std::endl;
        s << "  Note that the signature is as follows: CRSSparsity (nrow, ncol, col, rowind)." << std::endl;
        throw CasadiException(s.str());
      }
    }
    for (int k=0;k<col_.size();k++) {
      if (col_[k]>=ncol_ || col_[k] < 0) {
        std::stringstream s;
        s << "CRSSparsityNode:Compressed Row Storage is not sane. The following must hold:" << std::endl;
        s << "  0 <= col[i] < ncol for each i, but got   col[i] = " << col_[k] << "   and   ncol = "  << ncol_ << std::endl;
        s << "  Note that the signature is as follows: CRSSparsity (nrow, ncol, col, rowind)." << std::endl;
        throw CasadiException(s.str());
      }
    }
  
  }
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
  
  // Get the sparsity of the transpose in sparse triplet form
  const vector<int>& trans_row = col();
  vector<int> trans_col = getRow();

  // Create the sparsity pattern
  return sp_triplet(size2(),size1(),trans_row,trans_col,mapping);

}

#if 0
CRSSparsity CRSSparsity::patternUnion(const CRSSparsity& y, std::vector<unsigned char>& mapping) const{
  // Shorthand
  const CRSSparsity& x = *this;

  // Assert dimensions
  casadi_assert_message(x.numel()==1 || y.numel()==1 || x==y, "Dimensions does not match");
  
  // Return sparsity
  CRSSparsity ret;
  
  // Treat simple cases first
  if(x==y){
    ret = x;
    mapping.resize(x.size(),1 | 2 | 4 | 8 );
    mapping.back() ^=  4 | 8;
  } else if(x.size()==0){
    ret = y;
    mapping.resize(y.size(),2 | 8 );
    mapping.back() ^= 8;
  } else if(y.size()==0){
    ret = x;
    mapping.resize(x.size(),1 | 4);
    mapping.back() ^= 4;
/*  } else if(x.numel()==1){
    ret = 
    
    
    
    
    
    
    ret = y;
    if(x.size()==1){
      mapping.resize(y.size(),1 | 2 | 8);
    } else {
      mapping.resize(y.size(),2 | 
    */
    
  }

  
  
  
  
  
  
  
  // Create return object
  int sz1 = std::max(x.size1(),y.size1());
  int sz2 = std::max(x.size2(),y.size2());
  ret = CRSSparsity(sz1,sz2);

  // Get refences to the sparsity vectors
  vector<int>& r = ret.rowindRef();
  vector<int>& c = ret.colRef();
  
  // Prepare the assembly of the rowind vector below
  r.clear();
  r.push_back(0);
  
  // Clear the mapping
/*  mapping.clear();
  
  // If the first argument is a scalar
  if(x.numel()==1){
    if(y.dense()){
      ret = y;
      mapping.resize(y.size(), 1 | 2 | 8);
    } else if(f00_is_zero)
      
    
    
  }*/
  
  
//   // Quick return if the patterns are equal and f00_is_zero
//   if(x == y && (x.dense() || f00_is_zero)){
//     mapping.resize(x.size());
//     fill(mapping.begin(),mapping.end(), 1 | 2 | 4 | 8 );
//     mapping.back() ^= 4 | 8; // do not increase pointer when already at the end
//     return x;
//   }
//   
//   // Quick return if the first argument is a scalar and f00_is_zero
//   if(x.numel()==1 && (y.dense() || f00_is_zero)){
//     mapping.resize(y.size());
//     fill(mapping.begin(),mapping.end(), 1 | 2 | 8 );
//     mapping.back() ^= 8; // do not increase pointer when already at the end
//     return x;
//   }
//   
//   // Quick return if the second argument is a scalar and f00_is_zero
//   if(y.numel()==1 && (x.dense() || f00_is_zero)){
//     mapping.resize(y.size());
//     fill(mapping.begin(),mapping.end(), 1 | 2 | 4 );
//     mapping.back() ^= 4; // do not increase pointer when already at the end
//     return x;
//   }
  
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
      int col1 = el1<el1_last ? col(el1) : size2();
      int col2 = el2<el2_last ? y.col(el2) : size2();

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
#endif

CRSSparsity CRSSparsity::patternUnion(const CRSSparsity& y, vector<unsigned char>& mapping, bool f00_is_zero, bool f0x_is_zero, bool fx0_is_zero) const{

  // Assert dimensions
  casadi_assert_message(size1()==y.size1(), "The number of rows does not match");
  casadi_assert_message(size2()==y.size2(), "The number of columns does not match");

  // Return object
  CRSSparsity ret;

  // Quick intersection if the patterns are equal
  if(*this == y){
    ret = *this;
    mapping.resize(size());
    fill(mapping.begin(),mapping.end(), 1 | 2);
  } else {
  
    // Create return object
    ret = CRSSparsity(size1(),size2());
    
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
        int col1 = el1<el1_last ? col(el1) : size2();
        int col2 = el2<el2_last ? y.col(el2) : size2();

        // Add to the return matrix
        if(col1==col2){ //  both nonzero
          c.push_back(col1);
          mapping.push_back( 1 | 2);
          el1++; el2++;
        } else if(col1<col2){ //  only first argument is nonzero
          if(!fx0_is_zero){
            c.push_back(col1);
            mapping.push_back(1);
          } else {
            mapping.push_back(1 | 4);
          }
          el1++;
        } else { //  only second argument is nonzero
          if(!f0x_is_zero){
            c.push_back(col2);
            mapping.push_back(2);
          } else {
            mapping.push_back(2 | 4);
          }
          el2++;
        }
      }
      
      // Save the index of the last nonzero on the row
      r.push_back(c.size());
    }
  }
  
  // Check if we need to add extra nonzeros
  if(f00_is_zero || ret.dense()){
    // No nonzero elements in the return object or sparse entries evaluating to 0
    return ret;
  } else {
    // Create a new sparsity pattern with the remaining nonzeros added
    vector<unsigned char> mapping_dense(ret.numel(),0); // FIXME: this allocation can be avoided, just iterate in reverse order instead
    
    // Loop over rows
    for(int i=0; i<ret.size1(); ++i){
      // Loop over nonzeros
      for(int el=ret.rowind(i); el<ret.rowind(i+1); ++el){
        // Get column
        int j=ret.col(el);
        
        // Get the nonzero of the dense matrix
        int el_dense = j+i*ret.size2();
        
        // Save to mapping
        mapping_dense[el_dense] = mapping[el];
      }
    }
    
    // Use the dense mapping instead and return a dense sparsity
    mapping_dense.swap(mapping);
    return CRSSparsity(ret.size1(),ret.size2(),true);
  }
}

CRSSparsity CRSSparsity::patternIntersection(const CRSSparsity& y, vector<unsigned char>& mapping) const{
  return CRSSparsity();
}

CRSSparsity CRSSparsity::patternProduct(const CRSSparsity& y_trans) const{
  // Dimensions
  int x_nrow = size1();
  int y_ncol = y_trans.size1();

  // Quick return if both are dense
  if(dense() && y_trans.dense()){
    return CRSSparsity(x_nrow,y_ncol,true);
  }
  
  // return object
  CRSSparsity ret(x_nrow,y_ncol);
  
  // Get the vectors for the return pattern
  vector<int>& c = ret.colRef();
  vector<int>& r = ret.rowindRef();
  
  // Direct access to the arrays
  const vector<int> &x_col = col();
  const vector<int> &y_row = y_trans.col();
  const vector<int> &x_rowind = rowind();
  const vector<int> &y_colind = y_trans.rowind();

  // If the compiler supports C99, we shall use the long long datatype, which is 64 bit, otherwise long
#if __STDC_VERSION__ >= 199901L
  typedef unsigned long long int_t;
#else
  typedef unsigned long int_t;
#endif
  
  // Number of directions we can deal with at a time
  int nr = CHAR_BIT*sizeof(int_t); // the size of int_t in bits (CHAR_BIT is the number of bits per byte, usually 8)

  // Number of such groups needed
  int ng = x_nrow/nr;
  if(ng*nr != x_nrow) ng++;
  
  // Which columns exist in a row of the first factor
  vector<int_t> in_x_row(size2());
  vector<int_t> in_res_col(y_ncol);

  // Loop over the rows of the resulting matrix, nr rows at a time
  for(int rr=0; rr<ng; ++rr){

    // Mark the elements in the x row
    fill(in_x_row.begin(),in_x_row.end(),0); // NOTE: expensive?
    int_t b=1;
    for(int i=rr*nr; i<rr*nr+nr && i<x_nrow; ++i){
      for(int el1=x_rowind[i]; el1<x_rowind[i+1]; ++el1){
        in_x_row[x_col[el1]] |= b;
      }
      b <<= 1;
    }

    // Get the sparsity pattern for the set of rows
    fill(in_res_col.begin(),in_res_col.end(),0); // NOTE: expensive?
    for(int j=0; j<y_ncol; ++j){
      
      // Loop over the nonzeros of the column of the second factor
      for(int el2=y_colind[j]; el2<y_colind[j+1]; ++el2){
        
        // Get the row
        int i_y = y_row[el2];
        
        // Add nonzero if the element matches an element in the x row
        in_res_col[j] |= in_x_row[i_y];
      }
    }

    b = 1;
    for(int i=rr*nr; i<rr*nr+nr && i<x_nrow; ++i){
      
      // loop over the columns of the resulting matrix
      for(int j=0; j<y_ncol; ++j){
        
        // Save nonzero, if any
        if(in_res_col[j] & b){
          c.push_back(j);
        }
      }
      r[i+1] = c.size();
      b <<=1;
    }
  }
  return ret;
}

CRSSparsity CRSSparsity::patternProduct(const CRSSparsity& y_trans, vector< vector< pair<int,int> > >& mapping) const{
  // return object
  CRSSparsity ret = patternProduct(y_trans);
  
  // Get the vectors for the return pattern
  const vector<int>& c = ret.col();
  const vector<int>& r = ret.rowind();
  
  // Direct access to the arrays
  const vector<int> &x_col = col();
  const vector<int> &y_row = y_trans.col();
  const vector<int> &x_rowind = rowind();
  const vector<int> &y_colind = y_trans.rowind();

  // Clear the mapping
  mapping.resize(ret.size());

  // the entry of the matrix to be calculated
  vector< pair<int,int> > d;

  // loop over the row of the resulting matrix)
  for(int i=0; i<size1(); ++i){
    // Loop over nonzeros
    for(int el=r[i]; el<r[i+1]; ++el){
      int j = c[el];
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
      mapping[el] = d;
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

CRSSparsity CRSSparsity::scalarSparsity(1,1,true);

CRSSparsity CRSSparsity::scalarSparsitySparse(1,1,false);

CRSSparsity CRSSparsity::emptySparsity(0,0,true);

void CRSSparsity::enlarge(int nrow, int ncol, const vector<int>& ii, const vector<int>& jj){
  enlargeRows(nrow,ii);
  enlargeColumns(ncol,jj);
}

void CRSSparsity::enlargeRows(int nrow, const std::vector<int>& ii){
  // Assert dimensions
  casadi_assert(ii.size() == size1());

  // Update dimensions
  (*this)->nrow_ = nrow;

  // Sparsify the rows
  vector<int>& r = rowindRef();
  r.resize(nrow+1,size());
  int ik=ii.back(); // need only to update from the last new index
  int nz=size(); // number of nonzeros up till this row
  for(int i=ii.size()-1; i>=0; --i){
    // Update rowindex for new rows
    for(; ik>ii[i]; --ik){
      r[ik] = nz;
    }
    
    // Update non-zero counter
    nz = r[i];
    
    // Update rowindex for old rows
    r[ii[i]] = nz;
  }
  
  // Append zeros to the beginning
  for(; ik>=0; --ik){
    r[ik] = 0;
  }
}

void CRSSparsity::enlargeColumns(int ncol, const std::vector<int>& jj){
  // Assert dimensions
  casadi_assert(jj.size() == size2());
  
    // Update dimensions
  (*this)->ncol_ = ncol;

  // Begin by sparsify the columns
  vector<int>& c = colRef();
  for(int k=0; k<c.size(); ++k){
    c[k] = jj[c[k]];
  }
}

int CRSSparsity::depth_first_search(int j, int top, int *xi, int *pstack, const int *pinv){
  int i, p, p2, done, jnew, head = 0, *Gp, *Gi;
  Gp = getPtr(rowindRef()) ; Gi = getPtr(colRef());
  xi[0] = j;                // initialize the recursion stack
  while (head >= 0)
  {
    j = xi[head] ;         // get j from the top of the recursion stack
    jnew = pinv ? pinv[j] : j ;
    if (!marked(Gp, j))
    {
        mark(Gp, j) ;       // mark node j as visited
        pstack[head] = (jnew < 0) ? 0 : unflip(Gp[jnew]) ;
    }
    done = 1 ;                  // node j done if no unvisited neighbors
    p2 = (jnew < 0) ? 0 : unflip(Gp[jnew+1]) ;
    for (p = pstack[head] ; p < p2 ; p++)  // examine all neighbors of j 
    {
        i = Gi[p] ;            // consider neighbor node i 
        if (marked(Gp, i)) continue ;   // skip visited node i 
        pstack[head] = p ;     // pause depth-first search of node j 
        xi[++head] = i ;       // start dfs at node i 
        done = 0 ;              // node j is not done 
        break ;                 // break, to start dfs (i) 
    }
    if (done)               // depth-first search at node j is done 
    {
        head-- ;            // remove j from the recursion stack 
        xi[--top] = j ;    // and place in the output stack 
    }
  }
  return top;
}

CRSSparsity CRSSparsity::createDiagonal(int n){
  return createDiagonal(n,n);
}

CRSSparsity CRSSparsity::createDiagonal(int n, int m){
  CRSSparsity ret(n,m);
  
  // Set columns
  vector<int> &c = ret.colRef();
  c.resize(min(n,m));
  for(int i=0; i<c.size(); ++i)
    c[i] = i;
  
  // Set row indices
  vector<int> &r = ret.rowindRef();
  for(int i=0; i<n && i<m; ++i)
    r[i] = i;
  
  for(int i=min(n,m); i<n+1; ++i)
    r[i] = c.size();
  
  return ret;
}


void CRSSparsity::strongly_connected_components(){
  // not working
  #if 0
  int n, i, k, b, nb = 0, top, *xi, *pstack, *p, *r, *Ap, *ATp, *rcopy, *Blk;
    
  vector<int> AT_mapping;
  CRSSparsity AT = transpose(AT_mapping);
  n = size1() ; Ap = &rowindRef().front();
  
  vector<int> xi_data(2*n+1);
  xi = &xi_data.front();
  
  vector<int> Dp(n);
  vector<int> Dq(0);
  vector<int> Dr(n+6);
  vector<int> Ds(6);
  int Dnb;
  int Drr[5];
  int Dcc[5];
  
  Blk = xi ; rcopy = pstack = xi + n ;
  p = &Dp.front() ; r = &Dr.front() ; ATp = &AT.rowindRef().front();
  
  top = n ;
  for (i = 0 ; i < n ; i++){
      if (!marked(Ap, i)) top = depth_first_search(i, top, xi, pstack, NULL);
  }
  
  for (i = 0 ; i < n ; i++) 
    mark(Ap, i);
  top = n;
  nb = n;
  for (k = 0 ; k < n ; k++){
    i = xi [k] ;
    if (marked(ATp, i)) continue ;
    r [nb--] = top ;
    top = AT.depth_first_search(i, top, p, pstack, NULL) ;
  }
  r[nb] = 0;
  for (k = nb ; k <= n ; k++) r[k-nb] = r[k];
  Dnb = nb = n-nb;
  for (b = 0 ; b < nb ; b++){
      for (k = r [b] ; k < r [b+1] ; k++) Blk [p [k]] = b ;
  }
  for (b = 0 ; b <= nb ; b++) rcopy [b] = r [b] ;
  for (i = 0 ; i < n ; i++) p [rcopy [Blk [i]]++] = i ;
  
  cout << Dp << endl;
  cout << Dr << endl;
  cout << Dnb << endl;
  cout << rowind() << endl;
  cout << col() << endl;
  cout << "end" << endl;
  #endif //0
  
}

CRSSparsity CRSSparsity::makeDense(std::vector<int>& mapping) const{
  mapping.resize(size());
  for(int i=0; i<size1(); ++i){
    for(int el=rowind(i); el<rowind(i+1); ++el){
      int j = col(el);
      mapping[el] = j + i*size2();
    }
  }
  
  return CRSSparsity(size1(),size2(),true);
}

std::string CRSSparsity::dimString() 	const { 
  std::stringstream ss;
  ss << "(" << size1() << "x" << size2() << "=" << numel() << "|" << size() << ")";
  return ss.str();
}

CRSSparsity CRSSparsity::diag(std::vector<int>& mapping) const{
  if (size1()==size2()) {
    // Return object
    CRSSparsity ret(0,1);
    ret.reserve(std::min(size(),size1()),size1());
    
    // Mapping
    mapping.clear();
    
    // Loop over nonzero
    for(int i=0; i<size1(); ++i){
      
      // Enlarge the return matrix
      ret.resize(i+1,1);
    
      // Get to the right nonzero of the row
      int el = rowind(i);
      while(el<rowind(i+1) && col(el)<i){
        el++;
      }
      
      if (el>=size()) return ret;
      
      // Add element if nonzero on diagonal
      if(col(el)==i){
        ret.getNZ(i,0);
        mapping.push_back(el);
      }
    }
    
    return ret;
    
  } else if (size1()==1 || size2()==1) {
    std::vector<int> dummy;
    
    // Have a row vector
    const CRSSparsity &sp = (size1() == 1) ? (*this) : (*this).transpose(dummy);
    
    // Return object
    CRSSparsity ret(sp.size2(),sp.size2());
    ret.reserve(size(),sp.size2());
    
    mapping.clear();
    mapping.resize(size());
        
    // Loop over nonzero
    for(int k=0;k<size();k++) {
      mapping[k]=k; // mapping will just be a range(size())
      
      int i = sp.col()[k];
      
      ret.getNZ(i,i); // Create a nonzero into the ret sparsity pattern
    }
    
    return ret;
  } else {
    stringstream s;
    s << "diag: wrong argument shape. Expecting square matrix or vector-like, but got " << dimString() << " instead." <<  std::endl;
    throw CasadiException(s.str());
  }
}

std::vector<int> CRSSparsity::eliminationTree(bool ata) const{
  return (*this)->eliminationTree(ata);
}

std::vector<int> CRSSparsityNode::eliminationTree(bool ata) const{
  
  // Allocate result
  vector<int> parent(nrow_);
  
   // Allocate workspace 
  vector<int> ancestor(nrow_);
  vector<int> prev(ata ? ncol_ : 0, -1);
  
  // Loop over rows
  for(int k=0; k<nrow_; ++k){
    // Start with no parent or ancestor
    parent[k] = -1;
    ancestor[k] = -1;
    
    // Loop over nonzeros
    for(int p=rowind_[k]; p<rowind_[k+1]; ++p){
      
      // What is this?
      int i=ata ? (prev[col_[p]]) : (col_[p]);
      
      // Transverse from i to k
      while(i!=-1 && i<k){
        
        // Next i is the ancestor of i
        int inext = ancestor[i];

        // Path compression
        ancestor[i] = k;
        
        // No ancestor, parent is k
        if(inext==-1) 
          parent[i] = k;
        
        // Update i
        i=inext;
      }
      
      // What is this?
      if(ata){
        prev[col_[p]] = k;
      }
    }
  }
  
  return parent;
  
}

} // namespace CasADi


