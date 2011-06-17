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
  casadi_assert_message(i<size1() && j<size2(),"Indices out of bounds");
  
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

CRSSparsity CRSSparsity::scalarSparsity(1,1,true);

CRSSparsity CRSSparsity::scalarSparsitySparse(1,1,false);

CRSSparsity CRSSparsity::emptySparsity(0,0,true);

void CRSSparsity::enlarge(int nrow, int ncol, const vector<int>& ii, const vector<int>& jj){
  // Assert dimensions
  casadi_assert(ii.size() == size1());
  casadi_assert(jj.size() == size2());
  
    // Update dimensions
  (*this)->nrow_ = nrow;
  (*this)->ncol_ = ncol;

  // Begin by sparsify the columns
  vector<int>& c = colRef();
  for(int k=0; k<c.size(); ++k){
    c[k] = jj[c[k]];
  }
    
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

int CRSSparsity::depth_first_search(int j, int top, int *xi, int *pstack, const int *pinv){
  int i, p, p2, done, jnew, head = 0, *Gp, *Gi;
  Gp = &rowindRef()[0] ; Gi = &colRef()[0];
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


} // namespace CasADi


