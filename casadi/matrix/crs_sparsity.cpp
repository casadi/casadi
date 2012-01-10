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

#include "crs_sparsity_internal.hpp"
#include "sparsity_tools.hpp"
#include "../matrix/matrix.hpp"
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
  assignNode(new CRSSparsityInternal(nrow, ncol, col, rowind));
}

CRSSparsity::CRSSparsity(int nrow, int ncol, vector<int> col, vector<int> rowind){
  assignNode(new CRSSparsityInternal(nrow, ncol, col, rowind));
}
    
CRSSparsityInternal* CRSSparsity::operator->(){
  makeUnique();
  return (CRSSparsityInternal*)(SharedObject::operator->());
}

const CRSSparsityInternal* CRSSparsity::operator->() const{
  return (const CRSSparsityInternal*)(SharedObject::operator->());
}
  
bool CRSSparsity::checkNode() const{
  return dynamic_cast<const CRSSparsityInternal*>(get())!=0;
}

int CRSSparsity::size1() const{
  return (*this)->nrow_;
}
    
int CRSSparsity::size2() const{
  return (*this)->ncol_;
}
    
int CRSSparsity::numel() const{
  return (*this)->numel();
}
    
int CRSSparsity::size() const{
  return (*this)->size();
}
    
std::pair<int,int> CRSSparsity::shape() const{
  return (*this)->shape();
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
  makeUnique();
  (*this)->resize(nrow,ncol);
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

bool CRSSparsity::hasNZ(int i, int j) const {
  return (*this)->getNZ(i,j)!=-1;
}


int CRSSparsity::getNZ(int i, int j) const{
  return (*this)->getNZ(i,j);
}

CRSSparsity CRSSparsity::reshape(int n, int m) const{
  return (*this)->reshape(n,m);
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
  return (*this)->getNZ(ii,jj);
}

bool CRSSparsity::dense() const{
  return (*this)->dense();
}

bool CRSSparsity::diagonal() const{
  return (*this)->diagonal();
}

CRSSparsity CRSSparsity::getSub(const vector<int>& ii, const vector<int>& jj, vector<int>& mapping) const{
  return (*this)->getSub(ii,jj,mapping);
}

vector<int> CRSSparsity::erase(const vector<int>& ii, const vector<int>& jj){
  makeUnique();
  return (*this)->erase(ii,jj);
}

int CRSSparsity::sizeU() const{
  return (*this)->sizeU();
}

int CRSSparsity::sizeL() const{
  return (*this)->sizeL();
}


std::vector<int> CRSSparsity::getRow() const{
  return (*this)->getRow();
}

void CRSSparsity::getSparsityCRS(vector<int>& rowind, vector<int> &col) const{
  rowind = this->rowind();
  col = this->col();
}

void CRSSparsity::getSparsity(vector<int>& row, vector<int> &col) const{
  row = this->getRow();
  col = this->col();
}

CRSSparsity CRSSparsity::transpose(vector<int>& mapping) const{
  return (*this)->transpose(mapping);
}

CRSSparsity CRSSparsity::transpose() const{
  return (*this)->transpose();
}

CRSSparsity CRSSparsity::patternUnion(const CRSSparsity& y, vector<unsigned char>& mapping, bool f00_is_zero, bool f0x_is_zero, bool fx0_is_zero) const{
  return (*this)->patternUnion(y, mapping, f00_is_zero, f0x_is_zero, fx0_is_zero);
}

CRSSparsity CRSSparsity::patternIntersection(const CRSSparsity& y, vector<unsigned char>& mapping) const{
  return CRSSparsity();
}

CRSSparsity CRSSparsity::patternProduct(const CRSSparsity& y_trans) const{
  return (*this)->patternProduct(y_trans);
}

CRSSparsity CRSSparsity::patternProduct(const CRSSparsity& y_trans, vector< vector< pair<int,int> > >& mapping) const{
  return (*this)->patternProduct(y_trans,mapping);
}

bool CRSSparsity::operator==(const CRSSparsity& y) const{
  return (*this)->isEqual(y);
}

CRSSparsity CRSSparsity::operator+(const CRSSparsity& b) const {
  return (DMatrix(*this,1)+DMatrix(b,1)).sparsity();
}

void CRSSparsity::reserve(int nnz, int nrow){
  makeUnique();
  (*this)->reserve(nnz,nrow);
}

void CRSSparsity::append(const CRSSparsity& sp){
  casadi_assert(this!=&sp); // NOTE: this case needs to be handled
  makeUnique();
  (*this)->append(sp);
}

CRSSparsity CRSSparsity::scalarSparsity(1,1,true);

CRSSparsity CRSSparsity::scalarSparsitySparse(1,1,false);

CRSSparsity CRSSparsity::emptySparsity(0,0,true);

void CRSSparsity::enlarge(int nrow, int ncol, const vector<int>& ii, const vector<int>& jj){
  enlargeRows(nrow,ii);
  enlargeColumns(ncol,jj);
}

void CRSSparsity::enlargeRows(int nrow, const std::vector<int>& ii){
  makeUnique();
  (*this)->enlargeRows(nrow,ii);
}

void CRSSparsity::enlargeColumns(int ncol, const std::vector<int>& jj){
  makeUnique();
  (*this)->enlargeColumns(ncol,jj);
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

CRSSparsity CRSSparsity::makeDense(std::vector<int>& mapping) const{
  return (*this)->makeDense(mapping);
}

std::string CRSSparsity::dimString() 	const { 
  return (*this)->dimString();
}

CRSSparsity CRSSparsity::diag(std::vector<int>& mapping) const{
  return (*this)->diag(mapping);
}

std::vector<int> CRSSparsity::eliminationTree(bool ata) const{
  return (*this)->eliminationTree(ata);
}

int CRSSparsity::depthFirstSearch(int j, int top, std::vector<int>& xi, std::vector<int>& pstack, const std::vector<int>& pinv, std::vector<bool>& marked) const{
  return (*this)->depthFirstSearch(j,top,xi,pstack,pinv,marked);
}

int CRSSparsity::stronglyConnectedComponents(std::vector<int>& p, std::vector<int>& r) const{
  return (*this)->stronglyConnectedComponents(p,r);
}

int CRSSparsity::dulmageMendelsohn(std::vector<int>& rowperm, std::vector<int>& colperm, std::vector<int>& rowblock, std::vector<int>& colblock, std::vector<int>& coarse_rowblock, std::vector<int>& coarse_colblock, int seed) const{
  return (*this)->dulmageMendelsohn(rowperm, colperm, rowblock, colblock, coarse_rowblock, coarse_colblock, seed);
}

bool CRSSparsity::columnsSequential(bool strictly) const{
  return (*this)->columnsSequential(strictly);
}

void CRSSparsity::removeDuplicates(std::vector<int>& mapping){
  makeUnique();
  (*this)->removeDuplicates(mapping);
}

std::vector<int> CRSSparsity::getElementMapping() const{
  return (*this)->getElementMapping();
}

} // namespace CasADi
