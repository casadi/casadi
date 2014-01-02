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

  // Singletons
  class EmptySparsity : public CRSSparsity{  
  public:
    EmptySparsity(){
      vector<int> col,rowind(1,0);
      assignNode(new CRSSparsityInternal(0,0,col,rowind));
    }
  };

  class ScalarSparsity : public CRSSparsity{  
  public:
    ScalarSparsity(){
      vector<int> col(1,0),rowind(2);
      rowind[0] = 0;
      rowind[1] = 1;
      assignNode(new CRSSparsityInternal(1,1,col,rowind));
    }
  };

  class ScalarSparseSparsity : public CRSSparsity{  
  public:
    ScalarSparseSparsity(){
      vector<int> col,rowind(2,0);
      assignNode(new CRSSparsityInternal(1,1,col,rowind));
    }
  };
  
  CRSSparsity::CRSSparsity(int dummy){
    casadi_assert(dummy==0);
  }
  
  CRSSparsity CRSSparsity::create(CRSSparsityInternal *node){
    CRSSparsity ret;
    ret.assignNode(node);
    return ret;
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
 
    assignCached(nrow, ncol, col, rowind);
  }

  CRSSparsity::CRSSparsity(int nrow, int ncol, const vector<int>& col, const vector<int>& rowind){
    assignCached(nrow, ncol, col, rowind);
  }

  void CRSSparsity::reCache(){
    assignCached(size1(),size2(),col(),rowind());
  }
 
  CRSSparsityInternal* CRSSparsity::operator->(){
    makeUnique();
    return static_cast<CRSSparsityInternal*>(SharedObject::operator->());
  }

  const CRSSparsityInternal* CRSSparsity::operator->() const{
    return static_cast<const CRSSparsityInternal*>(SharedObject::operator->());
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

  bool CRSSparsity::empty() const{
    return (*this)->empty();
  }
  
  bool CRSSparsity::null() const{
    return (*this)->null();
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

  vector<int> CRSSparsity::getNZ(const vector<int>& ii, const vector<int>& jj) const{
    return (*this)->getNZ(ii,jj);
  }

  bool CRSSparsity::scalar(bool scalar_and_dense) const{
    return (*this)->scalar(scalar_and_dense);
  }

  bool CRSSparsity::dense() const{
    return (*this)->dense();
  }

  bool CRSSparsity::diagonal() const{
    return (*this)->diagonal();
  }

  bool CRSSparsity::square() const{
    return (*this)->square();
  }

  CRSSparsity CRSSparsity::sub(const vector<int>& ii, const vector<int>& jj, vector<int>& mapping) const{
    return (*this)->sub(ii,jj,mapping);
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

  int CRSSparsity::sizeD() const{
    return (*this)->sizeD();
  }

  std::vector<int> CRSSparsity::getRow() const{
    return (*this)->getRow();
  }

  void CRSSparsity::getSparsityCRS(vector<int>& rowind, vector<int> &col) const{
    rowind = this->rowind();
    col = this->col();
  }

  void CRSSparsity::getSparsityCCS(std::vector<int>& row, std::vector<int> &colind) const {
    transpose().getSparsityCRS(colind,row);
  }
    

  void CRSSparsity::getSparsity(vector<int>& row, vector<int> &col) const{
    row = this->getRow();
    col = this->col();
  }

  CRSSparsity CRSSparsity::transpose(vector<int>& mapping, bool invert_mapping) const{
    return (*this)->transpose(mapping,invert_mapping);
  }

  CRSSparsity CRSSparsity::transpose() const{
    return (*this)->transpose();
  }

  CRSSparsity CRSSparsity::patternCombine(const CRSSparsity& y, bool f0x_is_zero, bool fx0_is_zero, vector<unsigned char>& mapping) const{
    return (*this)->patternCombine(y, f0x_is_zero, fx0_is_zero, mapping);
  }

  CRSSparsity CRSSparsity::patternCombine(const CRSSparsity& y, bool f0x_is_zero, bool fx0_is_zero) const{
    return (*this)->patternCombine(y, f0x_is_zero, fx0_is_zero);
  }

  CRSSparsity CRSSparsity::patternUnion(const CRSSparsity& y, vector<unsigned char>& mapping) const{
    return (*this)->patternCombine(y, false, false, mapping);
  }

  CRSSparsity CRSSparsity::patternUnion(const CRSSparsity& y) const{
    return (*this)->patternCombine(y, false, false);
  }

  CRSSparsity CRSSparsity::patternIntersection(const CRSSparsity& y, vector<unsigned char>& mapping) const{
    return (*this)->patternCombine(y, true, true, mapping);
  }

  CRSSparsity CRSSparsity::patternIntersection(const CRSSparsity& y) const{
    return (*this)->patternCombine(y, true, true);
  }

  CRSSparsity CRSSparsity::patternProduct(const CRSSparsity& y_trans) const{
    return (*this)->patternProduct(y_trans);
  }

  CRSSparsity CRSSparsity::patternProduct(const CRSSparsity& y_trans, vector< vector< pair<int,int> > >& mapping) const{
    return (*this)->patternProduct(y_trans,mapping);
  }

  bool CRSSparsity::isEqual(const CRSSparsity& y) const{
    return (*this)->isEqual(y);
  }

  bool CRSSparsity::isEqual(int nrow, int ncol, const std::vector<int>& col, const std::vector<int>& rowind) const{
    return (*this)->isEqual(nrow,ncol,col,rowind);
  }

  CRSSparsity CRSSparsity::operator+(const CRSSparsity& b) const {
    return (DMatrix(*this,1)+DMatrix(b,1)).sparsity();
  }

  CRSSparsity CRSSparsity::operator*(const CRSSparsity& b) const {
    std::vector< unsigned char > mapping;
    return patternIntersection(b, mapping);
  }
  
  CRSSparsity CRSSparsity::patternInverse() const {
    return (*this)->patternInverse();
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

  CRSSparsity::CachingMap& CRSSparsity::getCache(){
    static CachingMap ret;
    return ret;
  }

  const CRSSparsity& CRSSparsity::getScalar(){
    static ScalarSparsity ret;
    return ret;
  }

  const CRSSparsity& CRSSparsity::getScalarSparse(){
    static ScalarSparseSparsity ret;
    return ret;
  }

  const CRSSparsity& CRSSparsity::getEmpty(){
    static EmptySparsity ret;
    return ret;
  }

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

  std::string CRSSparsity::dimString()         const { 
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

  std::vector<int> CRSSparsity::getElements(bool row_major) const{
    std::vector<int> loc;
    getElements(loc,row_major);
    return loc;
  }

  void CRSSparsity::getElements(std::vector<int>& loc, bool row_major) const{
    (*this)->getElements(loc,row_major);
  }

  void CRSSparsity::getNZInplace(std::vector<int>& indices) const{
    (*this)->getNZInplace(indices);
  }

  CRSSparsity CRSSparsity::unidirectionalColoring(const CRSSparsity& AT, int cutoff) const{
    if(AT.isNull()){
      return (*this)->unidirectionalColoring(transpose(),cutoff);
    } else {
      return (*this)->unidirectionalColoring(AT,cutoff);
    }
  }

  CRSSparsity CRSSparsity::starColoring(int ordering, int cutoff) const{
    return (*this)->starColoring(ordering,cutoff);
  }

  CRSSparsity CRSSparsity::starColoring2(int ordering, int cutoff) const{
    return (*this)->starColoring2(ordering,cutoff);
  }

  std::vector<int> CRSSparsity::largestFirstOrdering() const{
    return (*this)->largestFirstOrdering();
  }

  CRSSparsity CRSSparsity::pmult(const std::vector<int>& p, bool permute_rows, bool permute_columns, bool invert_permutation) const{
    return (*this)->pmult(p,permute_rows,permute_columns,invert_permutation);
  }

  void CRSSparsity::spyMatlab(const std::string& mfile) const{
    (*this)->spyMatlab(mfile);
  }

  void CRSSparsity::spy(std::ostream &stream) const {
    for (int i=0;i<size1();++i) {
      for (int j=0;j<size2();++j) {
        stream << (getNZ(i,j)==-1? "." : "*");
      }
      stream << std::endl;
    }
  }

  bool CRSSparsity::isTranspose(const CRSSparsity& y) const{
    return (*this)->isTranspose(*static_cast<const CRSSparsityInternal*>(y.get()));
  }

  std::size_t CRSSparsity::hash() const{
    return (*this)->hash();
  }

  void CRSSparsity::assignCached(int nrow, int ncol, const std::vector<int>& col, const std::vector<int>& rowind){

    // Scalars and empty patterns are handled separately
    if(nrow==0 && ncol==0){
      // If empty    
      *this = getEmpty();
      return;
    } else if(nrow==1 && ncol==1){
      if(col.empty()){        
        // If sparse scalar
        *this = getScalarSparse();
        return;
      } else {
        // If dense scalar
        *this = getScalar();
        return;
      }
    }

    // Hash the pattern
    std::size_t h = hash_sparsity(nrow,ncol,col,rowind);

    // Get a reference to the cache
    CachingMap& cache = getCache();
    
    // Record the current number of buckets (for garbage collection below)
#ifdef USE_CXX11
    int bucket_count_before = cache.bucket_count();
#endif // USE_CXX11

    // WORKAROUND, functions do not appear to work when bucket_count==0
#ifdef USE_CXX11
    if(bucket_count_before>0){
#endif // USE_CXX11

      // Find the range of patterns equal to the key (normally only zero or one)
      pair<CachingMap::iterator,CachingMap::iterator> eq = cache.equal_range(h);

      // Loop over maching patterns
      for(CachingMap::iterator i=eq.first; i!=eq.second; ++i){
      
        // Get a weak reference to the cached sparsity pattern
        WeakRef& wref = i->second;
      
        // Check if the pattern still exists
        if(wref.alive()){
        
          // Get an owning reference to the cached pattern
          CRSSparsity ref = shared_cast<CRSSparsity>(wref.shared());
        
          // Check if the pattern matches
          if(ref.isEqual(nrow,ncol,col,rowind)){
          
            // Found match!
            assignNode(ref.get());
            return;

          } else {
            // There are two options, either the pattern has changed or there is a hash collision, so let's rehash the pattern
            std::size_t h_ref = ref.hash();
          
            if(h_ref!=h){ // The sparsity pattern has changed (the most likely event)

              // Create a new pattern
              assignNode(new CRSSparsityInternal(nrow, ncol, col, rowind));

              // Cache this pattern instead of the old one
              wref = *this;

              // Recache the old sparsity pattern 
              // TODO: recache "ref"
              return;

            } else { // There is a hash colision (unlikely, but possible)
              // Leave the pattern alone, continue to the next matching pattern
              continue; 
            }
          }
        } else {
          // Check if one of the other cache entries indeed has a matching sparsity
          CachingMap::iterator j=i;
          j++; // Start at the next matching key
          for(; j!=eq.second; ++j){
            if(j->second.alive()){
            
              // Recover cached sparsity
              CRSSparsity ref = shared_cast<CRSSparsity>(j->second.shared());
            
              // Match found if sparsity matches
              if(ref.isEqual(nrow,ncol,col,rowind)){
                assignNode(ref.get());
                return;
              }
            }
          }

          // The cached entry has been deleted, create a new one
          assignNode(new CRSSparsityInternal(nrow, ncol, col, rowind));
        
          // Cache this pattern
          wref = *this;

          // Return
          return;
        }
      }

      // END WORKAROUND
#ifdef USE_CXX11
    }
#endif // USE_CXX11

    // No matching sparsity pattern could be found, create a new one
    assignNode(new CRSSparsityInternal(nrow, ncol, col, rowind));

    // Cache this pattern
    //cache.insert(eq.second,std::pair<std::size_t,WeakRef>(h,ret));
    cache.insert(std::pair<std::size_t,WeakRef>(h,*this));

    // Garbage collection (currently only supported for unordered_multimap)
#ifdef USE_CXX11
    int bucket_count_after = cache.bucket_count();
    
    // We we increased the number of buckets, take time to garbage-collect deleted references
    if(bucket_count_before!=bucket_count_after){
      CachingMap::const_iterator i=cache.begin();
      while(i!=cache.end()){
        if(!i->second.alive()){
          i = cache.erase(i);
        } else {
          i++;
        }
      }
    }
#endif // USE_CXX11    
  }

  void CRSSparsity::clearCache(){
    getCache().clear();
  }


} // namespace CasADi
