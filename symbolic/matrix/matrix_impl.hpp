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

#ifndef MATRIX_IMPL_HPP
#define MATRIX_IMPL_HPP

// The declaration of the class is in a separate file
#include "matrix.hpp"
#include "matrix_tools.hpp"
#include "sparsity_tools.hpp"

namespace CasADi{
  // Implementations

  template<class T>
  const T& Matrix<T>::elem(int rr, int cc) const{
    int ind = sparsity().getNZ(rr,cc);
    if(ind==-1)
      return casadi_limits<T>::zero;
    else
      return at(ind);
  }

  template<class T>
  int Matrix<T>::stream_precision_ = 6;
  template<class T>
  int Matrix<T>::stream_width_ = 0;
  template<class T>
  bool Matrix<T>::stream_scientific_ = false;

  template<class T>
  T& Matrix<T>::elem(int rr, int cc){
    int oldsize = sparsity().size();
    int ind = sparsityRef().getNZ(rr,cc);
    if(oldsize != sparsity().size())
      data().insert(begin()+ind,T(0));
    return at(ind);
  }

  template<class T>
  bool Matrix<T>::__nonzero__() const {
    if (isNull()) {casadi_error("Cannot determine truth value of null Matrix.");}
    if (numel()!=1) {casadi_error("Only scalar Matrix could have a truth value, but you provided a shape" << dimString());}
    return CasADi::__nonzero__(at(0));
  }

  template<class T>
  const Matrix<T> Matrix<T>::sub(int rr, int cc) const{
    return elem(rr,cc);
  }

  template<class T>
  const Matrix<T> Matrix<T>::sub(const std::vector<int>& jj, const std::vector<int>& ii) const{
    // Nonzero mapping from submatrix to full
    std::vector<int> mapping;
  
    // Get the sparsity pattern - does bounds checking
    Sparsity sp = sparsity().sub(jj,ii,mapping);

    // Create return object
    Matrix<T> ret(sp);
  
    // Copy nonzeros
    for(int k=0; k<mapping.size(); ++k)
      ret.data()[k] = data()[mapping[k]];
  
    // Return (RVO)
    return ret;
  }

  template<class T>
  const Matrix<T> Matrix<T>::sub(const Matrix<int>& k, const std::vector<int>& ii) const{
    std::vector< int > rows = range(size1());
    std::vector< Matrix<T> > temp;
  
    if (!inBounds(ii,size2())) {
      casadi_error("Slicing [ii,k] out of bounds. Your ii contains " << *std::min_element(ii.begin(),ii.end()) << " up to " << *std::max_element(ii.begin(),ii.end()) << ", which is outside of the matrix shape " << dimString() << ".");
    }

    for (int i=0;i<ii.size();++i) {
      Matrix<T> m = k;
      for (int j=0;j<m.size();++j) {
        m.data()[j] = elem(k.at(j),ii.at(i));
      }
      temp.push_back(m);
    }
  
    return horzcat(temp);
  }

  template<class T>
  const Matrix<T> Matrix<T>::sub(const std::vector<int>& jj, const Matrix<int>& k) const{
    std::vector< int > cols = range(size2());
    std::vector< Matrix<T> > temp;

    if (!inBounds(jj,size1())) {
      casadi_error("Slicing [ii,k] out of bounds. Your jj contains " << *std::min_element(jj.begin(),jj.end()) << " up to " << *std::max_element(jj.begin(),jj.end()) << ", which is outside of the matrix shape " << dimString() << ".");
    }
  
    for (int j=0;j<jj.size();++j) {
      Matrix<T> m = k;
      for (int i=0;i<m.size();++i) {
        m.data()[i] = elem(jj.at(j),k.at(i));
      }
      temp.push_back(m);
    }
  
    return vertcat(temp);
  }

  template<class T>
  const Matrix<T> Matrix<T>::sub(const Matrix<int>& j, const Matrix<int>& i) const {
    casadi_assert_message(i.sparsity()==j.sparsity(),"sub(Imatrix i, Imatrix j): sparsities must match. Got " << i.dimString() << " and " << j.dimString() << ".");

    Matrix<T> ret(i.sparsity());
    for (int k=0;k<i.size();++k) {
      ret.data()[k] = elem(j.at(k),i.at(k));
    }

    return ret;
  }

  template<class T>
  const Matrix<T> Matrix<T>::sub(const Sparsity& sp, int dummy) const {
    casadi_assert_message(size1()==sp.size1() && size2()==sp.size2(),"sub(Sparsity sp): shape mismatch. This matrix has shape " << size1() << " x " << size2() << ", but supplied sparsity index has shape " << sp.size1() << " x " << sp.size2() << "." );
    Matrix<T> ret(sp);

    std::vector<unsigned char> mapping; // Mapping that will be filled by patternunion
    sparsity().patternCombine(sp, false, true, mapping);

    int k = 0;     // Flat index into non-zeros of this matrix
    int j = 0;     // Flat index into non-zeros of the resultant matrix;
    for (int i=0;i<mapping.size();++i) {
      if (mapping[i] & 1) { // If the original matrix has a non-zero entry in the union
        if (!(mapping[i] & 4)) ret.at(j) = at(k); // If this non-zero entry appears in the intersection, add it to the mapping
        k++;                 // Increment the original matrix' non-zero index counter
      }
      if (mapping[i] & 2) j++;
    }

    return ret;
  }

  template<class T>
  void Matrix<T>::setSub(const Matrix<T>& m, int j, int i){
    if(m.isDense()){
      elem(j,i) = m.toScalar();
    } else {
      setSub(m,std::vector<int>(1,j),std::vector<int>(1,i));
    }
  }

  template<class T>
  void Matrix<T>::setSub(const Matrix<T>& m, const std::vector<int>& rr, const std::vector<int>& cc){
    casadi_assert_message(m.numel()==1 || (cc.size() == m.size2() && rr.size() == m.size1()),"Dimension mismatch." << std::endl << "lhs is " << cc.size() << " x " << rr.size() << ", while rhs is " << m.dimString());

    if (!inBounds(rr,size1())) {
      casadi_error("setSub[.,rr,cc] out of bounds. Your rr contains " << *std::min_element(rr.begin(),rr.end()) << " up to " << *std::max_element(rr.begin(),rr.end()) << ", which is outside of the matrix shape " << dimString() << ".");
    }
    if (!inBounds(cc,size2())) {
      casadi_error("setSub [.,rr,cc] out of bounds. Your cc contains " << *std::min_element(cc.begin(),cc.end()) << " up to " << *std::max_element(cc.begin(),cc.end()) << ", which is outside of the matrix shape " << dimString() << ".");
    }
  
    // If m is scalar
    if(m.numel() != cc.size() * rr.size()){
      setSub(Matrix<T>::repmat(m.toScalar(),rr.size(),cc.size()),rr,cc);
      return;
    }

    if(isDense() && m.isDense()){
      // Dense mode
      for(int i=0; i<cc.size(); ++i) {
        for(int j=0; j<rr.size(); ++j) {
          data()[cc[i]*size1() + rr[j]]=m.data()[i*m.size1()+j];
        }
      }
    } else {
      // Sparse mode

      // Remove submatrix to be replaced
      erase(rr,cc);

      // Extend el to the same dimension as this
      Matrix<T> el_ext = m;
      el_ext.enlarge(size1(),size2(),rr,cc);

      // Unite the sparsity patterns
      *this = unite(*this,el_ext);
    }
  }

  template<class T>
  void Matrix<T>::setSub(const Matrix<T>& m, const std::vector<int>& jj, const Matrix<int>& i) {
    // If el is scalar
    if(m.isScalar() && (jj.size() > 1 || i.size() > 1)){
      setSub(repmat(Matrix<T>(i.sparsity(),m.toScalar()),jj.size(),1),jj,i);
      return;
    }

    if (!inBounds(jj,size1())) {
      casadi_error("setSub[.,i,jj] out of bounds. Your jj contains " << *std::min_element(jj.begin(),jj.end()) << " up to " << *std::max_element(jj.begin(),jj.end()) << ", which is outside of the matrix shape " << dimString() << ".");
    }
  
    Sparsity result_sparsity = repmat(i,jj.size(),1).sparsity();
  
  
    casadi_assert_message(result_sparsity == m.sparsity(),"setSub(Imatrix" << i.dimString() << ",Ivector(length=" << jj.size() << "),Matrix<T>)::Dimension mismatch. The sparsity of repmat(IMatrix," << jj.size() << ",1) = " << result_sparsity.dimString()  << " must match the sparsity of Matrix<T> = "  << m.dimString() << ".");
  
  
    std::vector<int> slice_i = range(i.size2());
  
    for(int k=0; k<jj.size(); ++k) {
      Matrix<T> el_k = m(range(k*i.size1(),(k+1)*i.size1()),slice_i);
      for (int j=0;j<i.size();++j) {
        elem(jj[k],i.at(j))=el_k.at(j);
      }
    }
  
  }

  template<class T>
  void Matrix<T>::setSub(const Matrix<T>& m, const Matrix<int>& j, const std::vector<int>& ii) {
  
    // If el is scalar
    if(m.isScalar() && (ii.size() > 1 || j.size() > 1)){
      setSub(repmat(Matrix<T>(j.sparsity(),m.toScalar()),1,ii.size()),j,ii);
      return;
    }

    if (!inBounds(ii,size2())) {
      casadi_error("setSub[.,ii,j] out of bounds. Your ii contains " << *std::min_element(ii.begin(),ii.end()) << " up to " << *std::max_element(ii.begin(),ii.end()) << ", which is outside of the matrix shape " << dimString() << ".");
    }
  
    Sparsity result_sparsity = repmat(j,1,ii.size()).sparsity();
  
  
    casadi_assert_message(result_sparsity == m.sparsity(),"setSub(Ivector(length=" << ii.size() << "),Imatrix" << j.dimString() << ",Matrix<T>)::Dimension mismatch. The sparsity of repmat(Imatrix,1," << ii.size() << ") = " << result_sparsity.dimString() << " must match the sparsity of Matrix<T> = " << m.dimString() << ".");
  
    std::vector<int> slice_j = range(j.size1());
  
    for(int k=0; k<ii.size(); ++k) {
      Matrix<T> el_k = m(slice_j,range(k*j.size2(),(k+1)*j.size2()));
      for (int i=0;i<j.size();++i) {
        elem(j.at(i),ii[k])=el_k.at(i);
      }
    }
  
  }


  template<class T>
  void Matrix<T>::setSub(const Matrix<T>& m, const Matrix<int>& j, const Matrix<int>& i) {
    casadi_assert_message(i.sparsity()==j.sparsity(),"setSub(., Imatrix i, Imatrix j): sparsities must match. Got " << i.dimString() << " for i and " << j.dimString() << " for j.");

    // If m is scalar
    if(m.isScalar() && i.numel() > 1){
      setSub(Matrix<T>(i.sparsity(),m.toScalar()),j,i);
      return;
    }
  
    casadi_assert_message(m.sparsity()==i.sparsity(),"setSub(Matrix m, Imatrix i, Imatrix j): sparsities must match. Got " << m.dimString() << " for m and " << j.dimString() << " for i and j.");
  
    for(int k=0; k<i.size(); ++k) {
      elem(j.at(k),i.at(k)) = m.at(k); 
    }
  }

  template<class T>
  void Matrix<T>::setSub(const Matrix<T>& m, const Sparsity& sp, int dummy) {
    casadi_assert_message(size2()==sp.size2() && size1()==sp.size1(),"sub(Sparsity sp): shape mismatch. This matrix has shape " << size2() << " x " << size1() << ", but supplied sparsity index has shape " << sp.size2() << " x " << sp.size1() << "." );
    // TODO: optimize this for speed
    Matrix<T> elm;
    if (m.isScalar()) {
      elm = Matrix<T>(sp,m.at(0));
    } else {
      elm = m.sub(sp);
    }

    for(int i=0; i<sp.colind().size()-1; ++i){
      for(int k=sp.colind()[i]; k<sp.colind()[i+1]; ++k){
        int j=sp.row()[k];
        elem(j,i)=elm.data()[k];
      }
    }
  }

  template<class T>
  const Matrix<T> Matrix<T>::getNZ(const std::vector<int>& k) const{
    try{
      Matrix<T> ret = zeros(k.size());
      for(int el=0; el<k.size(); ++el)
        ret.data()[el] = data().at(k[el]);
  
      return ret;
    } catch(std::out_of_range& ex){
      std::stringstream ss;
      ss << "Out of range error in Matrix<>::getNZ: " << k << " not all in range [0," << size() << ")";
      throw CasadiException(ss.str());
    }
  }

  template<class T>
  const Matrix<T> Matrix<T>::getNZ(const Matrix<int>& k) const{
    try{
      Matrix<T> ret = zeros(k.sparsity());
      for(int el=0; el<k.size(); ++el)
        ret.data()[el] = data().at(k.at(el));
  
      return ret;
    } catch(std::out_of_range& ex){
      std::stringstream ss;
      ss << "Out of range error in Matrix<>::getNZ: " << k << " not all in range [0," << size() << ")";
      throw CasadiException(ss.str());
    }
  }

  template<class T>
  void Matrix<T>::setNZ(int k, const Matrix<T>& m){
    if (k<0) k+=size();
    at(k) = m.toScalar();
  }

  template<class T>
  void Matrix<T>::setNZ(const std::vector<int>& kk, const Matrix<T>& m){
    if (m.isScalar()){
      // Assign all elements with the same scalar
      for(int k=0; k<kk.size(); ++k){
        setNZ(kk[k],m);
      }
    } else {
      // Assignment elementwise
      casadi_assert_message(kk.size()==m.size(),"Matrix<T>::setNZ: length of non-zero indices (" << kk.size() << ") " << std::endl << "must match size of rhs (" << m.size() << ").");
      for(int k=0; k<kk.size(); ++k){
        setNZ(kk[k],m[k]);
      }
    }
  }

  template<class T>
  void Matrix<T>::setNZ(const Matrix<int>& kk, const Matrix<T>& m){
    if (m.isScalar()){
      // Assign all elements with the same scalar
      for(int k=0; k<kk.size(); ++k){
        setNZ(kk.at(k),m);
      }
    } else if (kk.isDense() && !m.isDense() && kk.size2()==m.size2() && kk.size1()==m.size1()) {
      const std::vector<int>& row = m.sparsity().row();
      const std::vector<int>& colind = m.sparsity().colind();
      for(int i=0; i<colind.size()-1; ++i){
        for(int k=colind[i]; k<colind[i+1]; ++k){
          int j=row[k];
          setNZ(kk.elem(j,i),m[k]);
        }
      }
    } else {
      casadi_assert_message(kk.sparsity()==m.sparsity(),"Matrix<T>::setNZ: sparsity of IMatrix index " << kk.dimString() << " " << std::endl << "must match sparsity of rhs " << m.dimString() << ".");
      for(int k=0; k<kk.size(); ++k){
        setNZ(kk.at(k),m[k]);
      }
    }
  }

  template<class T>
  void Matrix<T>::makeDense(int nrow, int ncol, const T& val){
    // Quick return if already dense
    if(ncol*nrow == size())
      return;
  
    if(size2()!=ncol || size1()!=nrow){
      // Also resize
      sparsity_ = Sparsity(nrow,ncol,true);
      std::fill(data().begin(),data().end(),val);
      data().resize(ncol*nrow, val);
    } else {
      // Create a new data vector
      data().resize(ncol*nrow,val);
    
      // Loop over the cols in reverse order
      for(int i=ncol-1; i>=0; --i){
        // Loop over nonzero elements in reverse order
        for(int el=colind(i+1)-1; el>=colind(i); --el){
          // Row
          int j = row(el);
        
          // Swap the old position with the new position
          if(el!=j+i*nrow){
            data()[j+i*nrow] = data()[el];
            data()[el] = val;
          }
        }
      }
      
      // Save the new sparsity pattern
      sparsity_ = Sparsity(nrow,ncol,true);
    }
  }

  template<class T>
  Matrix<T>::Matrix() : sparsity_(Sparsity(0,0,false)){
  }

  template<class T>
  Matrix<T>::Matrix(const Matrix<T>& m) : sparsity_(m.sparsity_), data_(m.data_){
  }

  template<class T>
  Matrix<T>::Matrix(const std::vector<T>& x) : sparsity_(Sparsity(x.size(),1,true)), data_(x){
  }

  template<class T>
  Matrix<T>::Matrix(const std::vector<T>& x, int nrow, int ncol) : sparsity_(Sparsity(nrow,ncol,true)), data_(x){
    casadi_assert_message(x.size() == nrow*ncol, "Dimension mismatch." << std::endl << "You supplied a vector of length " << x.size() << ", but " << nrow << " x " << ncol << " = " << nrow*ncol);
  }

  template<class T>
  Matrix<T>& Matrix<T>::operator=(const Matrix<T>& m){
    sparsity_ = m.sparsity_;
    data_ = m.data_;
    return *this;
  }

  template<class T>
  Matrix<T>::Matrix(int nrow, int ncol) : sparsity_(Sparsity(nrow,ncol,false)){
  }

  template<class T>
  Matrix<T>::Matrix(int nrow, int ncol, const T& val) : sparsity_(Sparsity(nrow,ncol,true)), data_(std::vector<T>(nrow*ncol, val)){
  }

  template<class T>
  void Matrix<T>::makeEmpty(int nrow, int ncol){
    sparsity_ = Sparsity(nrow,ncol,false);
    data().clear();
  }

  template<class T>
  std::string Matrix<T>::className(){ return matrixName<T>(); }

  template<class T>
  void Matrix<T>::printScalar(std::ostream &stream) const {
    casadi_assert_message(numel()==1, "Not a scalar");
  
    std::streamsize precision = stream.precision();
    std::streamsize width = stream.width();
    std::ios_base::fmtflags flags = stream.flags();
  
    stream.precision(stream_precision_);
    stream.width(stream_width_);
    if (stream_scientific_) {
      stream.setf(std::ios::scientific);
    } else {
      stream.unsetf(std::ios::scientific);
    }
    
    if(size()==0){
      stream << "00";
    } else {
      stream << toScalar();
    }
  
    stream << std::flush;
    stream.precision(precision);
    stream.width(width);
    stream.flags(flags); 
  }
  
  template<class T>
  void Matrix<T>::printVector(std::ostream &stream) const {
    casadi_assert_message(isVector(),"Not a vector");
  
    std::streamsize precision = stream.precision();
    std::streamsize width = stream.width();
    std::ios_base::fmtflags flags = stream.flags();
  
    stream.precision(stream_precision_);
    stream.width(stream_width_);
    if (stream_scientific_) {
      stream.setf(std::ios::scientific);
    } else {
      stream.unsetf(std::ios::scientific);
    }
  
    // Access data structures
    const std::vector<int>& r = row();
      
    // Nonzero
    int el=0;

    // Loop over rows
    stream << "[";
    for(int rr=0; rr<size1(); ++rr){
      // Add delimitor
      if(rr!=0) stream << ",";
      
      // Check if nonzero
      if(el<r.size() && rr==r[el]){
        stream << at(el++);
      } else {
        stream << "00";
      }
    }
    stream << "]"; 
    
    stream << std::flush;
    stream.precision(precision);
    stream.width(width);
    stream.flags(flags); 
  }

  template<class T>
  void Matrix<T>::printDense(std::ostream &stream) const{
    // Print as a single line
    bool oneliner=this->size1()<=1;
  
    std::streamsize precision = stream.precision();
    std::streamsize width = stream.width();
    std::ios_base::fmtflags flags = stream.flags();
  
    stream.precision(stream_precision_);
    stream.width(stream_width_);
    if (stream_scientific_) {
      stream.setf(std::ios::scientific);
    } else {
      stream.unsetf(std::ios::scientific);
    }

    // Index counter for each column
    std::vector<int> cind = colind();

    // Loop over rows
    for(int rr=0; rr<size1(); ++rr){
      // Beginning of row
      if(rr==0){
        if(!oneliner) stream << std::endl;
        stream << "[[";
      } else {
        stream << " [";
      }
      
      // Loop over columns
      for(int cc=0; cc<size2(); ++cc){
        // Separating comma
        if(cc>0) stream << ",";

        // Check if nonzero
        if(cind[cc]<colind(cc+1) && row(cind[cc])==rr){
          stream << data().at(cind[cc]++);
        } else {
          stream << "00";
        }
      }
    
      // End of row
      if(rr<size1()-1){
        stream << "],";
        if(!oneliner) stream << std::endl;
      } else {
        stream << "]]";
      }
    }
        
    stream << std::flush;
    stream.precision(precision);
    stream.width(width);
    stream.flags(flags);
  }

  template<class T>
  void Matrix<T>::printSparse(std::ostream &stream) const {
    if(size()==0){
      stream << "all zero sparse: " << size1() << "-by-" << size2();
    } else {
      stream << "sparse: " << size1() << "-by-" << size2() << ", " << size() << " nnz" << std::endl;
      for(int cc=0; cc<size2(); ++cc){
        for(int el=colind(cc); el<colind(cc+1); ++el){
          int rr=row(el);
          stream << " (" << rr << "," << cc << ") -> " << at(el) << std::endl;
        }
      }
    }
    stream << std::flush;
  }

  template<class T>
  void Matrix<T>::print(std::ostream &stream) const{
    if(isEmpty()){
      stream << "[]";
    } else if(numel()==1){
      printScalar(stream);
    } else if(isVector()){
      printVector(stream);
    } else if(std::max(size1(),size2())<=10 || double(size())/numel()>=0.5){ // if "small" or "dense"
      printDense(stream);
    } else {
      printSparse(stream);
    }
  }

  template<class T>
  void Matrix<T>::repr(std::ostream &stream) const{
    stream << className() << "(";
    print(stream);
    stream << ")" << std::flush;
  }

  template<class T>
  const std::vector<int>& Matrix<T>::row() const{
    return sparsity().row();
  }

  template<class T>
  const std::vector<int>& Matrix<T>::colind() const{
    return sparsity_.colind();
  }

  template<class T>
  int Matrix<T>::row(int el) const{
    return sparsity_.row(el);
  }

  template<class T>
  int Matrix<T>::colind(int col) const{
    return sparsity_.colind(col);
  }

  template<class T>
  void Matrix<T>::reserve(int nnz){
    reserve(nnz,size2());
  }

  template<class T>
  void Matrix<T>::reserve(int nnz, int ncol){
    data().reserve(nnz);
    sparsity_.reserve(nnz,ncol);
  }

  template<class T>
  void Matrix<T>::resize(int nrow, int ncol){
    sparsity_.resize(nrow,ncol);
  }

  template<class T>
  void Matrix<T>::clear(){
    sparsity_ = Sparsity(0,0,false);
    data().clear();
  }

  template<class T>
  Matrix<T>::Matrix(double val) : sparsity_(Sparsity(1,1,true)), data_(std::vector<T>(1,val)) {
  }

  template<class T>
  Matrix<T>::Matrix(int nrow, int ncol, const std::vector<int>& colind, const std::vector<int>& row, const std::vector<T>& d) : sparsity_(Sparsity(nrow,ncol,colind,row)), data_(d){
    if(data_.size() != sparsity_.size())
      data_.resize(sparsity_.size()); // Why not throw an error?
    sanityCheck(true);
  }

  template<class T>
  Matrix<T>::Matrix(const std::vector< std::vector<T> >& d){
    // Get dimensions
    int nrow=d.size();
    int ncol=d.empty() ? 1 : d.front().size();

    // Assert consistency
    for(int rr=0; rr<nrow; ++rr){
      casadi_assert_message(ncol==d[rr].size(), 
        "Matrix<T>::Matrix(const std::vector< std::vector<T> >& d): shape mismatch" << std::endl <<
        "Attempting to construct a matrix from a nested list." << std::endl <<
        "I got convinced that the desired size is ("<< nrow << " x " << ncol << " ), but now I encounter a vector of size (" << 
        d[rr].size() <<  " )" << std::endl);
    }

    // Form matrix
    sparsity_ = Sparsity(nrow,ncol,true);
    data().resize(nrow*ncol);
    typename std::vector<T>::iterator it=begin();
    for(int cc=0; cc<ncol; ++cc){
      for(int rr=0; rr<nrow; ++rr){
        *it++ = d[rr][cc];
      }
    }
  }

  template<class T>
  Matrix<T>::Matrix(const Sparsity& sparsity, const T& val) : sparsity_(sparsity), data_(std::vector<T>(sparsity.size(),val)){
  }

  template<class T>
  Matrix<T>::Matrix(const Sparsity& sparsity, const std::vector<T>& d) : sparsity_(sparsity), data_(d) {
    casadi_assert_message(sparsity.size()==d.size(),"Size mismatch." << std::endl << "You supplied a sparsity of " << sparsity.dimString() << ", but the supplied vector is of length " << d.size());
  }

  template<class T>
  void Matrix<T>::setZero(){
    setAll(0);
  }

  template<class T>
  void Matrix<T>::setAll(const T& val){
    std::fill(begin(),end(),val);
  }

  template<class T>
  Matrix<T> Matrix<T>::unary(int op, const Matrix<T> &x){
    // Return value
    Matrix<T> ret(x.sparsity());
  
    // Nonzeros
    std::vector<T>& ret_data = ret.data();
    const std::vector<T>& x_data = x.data();
  
    // Do the operation on all non-zero elements
    for(int el=0; el<x.size(); ++el){
      casadi_math<T>::fun(op,x_data[el],x_data[el],ret_data[el]);
    }

    // Check the value of the structural zero-entries, if there are any
    if(!x.isDense() && !operation_checker<F0XChecker>(op)){
      // Get the value for the structural zeros
      T fcn_0;
      casadi_math<T>::fun(op,0,0,fcn_0);
      if(!casadi_limits<T>::isZero(fcn_0)){ // Remove this if?
        ret.makeDense(ret.size1(),ret.size2(),fcn_0);
      }
    }
    
    return ret;
  }

  template<class T>
  Matrix<T> Matrix<T>::operator-() const{
    return unary(OP_NEG,*this);
  }

  template<class T>
  Matrix<T> Matrix<T>::operator+() const{
    return *this;
  }

  template<class T>
  Matrix<T> Matrix<T>::__add__(const Matrix<T> &y) const{
    return binary(OP_ADD,*this,y);
  }

  template<class T>
  Matrix<T> Matrix<T>::__sub__(const Matrix<T> &y) const{
    return binary(OP_SUB,*this,y);
  }

  template<class T>
  Matrix<T> Matrix<T>::__mul__(const Matrix<T> &y) const{
    return binary(OP_MUL,*this,y);
  }

  template<class T>
  Matrix<T> Matrix<T>::__div__(const Matrix<T> &y) const{
    return binary(OP_DIV,*this,y);
  }

  template<class T>
  Matrix<T> Matrix<T>::__lt__(const Matrix<T> &y) const{
    return binary(OP_LT,*this,y);
  }

  template<class T>
  Matrix<T> Matrix<T>::__le__(const Matrix<T> &y) const{
    return binary(OP_LE,*this,y);
  }

  template<class T>
  Matrix<T> Matrix<T>::__eq__(const Matrix<T> &y) const{
    return binary(OP_EQ,*this,y);
  }

  template<class T>
  Matrix<T> Matrix<T>::__ne__(const Matrix<T> &y) const{
    return binary(OP_NE,*this,y);
  }

  template<class T>
  Matrix<T> Matrix<T>::__mrdivide__(const Matrix<T>& b) const { if (b.numel()==1) return *this/b; throw CasadiException("mrdivide: Not implemented");}

  template<class T>
  Matrix<T> Matrix<T>::__mpower__(const Matrix<T>& b) const { if (b.numel()==1) return (*this).__pow__(b); throw CasadiException("mpower: Not implemented");}

  template<class T>
  Sparsity& Matrix<T>::sparsityRef(){
    sparsity_.makeUnique(); // NOTE: Remove?
    return sparsity_;
  }

  template<class T>
  void Matrix<T>::getBand(int kl, int ku, int ldres, T *res) const{
    // delete the content of the matrix
    for(int j=0; j<size1(); ++j) // loop over rows
      for(int s=0; s<kl+ku+1; ++s) // loop over the subdiagonals
        res[s + ldres*j] = 0;
  
    // loop over cols
    for(int i=0; i<size2(); ++i){ 
    
      // loop over the non-zero elements
      for(int el=colind(i); el<colind(i+1); ++el){ 
        int j=row(el);  // row
      
        // Check if we have not yet inside the band
        if(j<i-kl) continue;

        // Check if we are already outside the band
        if(j>i+ku) break;

        // Get the subdiagonal
        int s = i - j + ku;

        // Store the element
        res[s + ldres*j] = data()[el];
      }
    }
  }

  template<class T>
  void Matrix<T>::set(T val, SparsityType sp){
    std::fill(data().begin(),data().end(),val);
  }
    
  template<class T>
  void Matrix<T>::get(T& val, SparsityType sp) const{
    getArray(&val,1,DENSE);
  }

  template<class T>
  void Matrix<T>::set(const std::vector<T>& val, SparsityType sp){
    setArray(val.empty() ? 0 : &val.front(),val.size(),sp);
  }

  template<class T>
  void Matrix<T>::get(std::vector<T>& val, SparsityType sp) const{
    getArray(val.empty() ? 0 : &val.front(),val.size(),sp);
  }

  template<class T>
  void Matrix<T>::set(const Matrix<T>& val, SparsityType sp){
    sparsity().set(getPtr(data()),getPtr(val.data()),val.sparsity());
  }

  template<class T>
  void Matrix<T>::setBV(const Matrix<T>& val){
    bvec_t* bw_this = reinterpret_cast<bvec_t*>(getPtr(data()));
    const bvec_t* bw_val = reinterpret_cast<const bvec_t*>(getPtr(val.data()));
    sparsity().set(bw_this,bw_val,val.sparsity());
  }

  template<class T>
  void Matrix<T>::setZeroBV(){
    bvec_t* bw_this = reinterpret_cast<bvec_t*>(getPtr(data()));
    std::fill(bw_this,bw_this+size(),bvec_t(0));
  }

  template<class T>
  void Matrix<T>::borBV(const Matrix<T>& val){
    bvec_t* bw_this = reinterpret_cast<bvec_t*>(getPtr(data()));
    const bvec_t* bw_val = reinterpret_cast<const bvec_t*>(getPtr(val.data()));
    sparsity().bor(bw_this,bw_val,val.sparsity());
  }

  template<class T>
  void Matrix<T>::getArrayBV(bvec_t* val, int len) const{
    casadi_assert(len==size());
    const bvec_t* bw_this = reinterpret_cast<const bvec_t*>(getPtr(data()));
    std::copy(bw_this,bw_this+len,val);
  }

  template<class T>
  void Matrix<T>::setArrayBV(const bvec_t* val, int len){
    casadi_assert(len==size());
    bvec_t* bw_this = reinterpret_cast<bvec_t*>(getPtr(data()));
    std::copy(val,val+len,bw_this);
  }

  template<class T>
  void Matrix<T>::borArrayBV(const bvec_t* val, int len){
    casadi_assert(len==size());
    bvec_t* bw_this = reinterpret_cast<bvec_t*>(getPtr(data()));
    for(int i=0; i<len; ++i) *bw_this++ |= *val++;
  }

  template<class T>
  void Matrix<T>::get(Matrix<T>& val, SparsityType sp) const{
    val.set(*this,sp);
  }

  template<class T>
  void Matrix<T>::set(const T* val, SparsityType sp){
    int len = sp==SPARSE ? size() : sp==DENSE || sp==DENSETRANS ? numel() : sp==SPARSESYM ? sizeU() : -1;
    setArray(val,len,sp);
  }

  template<class T>
  void Matrix<T>::get(T* val, SparsityType sp) const{
    int len = sp==SPARSE ? size() : sp==DENSE || sp==DENSETRANS ? numel() : sp==SPARSESYM ? sizeU() : -1;
    getArray(val,len,sp);
  }

  template<class T>
  void Matrix<T>::getArray(T* val, int len, SparsityType sp) const{
    // Get references to data for quick access
    const std::vector<T> &data = this->data();
    const int size1 = this->size1();
    const int size2 = this->size2();
    const std::vector<int>& colind = this->colind();
    const std::vector<int>& row = this->row();
    
    if(sp==SPARSE || (sp==DENSE && isDense())){
      casadi_assert_message(len==size(),
                            "Matrix<T>::getArray: Dimension mismatch." << std::endl <<  
                            "Trying to fetch " << len << " elements from a " << dimString() << " matrix with " << size() << " non-zeros.");
      copy(data.begin(),data.end(),val);
    } else if(sp==DENSE){
      casadi_assert_message(len==numel(),
                            "Matrix<T>::getArray: Dimension mismatch." << std::endl <<  
                            "Trying to fetch " << len << " elements from a " << dimString() << " matrix with " << numel() << " entries.");
      // Begin with all zeros
      std::fill(val,val+len,0);
     
      // Add the nonzeros
      for(int cc=0; cc<size2; ++cc){ // loop over columns
        for(int el=colind[cc]; el<colind[cc+1]; ++el){ // loop over the non-zero elements
          int rr=row[el];
          val[rr+cc*size1] = data[el];
        }
      }
    } else if(sp==DENSETRANS){
      casadi_assert_message(len==numel(),
                            "Matrix<T>::getArray: Dimension mismatch." << std::endl <<  
                            "Trying to fetch " << len << " elements from a " << dimString() << " matrix with " << numel() << " entries.");
      // Begin with all zeros
      std::fill(val,val+len,0);
     
      // Add the nonzeros
      for(int cc=0; cc<size2; ++cc){ // loop over columns
        for(int el=colind[cc]; el<colind[cc+1]; ++el){ // loop over the non-zero elements
          int rr=row[el];
          val[cc+rr*size2] = data[el];
        }
      }
    } else if(sp==SPARSESYM){
      // copy to the result vector
      int nz = 0;
      for(int cc=0; cc<size2; ++cc){
        // Loop over the elements in the col
        for(int el=colind[cc]; el<colind[cc+1]; ++el){ // loop over the non-zero elements
          if(row[el] > cc) break; // break inner loop (only upper triangular part is used)
          val[nz++] = data[el];
        }
      }
    } else {
      casadi_error("Matrix<T>::getArray: not SPARSE, SPARSESYM, DENSE or DENSETRANS");
    }
  }

  /**
     Set stride to zero for unstrided acces
  */
  template<class T>
  void Matrix<T>::getStridedArray(T* val, int len, int stride1, int stride2, SparsityType sp) const{
    if (stride1==0 || stride2==0 || (sp==DENSE && stride2==1 && stride1==size1())) return getArray(val, len, sp);

    // Get references to data for quick access
    const std::vector<T> &data = this->data();
    const int size1 = this->size1();
    const int size2 = this->size2();
    const std::vector<int>& colind = this->colind();
    const std::vector<int>& row = this->row();
    
    if(sp==SPARSE){
      throw CasadiException("Matrix<T>::getArray: strided SPARSE not implemented");
    } else if(sp==DENSE && isDense()) {
      for(int cc=0; cc<size2; ++cc){ // loop over columns
        for(int el=colind[cc]; el<colind[cc+1]; ++el){ // loop over the non-zero elements
          int rr=row[el];
          val[rr*stride2 + cc*stride1] = data[el];
        }
      }
    } else if(sp==DENSETRANS && isDense()) {
      for(int cc=0; cc<size2; ++cc){ // loop over columns
        for(int el=colind[cc]; el<colind[cc+1]; ++el){ // loop over the non-zero elements
          int rr=row[el];
          val[cc*stride2 + rr*stride1] = data[el];
        }
      }
    } else if(sp==DENSE){
      throw CasadiException("Matrix<T>::getStridedArray: strided sparse DMatrix->dense not implemented");
    } else if(sp==SPARSESYM){
      throw CasadiException("Matrix<T>::getStridedArray: strided SPARSESYM not implemented");
    } else {
      throw CasadiException("Matrix<T>::getStridedArray: not SPARSE or DENSE");
    }

  }

  template<class T>
  void Matrix<T>::setArray(const T* val, int len, SparsityType sp){
    // Get references to data for quick access
    std::vector<T> &data = this->data();
    const int size1 = this->size1();
    const int size2 = this->size2();
    const std::vector<int>& colind = this->colind();
    const std::vector<int>& row = this->row();

    if(sp==SPARSE || (sp==DENSE && numel()==size())){
      casadi_assert_message(len==size(),
                            "Matrix<T>::setArray: Dimension mismatch." << std::endl <<  
                            "Trying to pass " << len << " elements to a " << dimString() << " matrix with " << size() << " non-zeros.");
      copy(val,val+len,data.begin());
    } else if(sp==DENSE){
      casadi_assert_message(len==numel(),
                            "Matrix<T>::setArray: Dimension mismatch." << std::endl <<  
                            "Trying to pass " << len << " elements to a " << dimString() << " matrix with " << numel() << " entries.");
      // Get the nonzeros
      for(int cc=0; cc<size2; ++cc){ // loop over columns
        for(int el=colind[cc]; el<colind[cc+1]; ++el){ // loop over the non-zero elements
          int rr=row[el];
          data[el] = val[rr+cc*size1];
        }
      }
    } else if(sp==DENSETRANS){
      casadi_assert_message(len==numel(),
                            "Matrix<T>::setArray: Dimension mismatch." << std::endl <<  
                            "Trying to pass " << len << " elements to a " << dimString() << " matrix with " << numel() << " entries.");
      // Get the nonzeros
      for(int cc=0; cc<size2; ++cc){ // loop over columns
        for(int el=colind[cc]; el<colind[cc+1]; ++el){ // loop over the non-zero elements
          int rr=row[el];
          data[el] = val[cc+rr*size2];
        }
      }
    } else if(sp==SPARSESYM) {
      // NOTE: Has to be rewritten! sparsity().transpose(...) involves memory allocation and is not threadsafe!!!
      // This routines is supposed to be used inside threadsafe code.
      std::vector<int> mapping;
      sparsity().transpose(mapping,false);
      // copy to the result vector
      int nz = 0;
      for(int cc=0; cc<size2; ++cc){
        // Loop over the elements in the col
        for(int el=colind[cc]; el<colind[cc+1]; ++el){ // loop over the non-zero elements
          if(row[el] > cc) break; // break inner loop (only lower triangular part is used)
          data[mapping[el]] = data[el] = val[nz++];
        }
      }
    } else {
      throw CasadiException("Matrix<T>::setArray: not SPARSE, SPARSESYM, DENSE or DENSETRANS");
    }
  }

  template<class T>
  void Matrix<T>::getArray(T* val) const{
    getArray(val,size(),SPARSE);
  }

  template<class T>
  void Matrix<T>::setArray(const T* val){
    setArray(val,size(),SPARSE);
  }

  template<class T>
  Matrix<T> Matrix<T>::__pow__(const Matrix<T>& y) const{
    return binary(OP_POW,*this,y);
  }

  template<class T>
  Matrix<T> Matrix<T>::__constpow__(const Matrix<T>& y) const{
    return binary(OP_CONSTPOW,*this,y);
  }

  template<class T>
  Matrix<T> Matrix<T>::sin() const{
    return unary(OP_SIN,*this);
  }

  template<class T>
  Matrix<T> Matrix<T>::cos() const{
    return unary(OP_COS,*this);
  }

  template<class T>
  Matrix<T> Matrix<T>::tan() const{
    return unary(OP_TAN,*this);
  }

  template<class T>
  Matrix<T> Matrix<T>::erf() const{
    return unary(OP_ERF,*this);
  }

  template<class T>
  Matrix<T> Matrix<T>::arcsin() const{
    return unary(OP_ASIN,*this);
  }

  template<class T>
  Matrix<T> Matrix<T>::arccos() const{
    return unary(OP_ACOS,*this);
  }

  template<class T>
  Matrix<T> Matrix<T>::arctan() const{
    return unary(OP_ATAN,*this);
  }

  template<class T>
  Matrix<T> Matrix<T>::sinh() const{
    return unary(OP_SINH,*this);
  }

  template<class T>
  Matrix<T> Matrix<T>::cosh() const{
    return unary(OP_COSH,*this);
  }

  template<class T>
  Matrix<T> Matrix<T>::tanh() const{
    return unary(OP_TANH,*this);
  }

  template<class T>
  Matrix<T> Matrix<T>::arcsinh() const{
    return unary(OP_ASINH,*this);
  }

  template<class T>
  Matrix<T> Matrix<T>::arccosh() const{
    return unary(OP_ACOSH,*this);
  }

  template<class T>
  Matrix<T> Matrix<T>::arctanh() const{
    return unary(OP_ATANH,*this);
  }

  template<class T>
  Matrix<T> Matrix<T>::exp() const{
    return unary(OP_EXP,*this);
  }

  template<class T>
  Matrix<T> Matrix<T>::log() const{
    return unary(OP_LOG,*this);
  }

  template<class T>
  Matrix<T> Matrix<T>::log10() const{
    return log()*(1/std::log(10.));
  }

  template<class T>
  Matrix<T> Matrix<T>::sqrt() const{
    return unary(OP_SQRT,*this);
  }

  template<class T>
  Matrix<T> Matrix<T>::floor() const{
    return unary(OP_FLOOR,*this);
  }

  template<class T>
  Matrix<T> Matrix<T>::ceil() const{
    return unary(OP_CEIL,*this);
  }

  template<class T>
  Matrix<T> Matrix<T>::fabs() const{
    return unary(OP_FABS,*this);
  }

  template<class T>
  Matrix<T> Matrix<T>::sign() const{
    return unary(OP_SIGN,*this);
  }

  template<class T>
  Matrix<T> Matrix<T>::__copysign__(const Matrix<T>& y) const{
    return binary(OP_COPYSIGN,*this,y);
  }

  template<class T>
  Matrix<T> Matrix<T>::erfinv() const{
    return unary(OP_ERFINV,*this);
  }

  template<class T>
  Matrix<T> Matrix<T>::fmin(const Matrix<T>& y) const{
    return binary(OP_FMIN,*this,y);
  }

  template<class T>
  Matrix<T> Matrix<T>::arctan2(const Matrix<T>& y) const{
    return binary(OP_ATAN2,*this,y);
  }

  template<class T>
  Matrix<T> Matrix<T>::fmax(const Matrix<T>& y) const{
    return binary(OP_FMAX,*this,y);
  }

  template<class T>
  Matrix<T> Matrix<T>::printme(const Matrix<T>& y) const{
    return binary(OP_PRINTME,*this,y);
  }

  template<class T>
  Matrix<T> Matrix<T>::logic_not() const{
    return unary(OP_NOT,*this);
  }

  template<class T>
  Matrix<T> Matrix<T>::logic_and(const Matrix<T>& y) const{
    return binary(OP_AND,*this,y);
  }

  template<class T>
  Matrix<T> Matrix<T>::logic_or(const Matrix<T>& y) const{
    return binary(OP_OR,*this,y);
  }

  template<class T>
  Matrix<T> Matrix<T>::if_else_zero(const Matrix<T>& y) const{
    return binary(OP_IF_ELSE_ZERO,*this,y);
  }

  template<class T>
  std::vector<T>& Matrix<T>::data(){
    return data_;  
  }
    
  template<class T>
  const std::vector<T>& Matrix<T>::data() const{
    return data_;  
  }

  template<class T>
  void Matrix<T>::erase(const std::vector<int>& rr, const std::vector<int>& cc){
    // Erase from sparsity pattern
    std::vector<int> mapping = sparsityRef().erase(rr,cc);
  
    // Update non-zero entries
    for(int k=0; k<mapping.size(); ++k)
      data()[k] = data()[mapping[k]];
    
    // Truncate nonzero vector
    data().resize(mapping.size());
  }


  template<class T>
  void Matrix<T>::remove(const std::vector<int>& rr, const std::vector<int>& cc) {
    if (!inBounds(rr,size1())) {
      casadi_error("Remove(rr,cc) out of bounds. Your rr contains " << *std::min_element(rr.begin(),rr.end()) << " up to " << *std::max_element(rr.begin(),rr.end()) << ", which is outside of the matrix shape " << dimString() << ".");
    }
    if (!inBounds(cc,size2())) {
      casadi_error("Remove(rr,cc) out of bounds. Your cc contains " << *std::min_element(cc.begin(),cc.end()) << " up to " << *std::max_element(cc.begin(),cc.end()) << ", which is outside of the matrix shape " << dimString() << ".");
    }
  
    // Remove by performing a complementary slice
    std::vector<int> rrc = complement(rr,size1());
    std::vector<int> ccc = complement(cc,size2());
  
    Matrix<T> ret = (*this)(rrc,ccc);
  
    operator=(ret);
  
  }

  template<class T>
  void Matrix<T>::enlarge(int nrow, int ncol, const std::vector<int>& rr, const std::vector<int>& cc){
    sparsityRef().enlarge(nrow,ncol,rr,cc);
  }

  template<class T>
  void Matrix<T>::sanityCheck(bool complete) const {
    sparsity_.sanityCheck(complete);
  
    if (data_.size()!=sparsity_.row().size()) {
      std::stringstream s;
      s << "Matrix:Compressed Col Storage is not sane. The following must hold:" << std::endl;
      s << "  data.size() = nrow, but got   row.size()  = " << data_.size() << "   and   nrow = "  << sparsity_.row().size() << std::endl;
      s << "  Note that the signature is as follows: DMatrix (ncol, nrow, row, colind, data)." << std::endl;
      casadi_error(s.str());
    }
  }

  template<class T>
  Matrix<T> Matrix<T>::mul(const Matrix<T> &y, const Sparsity& sp_z) const {
    return this->mul_smart(y, sp_z);
  }

  template<class T>
  Matrix<T> Matrix<T>::mul_full(const Matrix<T> &y, const Sparsity& sp_z) const{
    // First factor
    const Matrix<T>& x = *this;
  
    // Return object (assure RVO)
    Matrix<T> ret;

    // Matrix multiplication

    // Form the transpose of x
    Matrix<T> x_trans = x.trans();
  
    if (sp_z.isNull()) {
      // Create the sparsity pattern for the matrix-matrix product
      Sparsity spres = y.sparsity().patternProduct(x_trans.sparsity());

      // Create the return object
      ret = Matrix<T>::zeros(spres);
    } else {
      ret = Matrix<T>::zeros(sp_z);
    }

    // Carry out the matrix product
    mul_no_alloc_tn(x_trans,y,ret);
  
    return ret;
  }

  template<class T>
  void Matrix<T>::mul_no_alloc_nn(const Matrix<T> &x, const Matrix<T> &y, Matrix<T>& z){
    // Assert dimensions
    casadi_assert_message(x.size1()==z.size1(),"Dimension error. Got x=" << x.dimString() << " and z=" << z.dimString() << ".");
    casadi_assert_message(y.size2()==z.size2(),"Dimension error. Got y=" << y.dimString() << " and z=" << z.dimString() << ".");
    casadi_assert_message(y.size1()==x.size2(),"Dimension error. Got y=" << y.dimString() << " and x=" << x.dimString() << ".");
  
    // Direct access to the arrays
    const std::vector<int> &y_colind = y.colind();
    const std::vector<int> &y_row = y.row();
    const std::vector<T> &y_data = y.data();
    const std::vector<int> &x_colind = x.colind();
    const std::vector<int> &x_row = x.row();
    const std::vector<T> &x_data = x.data();
    const std::vector<int> &z_colind = z.colind();
    const std::vector<int> &z_row = z.row();
    std::vector<T> &z_data = z.data();

    // loop over the cols of the first argument
    for(int i=0; i<y_colind.size()-1; ++i){
      for(int el=y_colind[i]; el<y_colind[i+1]; ++el){ // loop over the non-zeros of the first argument
        int j = y_row[el];
        int el1 = z_colind[i];
        int el2 = x_colind[j];
        while(el1 < z_colind[i+1] && el2 < x_colind[j+1]){ // loop over matching non-zero elements
          int j1 = z_row[el1];
          int i2 = x_row[el2];      
          if(j1==i2){
            z_data[el1++] += y_data[el]*x_data[el2++];
          } else if(j1<i2) {
            el1++;
          } else {
            el2++;
          }
        }
      }
    }
  }

  template<class T>
  void Matrix<T>::mul_no_alloc_tn(const Matrix<T> &x_trans, const std::vector<T> &y, std::vector<T>& z){
    // Assert dimensions
    casadi_assert_message(x_trans.size2()==z.size(),"Dimension error. Got x_trans=" << x_trans.dimString() << " and z=" << z.size() << ".");
    casadi_assert_message(x_trans.size1()==y.size(),"Dimension error. Got x_trans=" << x_trans.dimString() << " and y=" << y.size() << ".");
    
    // Direct access to the arrays
    const std::vector<int> &x_rowind = x_trans.colind();
    const std::vector<int> &x_col = x_trans.row();
    const std::vector<T> &x_trans_data = x_trans.data();
    
    // loop over the columns of the matrix
    for(int i=0; i<x_rowind.size()-1; ++i){
      for(int el=x_rowind[i]; el<x_rowind[i+1]; ++el){ // loop over the non-zeros of the matrix
        int j = x_col[el];
        
        // Perform operation
        z[i] += x_trans_data[el] * y[j];
      }
    }
  }
  
  template<class T>
  void Matrix<T>::mul_no_alloc_nn(const Matrix<T>& x, const std::vector<T> &y, std::vector<T> &z){
    // Assert dimensions
    casadi_assert_message(x.size1()==z.size(),"Dimension error. Got x=" << x.dimString() << " and z=" << z.size() << ".");
    casadi_assert_message(x.size2()==y.size(),"Dimension error. Got x=" << x.dimString() << " and y=" << y.size() << ".");
    
    // Direct access to the arrays
    const std::vector<int> &x_colind = x.colind();
    const std::vector<int> &x_row = x.row();
    const std::vector<T> &x_data = x.data();
    
    // loop over the rows of the matrix
    for(int i=0; i<x_colind.size()-1; ++i){
      for(int el=x_colind[i]; el<x_colind[i+1]; ++el){ // loop over the non-zeros of the matrix
        int j = x_row[el];
        z[j] += x_data[el] * y[i];
      }
    }
  }
  
  template<class T>
  void Matrix<T>::mul_no_alloc_nt(const Matrix<T> &x, const Matrix<T>& y_trans, Matrix<T> &z){
    // Assert dimensions
    casadi_assert_message(y_trans.size1()==z.size2(),"Dimension error. Got y_trans=" << y_trans.dimString() << " and z=" << z.dimString() << ".");
    casadi_assert_message(x.size1()==z.size1(),"Dimension error. Got x=" << x.dimString() << " and z=" << z.dimString() << ".");
    casadi_assert_message(y_trans.size2()==x.size2(),"Dimension error. Got y_trans=" << y_trans.dimString() << " and x=" << x.dimString() << ".");
  
    // Direct access to the arrays
    const std::vector<int> &y_rowind = y_trans.colind();
    const std::vector<int> &y_col = y_trans.row();
    const std::vector<T> &y_trans_data = y_trans.data();
    const std::vector<int> &x_colind = x.colind();
    const std::vector<int> &x_row = x.row();
    const std::vector<T> &x_data = x.data();
    const std::vector<int> &z_colind = z.colind();
    const std::vector<int> &z_row = z.row();
    std::vector<T> &z_data = z.data();

    // loop over the rows of the first argument
    for(int i=0; i<y_rowind.size()-1; ++i){
      for(int el=y_rowind[i]; el<y_rowind[i+1]; ++el){ // loop over the non-zeros of the first argument
        int j = y_col[el];
        int el1 = x_colind[i];
        int el2 = z_colind[j];
        while(el1 < x_colind[i+1] && el2 < z_colind[j+1]){ // loop over matching non-zero elements
          int j1 = x_row[el1];
          int i2 = z_row[el2];      
          if(j1==i2){
            z_data[el2++] += y_trans_data[el] * x_data[el1++];
          } else if(j1<i2) {
            el1++;
          } else {
            el2++;
          }
        }
      }
    }
  }

  template<class T>
  void Matrix<T>::mul_no_alloc_tn(const Matrix<T> &x_trans, const Matrix<T> &y, Matrix<T>& z){
    // Assert dimensions
    casadi_assert_message(y.size2()==z.size2(),"Dimension error. Got y=" << y.dimString() << " and z=" << z.dimString() << ".");
    casadi_assert_message(x_trans.size2()==z.size1(),"Dimension error. Got x_trans=" << x_trans.dimString() << " and z=" << z.dimString() << ".");
    casadi_assert_message(y.size1()==x_trans.size1(),"Dimension error. Got y=" << y.dimString() << " and x_trans=" << x_trans.dimString() << ".");
  
    // Direct access to the arrays
    const std::vector<int> &y_colind = y.colind();
    const std::vector<int> &y_row = y.row();
    const std::vector<T> &y_data = y.data();
    const std::vector<int> &x_rowind = x_trans.colind();
    const std::vector<int> &x_col = x_trans.row();
    const std::vector<T> &x_trans_data = x_trans.data();
    const std::vector<int> &z_colind = z.colind();
    const std::vector<int> &z_row = z.row();
    std::vector<T> &z_data = z.data();
  
    // loop over the cols of the resulting matrix
    for(int i=0; i<z_colind.size()-1; ++i){
      for(int el=z_colind[i]; el<z_colind[i+1]; ++el){ // loop over the non-zeros of the resulting matrix
        int j = z_row[el];
        int el1 = y_colind[i];
        int el2 = x_rowind[j];
        while(el1 < y_colind[i+1] && el2 < x_rowind[j+1]){ // loop over non-zero elements
          int j1 = y_row[el1];
          int i2 = x_col[el2];      
          if(j1==i2){
            z_data[el] += y_data[el1++] * x_trans_data[el2++];
          } else if(j1<i2) {
            el1++;
          } else {
            el2++;
          }
        }
      }
    }
  }

  template<class T>
  template<bool Fwd>
  void Matrix<T>::mul_sparsity(Matrix<T> &x_trans, Matrix<T> &y, Matrix<T>& z){
    // Direct access to the arrays
    const std::vector<int> &z_row = z.row();
    const std::vector<int> &z_colind = z.colind();
    const std::vector<int> &y_row = y.row();
    const std::vector<int> &x_col = x_trans.row();
    const std::vector<int> &y_colind = y.colind();
    const std::vector<int> &x_rowind = x_trans.colind();

    // Convert data array to arrays of integers
    bvec_t *y_data = get_bvec_t(y.data());
    bvec_t *x_trans_data = get_bvec_t(x_trans.data());
    bvec_t *z_data = get_bvec_t(z.data());
  
    // loop over the cols of the resulting matrix)
    for(int i=0; i<z_colind.size()-1; ++i){
      for(int el=z_colind[i]; el<z_colind[i+1]; ++el){ // loop over the non-zeros of the resulting matrix
        int j = z_row[el];
        int el1 = y_colind[i];
        int el2 = x_rowind[j];
        while(el1 < y_colind[i+1] && el2 < x_rowind[j+1]){ // loop over non-zero elements
          int j1 = y_row[el1];
          int i2 = x_col[el2];      
          if(j1==i2){
            // | and not & since we are propagating dependencies
            if(Fwd){
              z_data[el] |= y_data[el1] | x_trans_data[el2];
            } else {
              y_data[el1] |= z_data[el];
              x_trans_data[el2] |= z_data[el];
            }
            el1++;
            el2++;
          } else if(j1<i2) {
            el1++;
          } else {
            el2++;
          }
        }
      }
    }
  }

  template<class T>
  T Matrix<T>::quad_form(const Matrix<T>& A, const std::vector<T>& x){
    // Assert dimensions
    casadi_assert_message(x.size()==A.size2() && x.size()==A.size1(),"Dimension mismatch. Got x=" << x.size() << " and A=" << A.dimString());
  
    // Access the internal data of A
    const std::vector<int> &A_colind = A.colind();
    const std::vector<int> &A_row = A.row();
    const std::vector<T> &A_data = A.data();
  
    // Return value
    T ret=0;

    // Loop over the cols of A
    for(int i=0; i<x.size(); ++i){
      // Loop over the nonzeros of A
      for(int el=A_colind[i]; el<A_colind[i+1]; ++el){
        // Get row
        int j = A_row[el];
      
        // Add contribution
        ret += x[i]*A_data[el]*x[j];
      }
    }
  
    return ret;
  }

  template<class T>
  Matrix<T> Matrix<T>::trans() const{
    // quick return if empty or scalar
    if((size1()==0 && size2()==0) || isScalar()) return *this;

    // Create the new sparsity pattern and the mapping
    std::vector<int> mapping;
    Sparsity s = sparsity().transpose(mapping);

    // create the return matrix
    Matrix<T> ret(s);
  
    // Copy the content
    for(int i=0; i<mapping.size(); ++i)
      ret.at(i) = at(mapping[i]);
  
    return ret;
  }

  // template<class T>
  // Matrix<T>::operator const T() const{
  //   return toScalar();
  // }

  template<class T>
  const T Matrix<T>::toScalar() const{
    // Make sure that the matrix is 1-by-1
    casadi_assert_message(isScalar(),"Can only convert 1-by-1 matrices to scalars");

    // return zero or the nonzero element
    if(size()==1)
      return data()[0];
    else
      return casadi_limits<T>::zero;
  }

  template<class T>
  Matrix<T> Matrix<T>::binary(int op, const Matrix<T> &x, const Matrix<T> &y){
    if(x.numel()==1)
      return scalar_matrix(op,x,y);
    else if(y.numel()==1)  
      return matrix_scalar(op,x,y);
    else
      return matrix_matrix(op,x,y);
  }

  template<class T>
  Matrix<T> Matrix<T>::scalar_matrix(int op, const Matrix<T> &x, const Matrix<T> &y){
    // Return value
    Matrix<T> ret(y.sparsity());
  
    // Nonzeros
    std::vector<T>& ret_data = ret.data();
    const std::vector<T>& x_data = x.data();
    const T& x_val = x_data.empty() ? casadi_limits<T>::zero : x.front();
    const std::vector<T>& y_data = y.data();
  
    // Do the operation on all non-zero elements
    for(int el=0; el<y.size(); ++el){
      casadi_math<T>::fun(op,x_val,y_data[el],ret_data[el]);
    }

    // Check the value of the structural zero-entries, if there are any
    if(!y.isDense() && !operation_checker<FX0Checker>(op)){
      // Get the value for the structural zeros
      T fcn_0;
      casadi_math<T>::fun(op,x_val,casadi_limits<T>::zero,fcn_0);
      if(!casadi_limits<T>::isZero(fcn_0)){ // Remove this if?
        ret.makeDense(ret.size1(),ret.size2(),fcn_0);
      }
    }
    
    return ret;
  }

  template<class T>
  Matrix<T> Matrix<T>::matrix_scalar(int op, const Matrix<T> &x, const Matrix<T> &y){
    // Return value
    Matrix<T> ret(x.sparsity());
  
    // Nonzeros
    std::vector<T>& ret_data = ret.data();
    const std::vector<T>& x_data = x.data();
    const std::vector<T>& y_data = y.data();
    const T& y_val = y_data.empty() ? casadi_limits<T>::zero : y.front();
  
    // Do the operation on all non-zero elements
    for(int el=0; el<x.size(); ++el){
      casadi_math<T>::fun(op,x_data[el],y_val,ret_data[el]);
    }

    // Check the value of the structural zero-entries, if there are any
    if(!x.isDense() && !operation_checker<F0XChecker>(op)){
      // Get the value for the structural zeros
      T fcn_0;
      casadi_math<T>::fun(op,casadi_limits<T>::zero,y_val,fcn_0);
      if(!casadi_limits<T>::isZero(fcn_0)){ // Remove this if?
        ret.makeDense(ret.size1(),ret.size2(),fcn_0);
      }
    }
    
    return ret;
  }

  template<class T>
  Matrix<T> Matrix<T>::matrix_matrix(int op, const Matrix<T> &x, const Matrix<T> &y){

    if (!(x.size2() == y.size2() && x.size1() == y.size1())) {
      std::stringstream ss;
      casadi_math<T>::print(op,ss,"lhs","rhs");
      casadi_error("matrix_matrix: dimension mismatch in element-wise matrix operation " << ss.str() <<"." << std::endl << "Left argument has shape " << x.dimString() << ", right has shape " << y.dimString() << ". They should be equal."
                   ); 
    }

    // Get the sparsity pattern of the result (ignoring structural zeros giving rise to nonzero result)
    const Sparsity& x_sp = x.sparsity();
    const Sparsity& y_sp = y.sparsity();
    Sparsity r_sp = x_sp.patternCombine(y_sp, operation_checker<F0XChecker>(op), operation_checker<FX0Checker>(op));

    // Return value
    Matrix<T> r(r_sp);
  
    // Perform the operations elementwise
    if(x_sp==y_sp){
      // Matching sparsities
      casadi_math<T>::fun(op,getPtr(x.data()),getPtr(y.data()),getPtr(r.data()),r_sp.size());
    } else if(y_sp==r_sp){
      // Project first argument
      Matrix<T> x_mod = x(r_sp);
      casadi_math<T>::fun(op,getPtr(x_mod.data()),getPtr(y.data()),getPtr(r.data()),r_sp.size());
    } else if(x_sp==r_sp){
      // Project second argument
      Matrix<T> y_mod = y(r_sp);
      casadi_math<T>::fun(op,getPtr(x.data()),getPtr(y_mod.data()),getPtr(r.data()),r_sp.size());
    } else {
      // Project both arguments
      Matrix<T> x_mod = x(r_sp);
      Matrix<T> y_mod = y(r_sp);
      casadi_math<T>::fun(op,getPtr(x_mod.data()),getPtr(y_mod.data()),getPtr(r.data()),r_sp.size());
    }

    // Handle structural zeros giving rise to nonzero result, e.g. cos(0) == 1
    if(!r.isDense() && !operation_checker<F00Checker>(op)){
      // Get the value for the structural zeros
      T fcn_0;
      casadi_math<T>::fun(op,casadi_limits<T>::zero,casadi_limits<T>::zero,fcn_0);
      r.makeDense(r.size1(),r.size2(),fcn_0);
    }
  
    return r;
  }

  template<class T>
  Matrix<T> Matrix<T>::triplet(const std::vector<int>& row, const std::vector<int>& col, const std::vector<T>& d) {
    return triplet(row,col,d,*std::max_element(row.begin(),row.end()),*std::max_element(col.begin(),col.end()));
  }

  template<class T>
  Matrix<T> Matrix<T>::triplet(const std::vector<int>& row, const std::vector<int>& col, const std::vector<T>& d, const std::pair<int,int>& rc) {
    return triplet(row,col,d,rc.first,rc.second);
  }

  template<class T>
  Matrix<T> Matrix<T>::triplet(const std::vector<int>& row, const std::vector<int>& col, const std::vector<T>& d, int nrow, int ncol) {
    casadi_assert_message(col.size()==row.size() && col.size()==d.size(),"Argument error in Matrix<T>::sparse(row,col,d): supplied lists must all be of equal length, but got: " << row.size() << ", " << col.size()  << " and " << d.size());
    std::vector<int> mapping;
    Sparsity sp = sp_triplet(nrow,ncol,row,col,mapping);
    std::vector<T> v(mapping.size());
    for(int k=0; k<v.size(); ++k) v[k] = d[mapping[k]];
    return Matrix<T>(sp,v);
  }

  template<class T>
  Matrix<T> Matrix<T>::repmat(const T& x, const Sparsity& sp){
    return Matrix<T>(sp,x);
  }

  template<class T>
  Matrix<T> Matrix<T>::repmat(const Matrix<T>& x, const Sparsity& sp){
    casadi_assert_message(x.isScalar(),"repmat(Matrix<T> x,Sparsity sp) only defined for scalar x");
    return Matrix<T>(sp,x.toScalar());
  }

  template<class T>
  Matrix<T> Matrix<T>::repmat(const Matrix<T>& x, const std::pair<int,int>& rc){
    return repmat(x,rc.first,rc.second);
  }

  template<class T>
  Matrix<T> Matrix<T>::repmat(const Matrix<T>& x, int nrow, int ncol){
    if(x.isScalar()){
      if(x.isDense()){
        return Matrix<T>(nrow,ncol,x.toScalar());
      } else {
        return sparse(nrow,ncol);
      }
    } else {
      return vertcat(std::vector< Matrix<T> >(nrow,horzcat(std::vector< Matrix<T> >(ncol,x))));
    }
  }

  template<class T>
  Matrix<T> Matrix<T>::eye(int n){
    return Matrix<T>(Sparsity::createDiagonal(n),1);
  }

  template<class T>
  Matrix<T> Matrix<T>::inf(const Sparsity& sp){
    casadi_assert_message(std::numeric_limits<T>::has_infinity,"Datatype cannot represent infinity");
    return Matrix<T>(sp,std::numeric_limits<T>::infinity());
  }


  template<class T>
  Matrix<T> Matrix<T>::inf(const std::pair<int,int>& rc){
    return inf(rc.first, rc.second);
  }

  template<class T>
  Matrix<T> Matrix<T>::inf(int nrow, int ncol){
    return inf(sp_dense(nrow,ncol));
  }

  template<class T>
  Matrix<T> Matrix<T>::nan(const Sparsity& sp){
    casadi_assert_message(std::numeric_limits<T>::has_quiet_NaN,"Datatype cannot represent not-a-number");
    return Matrix<T>(sp,std::numeric_limits<T>::quiet_NaN());
  }

  template<class T>
  Matrix<T> Matrix<T>::nan(const std::pair<int,int>& rc){
    return nan(rc.first, rc.second);
  }

  template<class T>
  Matrix<T> Matrix<T>::nan(int nrow, int ncol){
    return nan(sp_dense(nrow,ncol));
  }

  template<class T>
  void Matrix<T>::append(const Matrix<T>& y){
    // Quick return if expr is empty
    if(size2()==0 && size1()==0){
      *this=y;
      return;
    }

    // Quick return if empty
    if(y.size2()==0 && y.size1()==0) return;

    // Appending can be done efficeintly if vectors
    if(isVector()){
      // Append the sparsity pattern vertically
      sparsityRef().append(y.sparsity());
  
      // Add the non-zeros at the end
      data().insert(end(),y.begin(),y.end());
    } else {
      // Fall back on vertical concatenation
      *this = vertcat(*this,y);
    }
  }

  template<class T>
  void Matrix<T>::appendColumns(const Matrix<T>& y){

    // Quick return if expr is empty
    if(size2()==0 && size1()==0){
      *this=y;
      return;
    }

    // Quick return if empty
    if(y.size2()==0 && y.size1()==0) return;
  
    // Append the sparsity pattern
    sparsityRef().appendColumns(y.sparsity());
  
    // Add the non-zeros at the end
    data().insert(end(),y.begin(),y.end());
  }

  template<class T>
  NonZeroIterator<T>::NonZeroIterator(const Matrix<T> & m) 
    : m_(m) {
    nz.i  = 0;
    nz.j  = 0;
    nz.k  = 0;
  }

  template<class T>
  bool NonZeroIterator<T>::operator==(const NonZeroIterator<T>& rhs) {return (m_ == rhs.m_) && (nz.k==rhs.nz.k);}

  template<class T>
  NonZero<T>& NonZeroIterator<T>::operator*() {
  
    return nz;
  }
    
  template<class T>
  NonZeroIterator<T>& NonZeroIterator<T>::operator++() {
    nz.k ++;
  
    if (nz.k < m_.size()) {
      nz.j = m_.row()[nz.k];
      nz.el = m_.data()[nz.k];
      while (nz.k>=m_.colind(nz.i)) {nz.i++; }
    }
    return *this;
  }

  template<class T>
  NonZeroIterator<T> NonZeroIterator<T>::begin() {
    NonZeroIterator<T> it = NonZeroIterator<T>(m_);
    return it;

  }
  template<class T>
  NonZeroIterator<T> NonZeroIterator<T>::end() {
    NonZeroIterator<T> it = NonZeroIterator<T>(m_);
    it.nz.k = m_.size()-1;
    return it;
  }


} // namespace CasADi

#endif // MATRIX_IMPL_HPP

