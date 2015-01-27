/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010-2014 Joel Andersson, Joris Gillis, Moritz Diehl,
 *                            K.U. Leuven. All rights reserved.
 *    Copyright (C) 2011-2014 Greg Horn
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


#ifndef CASADI_MATRIX_IMPL_HPP
#define CASADI_MATRIX_IMPL_HPP

// The declaration of the class is in a separate file
#include "matrix.hpp"
#include "matrix_tools.hpp"
#include "../std_vector_tools.hpp"

/// \cond INTERNAL

namespace casadi {
  // Implementations

  template<typename DataType>
  const DataType& Matrix<DataType>::elem(int rr, int cc) const {
    int ind = sparsity().getNZ(rr, cc);
    if (ind==-1)
      return casadi_limits<DataType>::zero;
    else
      return at(ind);
  }

  template<typename DataType>
  int Matrix<DataType>::stream_precision_ = 6;
  template<typename DataType>
  int Matrix<DataType>::stream_width_ = 0;
  template<typename DataType>
  bool Matrix<DataType>::stream_scientific_ = false;

  template<typename DataType>
  DataType& Matrix<DataType>::elem(int rr, int cc) {
    int oldsize = sparsity().nnz();
    int ind = sparsityRef().addNZ(rr, cc);
    if (oldsize != sparsity().nnz())
      data().insert(begin()+ind, DataType(0));
    return at(ind);
  }

  template<typename DataType>
  bool Matrix<DataType>::__nonzero__() const {
    if (numel()!=1) {casadi_error("Only scalar Matrix could have a truth value, but you "
                                  "provided a shape" << dimString());}
    return at(0)!=0;
  }

  template<typename DataType>
  bool Matrix<DataType>::isSlice(bool ind1) const {
    throw CasadiException("\"isSlice\" not defined for instantiation");
    return false;
  }

  template<typename DataType>
  Slice Matrix<DataType>::toSlice(bool ind1) const {
    throw CasadiException("\"toSlice\" not defined for instantiation");
    return Slice();
  }

  template<typename DataType>
  const Matrix<DataType> Matrix<DataType>::getSub(bool ind1,
                                                  const Slice& rr, const Slice& cc) const {
    // Both are scalar
    if (rr.isScalar(size1()) && cc.isScalar(size2())) {
      int k = sparsity().getNZ(rr.toScalar(size1()), cc.toScalar(size2()));
      if (k>=0) {
        return at(k);
      } else {
        return Matrix<DataType>(1, 1);
      }
    }

    // Fall back on IMatrix-IMatrix
    return getSub(ind1, rr.getAll(size1(), ind1), cc.getAll(size2(), ind1));
  }

  template<typename DataType>
  const Matrix<DataType>
  Matrix<DataType>::getSub(bool ind1, const Slice& rr, const Matrix<int>& cc) const {
    // Fall back on IMatrix-IMatrix
    return getSub(ind1, rr.getAll(size1(), ind1), cc);
  }

  template<typename DataType>
  const Matrix<DataType>
  Matrix<DataType>::getSub(bool ind1, const Matrix<int>& rr, const Slice& cc) const {
    // Fall back on IMatrix-IMatrix
    return getSub(ind1, rr, cc.getAll(size2(), ind1));
  }

  template<typename DataType>
  const Matrix<DataType>
  Matrix<DataType>::getSub(bool ind1, const Matrix<int>& rr, const Matrix<int>& cc) const {
    // Scalar
    if (rr.isScalar(true) && cc.isScalar(true)) {
      return getSub(ind1, rr.toSlice(ind1), cc.toSlice(ind1));
    }

    // Row vector rr (e.g. in MATLAB) is transposed to column vector
    if (rr.size1()==1 && rr.size2()>1) {
      return getSub(ind1, rr.T(), cc);
    }

    // Row vector cc (e.g. in MATLAB) is transposed to column vector
    if (cc.size1()==1 && cc.size2()>1) {
      return getSub(ind1, rr, cc.T());
    }

    casadi_assert_message(rr.isDense() && rr.isVector(),
                          "Marix::getSub: First index must be a dense vector");
    casadi_assert_message(cc.isDense() && cc.isVector(),
                          "Marix::getSub: Second index must be a dense vector");

    // Get the sparsity pattern - does bounds checking
    std::vector<int> mapping;
    Sparsity sp = sparsity().sub(rr.data(), cc.data(), mapping, ind1);

    // Copy nonzeros and return
    Matrix<DataType> ret(sp);
    for (int k=0; k<mapping.size(); ++k) ret.at(k) = at(mapping[k]);

    // Return (RVO)
    return ret;
  }

  template<typename DataType>
  const Matrix<DataType> Matrix<DataType>::getSub(bool ind1, const Slice& rr) const {
    // Scalar
    if (rr.isScalar(numel())) {
      int r = rr.toScalar(numel());
      int k = sparsity().getNZ(r % size1(), r / size1());
      if (k>=0) {
        return at(k);
      } else {
        return Matrix<DataType>(1, 1);
      }
    }

    // Fall back on IMatrix
    return getSub(ind1, rr.getAll(numel(), ind1));
  }

  template<typename DataType>
  const Matrix<DataType> Matrix<DataType>::getSub(bool ind1, const Matrix<int>& rr) const {
    // Scalar
    if (rr.isScalar(true)) {
      return getSub(ind1, rr.toSlice(ind1));
    }

    // If the indexed matrix is dense, use nonzero indexing
    if (isDense()) {
      return getNZ(ind1, rr);
    }

    // Get the sparsity pattern - does bounds checking
    std::vector<int> mapping;
    Sparsity sp = sparsity().sub(rr.data(), rr.sparsity(), mapping, ind1);

    // Copy nonzeros and return
    Matrix<DataType> ret(sp);
    for (int k=0; k<mapping.size(); ++k) ret.at(k) = at(mapping[k]);

    // Return (RVO)
    return ret;
  }

  template<typename DataType>
  const Matrix<DataType> Matrix<DataType>::getSub(bool ind1, const Sparsity& sp) const {
    casadi_assert_message(shape()==sp.shape(),
                          "getSub(Sparsity sp): shape mismatch. This matrix has shape "
                          << shape() << ", but supplied sparsity index has shape "
                          << sp.shape() << ".");
    return setSparse(sp);
  }

  template<typename DataType>
  void Matrix<DataType>::setSub(const Matrix<DataType>& m, bool ind1,
                                const Slice& rr, const Slice& cc) {
    // Both are scalar
    if (rr.isScalar(size1()) && cc.isScalar(size2()) && m.isDense()) {
      elem(rr.toScalar(size1()), cc.toScalar(size2())) = m.toScalar();
      return;
    }

    // Fall back on (IMatrix, IMatrix)
    setSub(m, ind1, rr.getAll(size1(), ind1), cc.getAll(size2(), ind1));
  }

  template<typename DataType>
  void Matrix<DataType>::setSub(const Matrix<DataType>& m, bool ind1,
                                const Slice& rr, const Matrix<int>& cc) {
    // Fall back on (IMatrix, IMatrix)
    setSub(m, ind1, rr.getAll(size1(), ind1), cc);
  }

  template<typename DataType>
  void Matrix<DataType>::setSub(const Matrix<DataType>& m, bool ind1,
                                const Matrix<int>& rr, const Slice& cc) {
    // Fall back on (IMatrix, IMatrix)
    setSub(m, ind1, rr, cc.getAll(size2(), ind1));
  }

  template<typename DataType>
  void Matrix<DataType>::setSub(const Matrix<DataType>& m, bool ind1,
                                const Matrix<int>& rr, const Matrix<int>& cc) {
    // Scalar
    if (rr.isScalar(true) && cc.isScalar(true) && m.isDense()) {
      return setSub(m, ind1, rr.toSlice(ind1), cc.toSlice(ind1));
    }

    // Row vector rr (e.g. in MATLAB) is transposed to column vector
    if (rr.size1()==1 && rr.size2()>1) {
      return setSub(m, ind1, rr.T(), cc);
    }

    // Row vector cc (e.g. in MATLAB) is transposed to column vector
    if (cc.size1()==1 && cc.size2()>1) {
      return setSub(m, ind1, rr, cc.T());
    }

    // Make sure rr and cc are dense vectors
    casadi_assert_message(rr.isDense() && rr.isVector(),
                          "Matrix::setSub: First index not dense vector");
    casadi_assert_message(cc.isDense() && cc.isVector(),
                          "Matrix::setSub: Second index not dense vector");

    // Assert dimensions of assigning matrix
    if (rr.size1() != m.size1() || cc.size1() != m.size2()) {
      if (m.isScalar()) {
        // m scalar means "set all"
        return setSub(repmat(m, rr.size1(), cc.size1()), ind1, rr, cc);
      } else if (rr.size1() == m.size2() && cc.size1() == m.size1()
                 && std::min(m.size1(), m.size2()) == 1) {
        // m is transposed if necessary
        return setSub(m.T(), ind1, rr, cc);
      } else {
        // Error otherwise
        casadi_error("Dimension mismatch." << "lhs is " << rr.size1() << "-by-"
                     << cc.size1() << ", while rhs is " << m.shape());
      }
    }

    // Dimensions
    int sz1 = size1(), sz2 = size2();

    // Report out-of-bounds
    if (!inBounds(rr.data(), -sz1+ind1, sz1+ind1)) {
      casadi_error("setSub[., r, c] out of bounds. Your rr contains "
                   << *std::min_element(rr.begin(), rr.end()) << " up to "
                   << *std::max_element(rr.begin(), rr.end())
                   << ", which is outside the range [" << -sz1+ind1 << ","<< sz1+ind1 <<  ").");
    }
    if (!inBounds(cc.data(), -sz2+ind1, sz2+ind1)) {
      casadi_error("setSub [., r, c] out of bounds. Your cc contains "
                   << *std::min_element(cc.begin(), cc.end()) << " up to "
                   << *std::max_element(cc.begin(), cc.end())
                   << ", which is outside the range [" << -sz2+ind1 << ","<< sz2+ind1 <<  ").");
    }

    // If we are assigning with something sparse, first remove existing entries
    if (!m.isDense()) {
      erase(rr.data(), cc.data(), ind1);
    }

    // Collect all assignments
    IMatrix el(m.sparsity());
    for (int j=0; j<el.size2(); ++j) { // Loop over columns of m
      int this_j = cc.at(j) - ind1; // Corresponding column in this
      if (this_j<0) this_j += sz2;
      for (int k=el.colind(j); k<el.colind(j+1); ++k) { // Loop over rows of m
        int i = m.row(k);
        int this_i = rr.at(i) - ind1; // Corresponding row in this
        if (this_i<0) this_i += sz1;
        el.at(k) = this_i + this_j*sz1;
      }
    }
    return setSub(m, false, el);
  }

  template<typename DataType>
  void Matrix<DataType>::setSub(const Matrix<DataType>& m, bool ind1, const Slice& rr) {
    // Scalar
    if (rr.isScalar(numel()) && m.isDense()) {
      int r = rr.toScalar(numel());
      elem(r % size1(), r / size1()) = m.toScalar();
      return;
    }

    // Fall back on IMatrix
    setSub(m, ind1, rr.getAll(numel(), ind1));
  }

  template<typename DataType>
  void Matrix<DataType>::setSub(const Matrix<DataType>& m, bool ind1, const Matrix<int>& rr) {
    // Scalar
    if (rr.isScalar(true) && m.isDense()) {
      return setSub(m, ind1, rr.toSlice(ind1));
    }

    // Assert dimensions of assigning matrix
    if (rr.sparsity() != m.sparsity()) {
      if (rr.shape() == m.shape()) {
        // Remove submatrix to be replaced
        erase(rr.data(), ind1);

        // Find the intersection between rr's and m's sparsity patterns
        Sparsity sp = rr.sparsity() * m.sparsity();

        // Project both matrices to this sparsity
        return setSub(m.setSparse(sp), ind1, rr.setSparse(sp));
      } else if (m.isScalar()) {
        // m scalar means "set all"
        if (m.isDense()) {
          return setSub(repmat(m, rr.sparsity()), ind1, rr);
        } else {
          return setSub(Matrix<DataType>(rr.shape()), ind1, rr);
        }
      } else if (rr.size1() == m.size2() && rr.size2() == m.size1()
                 && std::min(m.size1(), m.size2()) == 1) {
        // m is transposed if necessary
        return setSub(m.T(), ind1, rr);
      } else {
        // Error otherwise
        casadi_error("Dimension mismatch." << "lhs is " << rr.shape()
                     << ", while rhs is " << m.shape());
      }
    }

    // Dimensions of this
    int sz1 = size1(), sz2 = size2(), sz = nnz(), nel = numel(), rrsz = rr.nnz();

    // Quick return if nothing to set
    if (rrsz==0) return;

    // Check bounds
    if (!inBounds(rr.data(), -nel+ind1, nel+ind1)) {
      casadi_error("setSub[rr] out of bounds. Your rr contains "
                   << *std::min_element(rr.begin(), rr.end()) << " up to "
                   << *std::max_element(rr.begin(), rr.end())
                   << ", which is outside the range [" << -nel+ind1 << ","<< nel+ind1 <<  ").");
    }

    // Dense mode
    if (isDense() && m.isDense()) {
      return setNZ(m, ind1, rr);
    }

    // Construct new sparsity pattern
    std::vector<int> new_row, new_col, nz(rr.data());
    new_row.reserve(sz+rrsz);
    new_col.reserve(sz+rrsz);
    nz.reserve(rrsz);
    sparsity().getTriplet(new_row, new_col);
    for (std::vector<int>::iterator i=nz.begin(); i!=nz.end(); ++i) {
      if (ind1) (*i)--;
      if (*i<0) *i += nel;
      new_row.push_back(*i % sz1);
      new_col.push_back(*i / sz1);
    }
    Sparsity sp = Sparsity::triplet(sz1, sz2, new_row, new_col);

    // If needed, update pattern
    if (sp != sparsity()) *this = setSparse(sp);

    // Find the nonzeros corresponding to rr
    sparsity().getNZ(nz);

    // Carry out the assignments
    for (int i=0; i<nz.size(); ++i) {
      at(nz[i]) = m.at(i);
    }
  }

  template<typename DataType>
  void Matrix<DataType>::setSub(const Matrix<DataType>& m, bool ind1, const Sparsity& sp) {
    casadi_assert_message(shape()==sp.shape(),
                          "setSub(Sparsity sp): shape mismatch. This matrix has shape "
                          << shape() << ", but supplied sparsity index has shape "
                          << sp.shape() << ".");
    std::vector<int> ii = sp.find();
    if (m.isScalar()) {
      (*this)(ii) = densify(m);
    } else {
      (*this)(ii) = densify(m(ii));
    }
  }

  template<typename DataType>
  const Matrix<DataType> Matrix<DataType>::getNZ(bool ind1, const Slice& kk) const {
    // Scalar
    if (kk.isScalar(nnz())) {
      return at(kk.toScalar(nnz()));
    }

    // Fall back on IMatrix
    return getNZ(ind1, kk.getAll(nnz(), ind1));
  }

  template<typename DataType>
  const Matrix<DataType> Matrix<DataType>::getNZ(bool ind1, const Matrix<int>& kk) const {
    // Scalar
    if (kk.isScalar(true)) {
      return getNZ(ind1, kk.toSlice(ind1));
    }

    // Get nonzeros of kk
    const std::vector<int>& k = kk.data();
    int sz = nnz();

    // Check bounds
    if (!inBounds(k, -sz+ind1, sz+ind1)) {
      casadi_error("getNZ[kk] out of bounds. Your kk contains "
                   << *std::min_element(k.begin(), k.end()) << " up to "
                   << *std::max_element(k.begin(), k.end())
                   << ", which is outside the range [" << -sz+ind1 << ","<< sz+ind1 <<  ").");
    }

    // Copy nonzeros
    Matrix<DataType> ret = zeros(kk.sparsity());
    for (int el=0; el<k.size(); ++el) {
      int k_el = k[el]-ind1;
      ret.at(el) = at(k_el>=0 ? k_el : k_el+sz);
    }
    return ret;
  }

  template<typename DataType>
  void Matrix<DataType>::setNZ(const Matrix<DataType>& m, bool ind1, const Slice& kk) {
    // Scalar
    if (kk.isScalar(nnz())) {
      at(kk.toScalar(nnz())) = m.toScalar();
      return;
    }

    // Fallback on IMatrix
    setNZ(m, ind1, kk.getAll(nnz(), ind1));
  }

  template<typename DataType>
  void Matrix<DataType>::setNZ(const Matrix<DataType>& m, bool ind1, const Matrix<int>& kk) {
    // Scalar
    if (kk.isScalar(true)) {
      return setNZ(m, ind1, kk.toSlice(ind1));
    }

    // Assert dimensions of assigning matrix
    if (kk.sparsity() != m.sparsity()) {
      if (m.isScalar()) {
        // m scalar means "set all"
        if (!m.isDense()) return; // Nothing to set
        return setNZ(repmat(m, kk.sparsity()), ind1, kk);
      } else if (kk.shape() == m.shape()) {
        // Project sparsity if needed
        return setNZ(m.setSparse(kk.sparsity()), ind1, kk);
      } else if (kk.size1() == m.size2() && kk.size2() == m.size1()
                 && std::min(m.size1(), m.size2()) == 1) {
        // m is transposed if necessary
        return setNZ(m.T(), ind1, kk);
      } else {
        // Error otherwise
        casadi_error("Dimension mismatch." << "lhs is " << kk.shape()
                     << ", while rhs is " << m.shape());
      }
    }

    // Get nonzeros
    const std::vector<int>& k = kk.data();
    int sz = nnz();

    // Check bounds
    if (!inBounds(k, -sz+ind1, sz+ind1)) {
      casadi_error("setNZ[kk] out of bounds. Your kk contains "
                   << *std::min_element(k.begin(), k.end()) << " up to "
                   << *std::max_element(k.begin(), k.end())
                   << ", which is outside the range [" << -sz+ind1 << ","<< sz+ind1 <<  ").");
    }

    // Set nonzeros, ignoring negative indices
    for (int el=0; el<k.size(); ++el) {
      int k_el = k[el]-ind1;
      at(k_el>=0 ? k_el : k_el+sz) = m.at(el);
    }
  }

  template<typename DataType>
  void Matrix<DataType>::makeDense(const DataType& val) {
    // Quick return if possible
    if (isDense()) return;

    // Get sparsity pattern
    int nrow = size1();
    int ncol = size2();
    const int* colind = this->colindPtr();
    const int* row = this->rowPtr();

    // Resize data and copy
    data_.resize(nrow*ncol, val);

    // Loop over the columns in reverse order
    for (int cc=ncol-1; cc>=0; --cc) {
      // Loop over nonzero elements of the column in reverse order
      for (int el=colind[cc+1]-1; el>=colind[cc]; --el) {
        int rr = row[el];
        int new_el = cc*nrow + rr;
        if (el==new_el) break; // Already done, the rest of the elements must be in the same place
        std::swap(data_[new_el], data_[el]);
      }
    }

    // Update the sparsity pattern
    sparsity_ = Sparsity::dense(shape());
  }

  template<typename DataType>
  void Matrix<DataType>::makeSparse(double tol) {
    // Quick return if there are no entries to be removed
    bool remove_nothing = true;
    for (typename std::vector<DataType>::iterator it=begin(); it!=end() && remove_nothing; ++it) {
      remove_nothing = !casadi_limits<DataType>::isAlmostZero(*it, tol);
    }
    if (remove_nothing) return;

    // Get the current sparsity pattern
    int size1 = this->size1();
    int size2 = this->size2();
    const int* colind = this->colindPtr();
    const int* row = this->rowPtr();

    // Construct the new sparsity pattern
    std::vector<int> new_colind(1, 0), new_row;

    // Loop over the columns
    for (int cc=0; cc<size2; ++cc) {
      // Loop over existing nonzeros
      for (int el=colind[cc]; el<colind[cc+1]; ++el) {
        // If it is not known to be a zero
        if (!casadi_limits<DataType>::isAlmostZero(data_[el], tol)) {
          // Save the nonzero in its new location
          data_[new_row.size()] = data_[el];

          // Add to pattern
          new_row.push_back(row[el]);
        }
      }
      // Save the new column offset
      new_colind.push_back(new_row.size());
    }

    // Trim the data vector
    data_.resize(new_row.size());

    // Update the sparsity pattern
    sparsity_ = Sparsity(size1, size2, new_colind, new_row);
  }

  template<typename DataType>
  Matrix<DataType>::Matrix() : sparsity_(Sparsity(0, 0)) {
  }

  template<typename DataType>
  Matrix<DataType>::Matrix(const Matrix<DataType>& m) : sparsity_(m.sparsity_), data_(m.data_) {
  }

  template<typename DataType>
  Matrix<DataType>::Matrix(const std::vector<DataType>& x) :
      sparsity_(Sparsity::dense(x.size(), 1)), data_(x) {
  }

  template<typename DataType>
  Matrix<DataType>::Matrix(const std::vector<DataType>& x, int nrow, int ncol) :
      sparsity_(Sparsity::dense(nrow, ncol)), data_(x) {
    casadi_assert_message(x.size() == nrow*ncol, "Dimension mismatch." << std::endl
                          << "You supplied a vector of length " << x.size() << ", but " << nrow
                          << " x " << ncol << " = " << nrow*ncol);
  }

  template<typename DataType>
  Matrix<DataType>& Matrix<DataType>::operator=(const Matrix<DataType>& m) {
    sparsity_ = m.sparsity_;
    data_ = m.data_;
    return *this;
  }

  template<typename DataType>
  std::string Matrix<DataType>::className() { return matrixName<DataType>(); }

  template<typename DataType>
  void Matrix<DataType>::printScalar(std::ostream &stream, bool trailing_newline) const {
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

    if (nnz()==0) {
      stream << "00";
    } else {
      stream << toScalar();
    }

    if (trailing_newline) stream << std::endl;
    stream << std::flush;
    stream.precision(precision);
    stream.width(width);
    stream.flags(flags);
  }

  template<typename DataType>
  void Matrix<DataType>::printVector(std::ostream &stream, bool trailing_newline) const {
    casadi_assert_message(isVector(), "Not a vector");

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
    const int* r = rowPtr();
    int sz = nnz();

    // Nonzero
    int el=0;

    // Loop over rows
    stream << "[";
    for (int rr=0; rr<size1(); ++rr) {
      // Add delimiter
      if (rr!=0) stream << ", ";

      // Check if nonzero
      if (el<sz && rr==r[el]) {
        stream << at(el++);
      } else {
        stream << "00";
      }
    }
    stream << "]";

    if (trailing_newline) stream << std::endl;
    stream << std::flush;
    stream.precision(precision);
    stream.width(width);
    stream.flags(flags);
  }

  template<typename DataType>
  void Matrix<DataType>::printDense(std::ostream &stream, bool trailing_newline) const {
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
    const int* cptr = this->colindPtr();
    int ncol = size2();
    std::vector<int> cind(cptr, cptr+ncol+1);

    // Loop over rows
    for (int rr=0; rr<size1(); ++rr) {
      // Beginning of row
      if (rr==0) {
        if (!oneliner) stream << std::endl;
        stream << "[[";
      } else {
        stream << " [";
      }

      // Loop over columns
      for (int cc=0; cc<ncol; ++cc) {
        // Separating comma
        if (cc>0) stream << ", ";

        // Check if nonzero
        if (cind[cc]<colind(cc+1) && row(cind[cc])==rr) {
          stream << data().at(cind[cc]++);
        } else {
          stream << "00";
        }
      }

      // End of row
      if (rr<size1()-1) {
        stream << "], ";
        if (!oneliner) stream << std::endl;
      } else {
        stream << "]]";
      }
    }

    if (trailing_newline) stream << std::endl;
    stream << std::flush;
    stream.precision(precision);
    stream.width(width);
    stream.flags(flags);
  }

  template<typename DataType>
  void Matrix<DataType>::printSparse(std::ostream &stream, bool trailing_newline) const {
    if (nnz()==0) {
      stream << "all zero sparse: " << size1() << "-by-" << size2();
    } else {
      stream << "sparse: " << size1() << "-by-" << size2() << ", " << nnz() << " nnz";
      for (int cc=0; cc<size2(); ++cc) {
        for (int el=colind(cc); el<colind(cc+1); ++el) {
          int rr=row(el);
          stream << std::endl << " (" << rr << ", " << cc << ") -> " << at(el);
        }
      }
    }
    if (trailing_newline) stream << std::endl;
    stream << std::flush;
  }

  template<typename DataType>
  void Matrix<DataType>::print(std::ostream &stream, bool trailing_newline) const {
    if (isEmpty()) {
      stream << "[]";
    } else if (numel()==1) {
      printScalar(stream, false);
    } else if (isVector()) {
      printVector(stream, false);
    } else if (std::max(size1(), size2())<=10 || static_cast<double>(nnz())/numel()>=0.5) {
      // if "small" or "dense"
      printDense(stream, false);
    } else {
      printSparse(stream, false);
    }
    if (trailing_newline) stream << std::endl;
  }

  template<typename DataType>
  void Matrix<DataType>::repr(std::ostream &stream, bool trailing_newline) const {
    stream << className() << "(";
    print(stream, false);
    stream << ")";
    if (trailing_newline) stream << std::endl;
    stream << std::flush;
  }

  template<typename DataType>
  void Matrix<DataType>::reserve(int nnz) {
    reserve(nnz, size2());
  }

  template<typename DataType>
  void Matrix<DataType>::reserve(int nnz, int ncol) {
    data().reserve(nnz);
    sparsity_.reserve(nnz, ncol);
  }

  template<typename DataType>
  void Matrix<DataType>::resize(int nrow, int ncol) {
    sparsity_.resize(nrow, ncol);
  }

  template<typename DataType>
  void Matrix<DataType>::clear() {
    sparsity_ = Sparsity(0, 0);
    data().clear();
  }

  template<typename DataType>
  Matrix<DataType>::Matrix(double val) :
      sparsity_(Sparsity::dense(1, 1)), data_(std::vector<DataType>(1, val)) {
  }

  template<typename DataType>
  Matrix<DataType>::Matrix(const std::vector< std::vector<DataType> >& d) {
    // Get dimensions
    int nrow=d.size();
    int ncol=d.empty() ? 1 : d.front().size();

    // Assert consistency
    for (int rr=0; rr<nrow; ++rr) {
      casadi_assert_message(ncol==d[rr].size(),
        "Matrix<DataType>::Matrix(const std::vector< std::vector<DataType> >& d): "
        "shape mismatch" << std::endl
        << "Attempting to construct a matrix from a nested list." << std::endl
        << "I got convinced that the desired size is ("<< nrow << " x " << ncol
        << " ), but now I encounter a vector of size ("
        << d[rr].size() <<  " )" << std::endl);
    }

    // Form matrix
    sparsity_ = Sparsity::dense(nrow, ncol);
    data().resize(nrow*ncol);
    typename std::vector<DataType>::iterator it=begin();
    for (int cc=0; cc<ncol; ++cc) {
      for (int rr=0; rr<nrow; ++rr) {
        *it++ = d[rr][cc];
      }
    }
  }

  template<typename DataType>
  Matrix<DataType>::Matrix(const Sparsity& sp) :
      sparsity_(sp), data_(sp.nnz(), 0) {
  }


  template<typename DataType>
  Matrix<DataType>::Matrix(int nrow, int ncol) : sparsity_(nrow, ncol) {
  }

  template<typename DataType>
  Matrix<DataType>::Matrix(const std::pair<int, int>& rc) : sparsity_(rc) {
  }

  template<typename DataType>
  Matrix<DataType>::Matrix(const Sparsity& sp, const DataType& val, bool dummy) :
      sparsity_(sp), data_(sp.nnz(), val) {
  }

  template<typename DataType>
  Matrix<DataType>::Matrix(const Sparsity& sp, const std::vector<DataType>& d, bool dummy) :
      sparsity_(sp), data_(d) {
    casadi_assert_message(sp.nnz()==d.size(), "Size mismatch." << std::endl
                          << "You supplied a sparsity of " << sp.dimString()
                          << ", but the supplied vector is of length " << d.size());
  }

  template<typename DataType>
  Matrix<DataType>::Matrix(const Sparsity& sp, const Matrix<DataType>& d) {
    if (d.isScalar()) {
      *this = Matrix<DataType>(sp, d.toScalar(), false);
    } else if (d.isVector() || d.size1()==1) {
      casadi_assert(sp.nnz()==d.numel());
      if (d.isDense()) {
        *this = Matrix<DataType>(sp, d.data(), false);
      } else {
        *this = Matrix<DataType>(sp, densify(d).data(), false);
      }
    } else {
      casadi_error("Matrix(Sparsisty, Matrix): Only allowed for scalars and vectors");
    }
  }

  template<typename DataType>
  void Matrix<DataType>::setZero() {
    setAll(0);
  }

  template<typename DataType>
  void Matrix<DataType>::setAll(const DataType& val) {
    std::fill(begin(), end(), val);
  }

  template<typename DataType>
  Matrix<DataType> Matrix<DataType>::unary(int op, const Matrix<DataType> &x) {
    // Return value
    Matrix<DataType> ret(x.sparsity());

    // Nonzeros
    std::vector<DataType>& ret_data = ret.data();
    const std::vector<DataType>& x_data = x.data();

    // Do the operation on all non-zero elements
    for (int el=0; el<x.nnz(); ++el) {
      casadi_math<DataType>::fun(op, x_data[el], x_data[el], ret_data[el]);
    }

    // Check the value of the structural zero-entries, if there are any
    if (!x.isDense() && !operation_checker<F0XChecker>(op)) {
      // Get the value for the structural zeros
      DataType fcn_0;
      casadi_math<DataType>::fun(op, 0, 0, fcn_0);
      if (!casadi_limits<DataType>::isZero(fcn_0)) { // Remove this if?
        ret.makeDense(fcn_0);
      }
    }

    return ret;
  }

  template<typename DataType>
  Matrix<DataType> Matrix<DataType>::operator-() const {
    return unary(OP_NEG, *this);
  }

  template<typename DataType>
  Matrix<DataType> Matrix<DataType>::operator+() const {
    return *this;
  }

  template<typename DataType>
  Matrix<DataType> Matrix<DataType>::zz_plus(const Matrix<DataType> &y) const {
    return binary(OP_ADD, *this, y);
  }

  template<typename DataType>
  Matrix<DataType> Matrix<DataType>::zz_minus(const Matrix<DataType> &y) const {
    return binary(OP_SUB, *this, y);
  }

  template<typename DataType>
  Matrix<DataType> Matrix<DataType>::zz_times(const Matrix<DataType> &y) const {
    return binary(OP_MUL, *this, y);
  }

  template<typename DataType>
  Matrix<DataType> Matrix<DataType>::zz_rdivide(const Matrix<DataType> &y) const {
    return binary(OP_DIV, *this, y);
  }

  template<typename DataType>
  Matrix<DataType> Matrix<DataType>::zz_lt(const Matrix<DataType> &y) const {
    return binary(OP_LT, *this, y);
  }

  template<typename DataType>
  Matrix<DataType> Matrix<DataType>::zz_le(const Matrix<DataType> &y) const {
    return binary(OP_LE, *this, y);
  }

  template<typename DataType>
  Matrix<DataType> Matrix<DataType>::zz_eq(const Matrix<DataType> &y) const {
    return binary(OP_EQ, *this, y);
  }

  template<typename DataType>
  Matrix<DataType> Matrix<DataType>::zz_ne(const Matrix<DataType> &y) const {
    return binary(OP_NE, *this, y);
  }

  template<typename DataType>
  Matrix<DataType> Matrix<DataType>::__mrdivide__(const Matrix<DataType>& b) const {
    if (b.numel()==1) return *this/b;
    throw CasadiException("mrdivide: Not implemented");
  }

  template<typename DataType>
  Matrix<DataType> Matrix<DataType>::zz_mpower(const Matrix<DataType>& b) const {
    if (b.numel()==1) return pow(*this, b);
    throw CasadiException("mpower: Not implemented");
  }

  template<typename DataType>
  Sparsity& Matrix<DataType>::sparsityRef() {
    sparsity_.makeUnique(); // NOTE: Remove?
    return sparsity_;
  }

  template<typename DataType>
  void Matrix<DataType>::getBand(int kl, int ku, int ldres, DataType *res) const {
    // delete the content of the matrix
    for (int j=0; j<size1(); ++j) // loop over rows
      for (int s=0; s<kl+ku+1; ++s) // loop over the subdiagonals
        res[s + ldres*j] = 0;

    // loop over cols
    for (int i=0; i<size2(); ++i) {

      // loop over the non-zero elements
      for (int el=colind(i); el<colind(i+1); ++el) {
        int j=row(el);  // row

        // Check if we have not yet inside the band
        if (j<i-kl) continue;

        // Check if we are already outside the band
        if (j>i+ku) break;

        // Get the subdiagonal
        int s = i - j + ku;

        // Store the element
        res[s + ldres*j] = data()[el];
      }
    }
  }

  template<typename DataType>
  void Matrix<DataType>::set(DataType val, SparsityType sp) {
    std::fill(data().begin(), data().end(), val);
  }

  template<typename DataType>
  void Matrix<DataType>::get(DataType& val, SparsityType sp) const {
    getArray(&val, 1, SP_DENSE);
  }

  template<typename DataType>
  void Matrix<DataType>::set(const std::vector<DataType>& val, SparsityType sp) {
    setArray(val.empty() ? 0 : &val.front(), val.size(), sp);
  }

  template<typename DataType>
  void Matrix<DataType>::get(std::vector<DataType>& val, SparsityType sp) const {
    getArray(val.empty() ? 0 : &val.front(), val.size(), sp);
  }

  template<typename DataType>
  void Matrix<DataType>::set(const Matrix<DataType>& val, SparsityType sp) {
    sparsity().set(getPtr(data()), getPtr(val.data()), val.sparsity());
  }

  template<typename DataType>
  void Matrix<DataType>::setBV(const Matrix<DataType>& val) {
    bvec_t* bw_this = reinterpret_cast<bvec_t*>(getPtr(data()));
    const bvec_t* bw_val = reinterpret_cast<const bvec_t*>(getPtr(val.data()));
    sparsity().set(bw_this, bw_val, val.sparsity());
  }

  template<typename DataType>
  void Matrix<DataType>::setZeroBV() {
    bvec_t* bw_this = reinterpret_cast<bvec_t*>(getPtr(data()));
    std::fill(bw_this, bw_this+nnz(), bvec_t(0));
  }

  template<typename DataType>
  void Matrix<DataType>::borBV(const Matrix<DataType>& val) {
    bvec_t* bw_this = reinterpret_cast<bvec_t*>(getPtr(data()));
    const bvec_t* bw_val = reinterpret_cast<const bvec_t*>(getPtr(val.data()));
    sparsity().bor(bw_this, bw_val, val.sparsity());
  }

  template<typename DataType>
  void Matrix<DataType>::getArrayBV(bvec_t* val, int len) const {
    casadi_assert(len==nnz());
    const bvec_t* bw_this = reinterpret_cast<const bvec_t*>(getPtr(data()));
    std::copy(bw_this, bw_this+len, val);
  }

  template<typename DataType>
  void Matrix<DataType>::setArrayBV(const bvec_t* val, int len) {
    casadi_assert(len==nnz());
    bvec_t* bw_this = reinterpret_cast<bvec_t*>(getPtr(data()));
    std::copy(val, val+len, bw_this);
  }

  template<typename DataType>
  void Matrix<DataType>::borArrayBV(const bvec_t* val, int len) {
    casadi_assert(len==nnz());
    bvec_t* bw_this = reinterpret_cast<bvec_t*>(getPtr(data()));
    for (int i=0; i<len; ++i) *bw_this++ |= *val++;
  }

  template<typename DataType>
  void Matrix<DataType>::get(Matrix<DataType>& val, SparsityType sp) const {
    val.set(*this, sp);
  }

  template<typename DataType>
  void Matrix<DataType>::set(const DataType* val, SparsityType sp) {
    int len = sp==SP_SPARSE ? nnz() :
        sp==SP_DENSE || sp==SP_DENSETRANS ? numel() :
        sp==SP_SPARSESYM ? sizeU() : -1;
    setArray(val, len, sp);
  }

  template<typename DataType>
  void Matrix<DataType>::get(DataType* val, SparsityType sp) const {
    int len = sp==SP_SPARSE ? nnz() :
        sp==SP_DENSE || sp==SP_DENSETRANS ? numel() :
        sp==SP_SPARSESYM ? sizeU() : -1;
    getArray(val, len, sp);
  }

  template<typename DataType>
  void Matrix<DataType>::getArray(DataType* val, int len, SparsityType sp) const {
    // Get references to data for quick access
    const std::vector<DataType> &data = this->data();
    const int size1 = this->size1();
    const int size2 = this->size2();
    const int* colind = this->colindPtr();
    const int* row = this->rowPtr();

    if (sp==SP_SPARSE || (sp==SP_DENSE && isDense())) {
      casadi_assert_message(len==nnz(),
                            "Matrix<DataType>::getArray: Dimension mismatch." << std::endl <<
                            "Trying to fetch " << len << " elements from a " << dimString()
                            << " matrix with " << nnz() << " non-zeros.");
      copy(data.begin(), data.end(), val);
    } else if (sp==SP_DENSE) {
      casadi_assert_message(len==numel(),
                            "Matrix<DataType>::getArray: Dimension mismatch." << std::endl <<
                            "Trying to fetch " << len << " elements from a " << dimString()
                            << " matrix with " << numel() << " entries.");
      // Begin with all zeros
      std::fill(val, val+len, 0);

      // Add the nonzeros
      for (int cc=0; cc<size2; ++cc) { // loop over columns
        for (int el=colind[cc]; el<colind[cc+1]; ++el) { // loop over the non-zero elements
          int rr=row[el];
          val[rr+cc*size1] = data[el];
        }
      }
    } else if (sp==SP_DENSETRANS) {
      casadi_assert_message(len==numel(),
                            "Matrix<DataType>::getArray: Dimension mismatch." << std::endl <<
                            "Trying to fetch " << len << " elements from a " << dimString()
                            << " matrix with " << numel() << " entries.");
      // Begin with all zeros
      std::fill(val, val+len, 0);

      // Add the nonzeros
      for (int cc=0; cc<size2; ++cc) { // loop over columns
        for (int el=colind[cc]; el<colind[cc+1]; ++el) { // loop over the non-zero elements
          int rr=row[el];
          val[cc+rr*size2] = data[el];
        }
      }
    } else if (sp==SP_SPARSESYM) {
      // copy to the result vector
      int nz = 0;
      for (int cc=0; cc<size2; ++cc) {
        // Loop over the elements in the col
        for (int el=colind[cc]; el<colind[cc+1]; ++el) { // loop over the non-zero elements
          if (row[el] > cc) break; // break inner loop (only upper triangular part is used)
          val[nz++] = data[el];
        }
      }
    } else {
      casadi_error("Matrix<DataType>::getArray: not SPARSE, SP_SPARSESYM, DENSE or SP_DENSETRANS");
    }
  }

  /**
     Set stride to zero for unstrided acces
  */
  template<typename DataType>
  void Matrix<DataType>::getStridedArray(DataType* val, int len, int stride1, int stride2,
                                         SparsityType sp) const {
    if (stride1==0 || stride2==0 || (sp==SP_DENSE && stride2==1 && stride1==size1()))
        return getArray(val, len, sp);

    // Get references to data for quick access
    const std::vector<DataType> &data = this->data();
    //const int size1 = this->size1();
    const int size2 = this->size2();
    const int* colind = this->colindPtr();
    const int* row = this->rowPtr();

    if (sp==SP_SPARSE) {
      throw CasadiException("Matrix<DataType>::getArray: strided SPARSE not implemented");
    } else if (sp==SP_DENSE && isDense()) {
      for (int cc=0; cc<size2; ++cc) { // loop over columns
        for (int el=colind[cc]; el<colind[cc+1]; ++el) { // loop over the non-zero elements
          int rr=row[el];
          val[rr*stride2 + cc*stride1] = data[el];
        }
      }
    } else if (sp==SP_DENSETRANS && isDense()) {
      for (int cc=0; cc<size2; ++cc) { // loop over columns
        for (int el=colind[cc]; el<colind[cc+1]; ++el) { // loop over the non-zero elements
          int rr=row[el];
          val[cc*stride2 + rr*stride1] = data[el];
        }
      }
    } else if (sp==SP_DENSE) {
      casadi_error("Matrix<DataType>::getStridedArray: "
                   "strided sparse DMatrix->dense not implemented");
    } else if (sp==SP_SPARSESYM) {
      casadi_error("Matrix<DataType>::getStridedArray: strided SP_SPARSESYM not implemented");
    } else {
      casadi_error("Matrix<DataType>::getStridedArray: not SPARSE or DENSE");
    }

  }

  template<typename DataType>
  void Matrix<DataType>::setArray(const DataType* val, int len, SparsityType sp) {
    // Get references to data for quick access
    std::vector<DataType> &data = this->data();
    const int size1 = this->size1();
    const int size2 = this->size2();
    const int* colind = this->colindPtr();
    const int* row = this->rowPtr();

    if (sp==SP_SPARSE || (sp==SP_DENSE && numel()==nnz())) {
      casadi_assert_message(len==nnz(),
                            "Matrix<DataType>::setArray: Dimension mismatch." << std::endl <<
                            "Trying to pass " << len << " elements to a " << dimString()
                            << " matrix with " << nnz() << " non-zeros.");
      copy(val, val+len, data.begin());
    } else if (sp==SP_DENSE) {
      casadi_assert_message(len==numel(),
                            "Matrix<DataType>::setArray: Dimension mismatch." << std::endl <<
                            "Trying to pass " << len << " elements to a " << dimString()
                            << " matrix with " << numel() << " entries.");
      // Get the nonzeros
      for (int cc=0; cc<size2; ++cc) { // loop over columns
        for (int el=colind[cc]; el<colind[cc+1]; ++el) { // loop over the non-zero elements
          int rr=row[el];
          data[el] = val[rr+cc*size1];
        }
      }
    } else if (sp==SP_DENSETRANS) {
      casadi_assert_message(len==numel(),
                            "Matrix<DataType>::setArray: Dimension mismatch." << std::endl <<
                            "Trying to pass " << len << " elements to a " << dimString()
                            << " matrix with " << numel() << " entries.");
      // Get the nonzeros
      for (int cc=0; cc<size2; ++cc) { // loop over columns
        for (int el=colind[cc]; el<colind[cc+1]; ++el) { // loop over the non-zero elements
          int rr=row[el];
          data[el] = val[cc+rr*size2];
        }
      }
    } else if (sp==SP_SPARSESYM) {
      // NOTE: Has to be rewritten! sparsity().transpose(...) involves memory allocation
      // and is not threadsafe!!!
      // These routines is supposed to be used inside threadsafe code.
      std::vector<int> mapping;
      sparsity().transpose(mapping, false);
      // copy to the result vector
      int nz = 0;
      for (int cc=0; cc<size2; ++cc) {
        // Loop over the elements in the col
        for (int el=colind[cc]; el<colind[cc+1]; ++el) { // loop over the non-zero elements
          if (row[el] > cc) break; // break inner loop (only lower triangular part is used)
          data[mapping[el]] = data[el] = val[nz++];
        }
      }
    } else {
      throw CasadiException("Matrix<DataType>::setArray: "
                            "not SPARSE, SP_SPARSESYM, DENSE or SP_DENSETRANS");
    }
  }

  template<typename DataType>
  void Matrix<DataType>::getArray(DataType* val) const {
    getArray(val, nnz(), SP_SPARSE);
  }

  template<typename DataType>
  void Matrix<DataType>::setArray(const DataType* val) {
    setArray(val, nnz(), SP_SPARSE);
  }

  template<typename DataType>
  Matrix<DataType> Matrix<DataType>::zz_power(const Matrix<DataType>& y) const {
    return binary(OP_POW, *this, y);
  }

  template<typename DataType>
  Matrix<DataType> Matrix<DataType>::__constpow__(const Matrix<DataType>& y) const {
    return binary(OP_CONSTPOW, *this, y);
  }

  template<typename DataType>
  Matrix<DataType> Matrix<DataType>::zz_sin() const {
    return unary(OP_SIN, *this);
  }

  template<typename DataType>
  Matrix<DataType> Matrix<DataType>::zz_cos() const {
    return unary(OP_COS, *this);
  }

  template<typename DataType>
  Matrix<DataType> Matrix<DataType>::zz_tan() const {
    return unary(OP_TAN, *this);
  }

  template<typename DataType>
  Matrix<DataType> Matrix<DataType>::zz_erf() const {
    return unary(OP_ERF, *this);
  }

  template<typename DataType>
  Matrix<DataType> Matrix<DataType>::zz_asin() const {
    return unary(OP_ASIN, *this);
  }

  template<typename DataType>
  Matrix<DataType> Matrix<DataType>::zz_acos() const {
    return unary(OP_ACOS, *this);
  }

  template<typename DataType>
  Matrix<DataType> Matrix<DataType>::zz_atan() const {
    return unary(OP_ATAN, *this);
  }

  template<typename DataType>
  Matrix<DataType> Matrix<DataType>::zz_sinh() const {
    return unary(OP_SINH, *this);
  }

  template<typename DataType>
  Matrix<DataType> Matrix<DataType>::zz_cosh() const {
    return unary(OP_COSH, *this);
  }

  template<typename DataType>
  Matrix<DataType> Matrix<DataType>::zz_tanh() const {
    return unary(OP_TANH, *this);
  }

  template<typename DataType>
  Matrix<DataType> Matrix<DataType>::zz_asinh() const {
    return unary(OP_ASINH, *this);
  }

  template<typename DataType>
  Matrix<DataType> Matrix<DataType>::zz_acosh() const {
    return unary(OP_ACOSH, *this);
  }

  template<typename DataType>
  Matrix<DataType> Matrix<DataType>::zz_atanh() const {
    return unary(OP_ATANH, *this);
  }

  template<typename DataType>
  Matrix<DataType> Matrix<DataType>::zz_exp() const {
    return unary(OP_EXP, *this);
  }

  template<typename DataType>
  Matrix<DataType> Matrix<DataType>::zz_log() const {
    return unary(OP_LOG, *this);
  }

  template<typename DataType>
  Matrix<DataType> Matrix<DataType>::zz_log10() const {
    return log(*this)*(1/std::log(10.));
  }

  template<typename DataType>
  Matrix<DataType> Matrix<DataType>::zz_sqrt() const {
    return unary(OP_SQRT, *this);
  }

  template<typename DataType>
  Matrix<DataType> Matrix<DataType>::zz_floor() const {
    return unary(OP_FLOOR, *this);
  }

  template<typename DataType>
  Matrix<DataType> Matrix<DataType>::zz_ceil() const {
    return unary(OP_CEIL, *this);
  }

  template<typename DataType>
  Matrix<DataType> Matrix<DataType>::zz_mod(const Matrix<DataType>& y) const {
    return binary(OP_FMOD, *this, y);
  }

  template<typename DataType>
  Matrix<DataType> Matrix<DataType>::zz_abs() const {
    return unary(OP_FABS, *this);
  }

  template<typename DataType>
  Matrix<DataType> Matrix<DataType>::zz_sign() const {
    return unary(OP_SIGN, *this);
  }

  template<typename DataType>
  Matrix<DataType> Matrix<DataType>::__copysign__(const Matrix<DataType>& y) const {
    return binary(OP_COPYSIGN, *this, y);
  }

  template<typename DataType>
  Matrix<DataType> Matrix<DataType>::zz_erfinv() const {
    return unary(OP_ERFINV, *this);
  }

  template<typename DataType>
  Matrix<DataType> Matrix<DataType>::zz_min(const Matrix<DataType>& y) const {
    return binary(OP_FMIN, *this, y);
  }

  template<typename DataType>
  Matrix<DataType> Matrix<DataType>::zz_atan2(const Matrix<DataType>& y) const {
    return binary(OP_ATAN2, *this, y);
  }

  template<typename DataType>
  Matrix<DataType> Matrix<DataType>::zz_max(const Matrix<DataType>& y) const {
    return binary(OP_FMAX, *this, y);
  }

  template<typename DataType>
  Matrix<DataType> Matrix<DataType>::printme(const Matrix<DataType>& y) const {
    return binary(OP_PRINTME, *this, y);
  }

  template<typename DataType>
  Matrix<DataType> Matrix<DataType>::zz_not() const {
    return unary(OP_NOT, *this);
  }

  template<typename DataType>
  Matrix<DataType> Matrix<DataType>::zz_and(const Matrix<DataType>& y) const {
    return binary(OP_AND, *this, y);
  }

  template<typename DataType>
  Matrix<DataType> Matrix<DataType>::zz_or(const Matrix<DataType>& y) const {
    return binary(OP_OR, *this, y);
  }

  template<typename DataType>
  Matrix<DataType> Matrix<DataType>::zz_if_else_zero(const Matrix<DataType>& y) const {
    return binary(OP_IF_ELSE_ZERO, *this, y);
  }

  template<typename DataType>
  std::vector<DataType>& Matrix<DataType>::data() {
    return data_;
  }

  template<typename DataType>
  const std::vector<DataType>& Matrix<DataType>::data() const {
    return data_;
  }

  template<typename DataType>
  void Matrix<DataType>::erase(const std::vector<int>& rr, const std::vector<int>& cc, bool ind1) {
    // Erase from sparsity pattern
    std::vector<int> mapping = sparsityRef().erase(rr, cc, ind1);

    // Update non-zero entries
    for (int k=0; k<mapping.size(); ++k)
      data()[k] = data()[mapping[k]];

    // Truncate nonzero vector
    data().resize(mapping.size());
  }

  template<typename DataType>
  void Matrix<DataType>::erase(const std::vector<int>& rr, bool ind1) {
    // Erase from sparsity pattern
    std::vector<int> mapping = sparsityRef().erase(rr, ind1);

    // Update non-zero entries
    for (int k=0; k<mapping.size(); ++k)
      data()[k] = data()[mapping[k]];

    // Truncate nonzero vector
    data().resize(mapping.size());
  }

  template<typename DataType>
  void Matrix<DataType>::remove(const std::vector<int>& rr, const std::vector<int>& cc) {
    if (!inBounds(rr, size1())) {
      casadi_error("Remove(rr, cc) out of bounds. Your rr contains "
                   << *std::min_element(rr.begin(), rr.end()) << " up to "
                   << *std::max_element(rr.begin(), rr.end())
                   << ", which is outside of the matrix shape " << dimString() << ".");
    }
    if (!inBounds(cc, size2())) {
      casadi_error("Remove(rr, cc) out of bounds. Your cc contains "
                   << *std::min_element(cc.begin(), cc.end()) << " up to "
                   << *std::max_element(cc.begin(), cc.end())
                   << ", which is outside of the matrix shape " << dimString() << ".");
    }

    // Remove by performing a complementary slice
    std::vector<int> rrc = complement(rr, size1());
    std::vector<int> ccc = complement(cc, size2());

    Matrix<DataType> ret = (*this)(rrc, ccc);

    operator=(ret);

  }

  template<typename DataType>
  void Matrix<DataType>::enlarge(int nrow, int ncol, const std::vector<int>& rr,
                                 const std::vector<int>& cc, bool ind1) {
    sparsityRef().enlarge(nrow, ncol, rr, cc, ind1);
  }

  template<typename DataType>
  void Matrix<DataType>::sanityCheck(bool complete) const {
    sparsity_.sanityCheck(complete);
    if (data_.size()!=sparsity_.nnz()) {
      std::stringstream s;
      s << "Matrix is not sane. The following must hold:" << std::endl;
      s << "  data().size() = sparsity().nnz(), but got data().size()  = " << data_.size()
        << "   and sparsity().nnz() = "  << sparsity_.nnz() << std::endl;
      casadi_error(s.str());
    }
  }

  template<typename DataType>
  Matrix<DataType> Matrix<DataType>::zz_mtimes(const Matrix<DataType> &y) const {
    Matrix<DataType> z = Matrix<DataType>::zeros(mul(sparsity(), y.sparsity()));
    return zz_mtimes(y, z);
  }

  template<typename DataType>
  Matrix<DataType> Matrix<DataType>::zz_mtimes(const Matrix<DataType> &y,
                                               const Matrix<DataType> &z) const {
    if (isScalar() || y.isScalar()) {
      // Use element-wise multiplication if at least one factor scalar
      return z + *this*y;
    }

    // Check matching dimensions
    casadi_assert_message(size2()==y.size1(),
                          "Matrix product with incompatible dimensions. Lhs is "
                          << dimString() << " and rhs is " << y.dimString() << ".");

    // Check if we can simplify the product
    if (isIdentity()) {
      return y + z;
    } else if (y.isIdentity()) {
      return *this + z;
    } else if (isZero() || y.isZero()) {
      return z;
    } else {
      // Carry out the matrix product
      Matrix<DataType> ret = z;
      std::vector<DataType> work(size1());
      mul_no_alloc(*this, y, ret, work);
      return ret;
    }
  }

  template<typename DataType>
  void Matrix<DataType>::mul_no_alloc(const Matrix<DataType> &x, const Matrix<DataType> &y,
                                      Matrix<DataType>& z, std::vector<DataType>& work,
                                      bool transpose_x) {

    // Assert dimensions
    if (transpose_x) {
      casadi_assert_message(z.size1()==x.size2() && x.size1()==y.size1() && y.size2()==z.size2(),
                            "Dimension error. Got x=" << x.dimString() << ", y=" << y.dimString()
                            << " and z=" << z.dimString() << ".");
      casadi_assert_message(work.size()>=y.size1(),
                            "Work vector too small: " << work.size() << " < " << y.size1());
    } else {
      casadi_assert_message(z.size1()==x.size1() && x.size2()==y.size1() && y.size2()==z.size2(),
                            "Dimension error. Got x=" << x.dimString() << ", y=" << y.dimString()
                            << " and z=" << z.dimString() << ".");
      casadi_assert_message(work.size()>=z.size1(),
                            "Work vector too small: " << work.size() << " < " << z.size1());
    }

    // Direct access to the arrays
    const int* y_colind = y.colindPtr();
    const int* y_row = y.rowPtr();
    const std::vector<DataType> &y_data = y.data();
    const int* x_colind = x.colindPtr();
    const int* x_row = x.rowPtr();
    const std::vector<DataType> &x_data = x.data();
    const int* z_colind = z.colindPtr();
    const int* z_row = z.rowPtr();
    std::vector<DataType> &z_data = z.data();

    // Loop over the columns of y and z
    int ncol = z.size2();
    for (int cc=0; cc<ncol; ++cc) {
      if (transpose_x) { // Transposed variant, loop over z

        // Get the dense column of y
        for (int kk=y_colind[cc]; kk<y_colind[cc+1]; ++kk) {
          work[y_row[kk]] = y_data[kk];
        }

        // Loop over the nonzeros of z
        for (int kk=z_colind[cc]; kk<z_colind[cc+1]; ++kk) {
          int rr = z_row[kk];

          // Loop over corresponding columns of x
          for (int kk1=x_colind[rr]; kk1<x_colind[rr+1]; ++kk1) {
            z_data[kk] += x_data[kk1] * work[x_row[kk1]];
          }
        }

      } else { // Non-transposed variant, loop over y

        // Get the dense column of z
        for (int kk=z_colind[cc]; kk<z_colind[cc+1]; ++kk) {
          work[z_row[kk]] = z_data[kk];
        }

        // Loop over the nonzeros of y
        for (int kk=y_colind[cc]; kk<y_colind[cc+1]; ++kk) {
          int rr = y_row[kk];

          // Loop over corresponding columns of x
          for (int kk1=x_colind[rr]; kk1<x_colind[rr+1]; ++kk1) {
            work[x_row[kk1]] += x_data[kk1] * y_data[kk];
          }
        }

        // Get the sparse column of z
        for (int kk=z_colind[cc]; kk<z_colind[cc+1]; ++kk) {
          z_data[kk] = work[z_row[kk]];
        }
      }
    }
  }

  template<typename DataType>
  void Matrix<DataType>::mul_no_alloc_nn(const Matrix<DataType> &x, const Matrix<DataType> &y,
                                         Matrix<DataType>& z) {
    // Assert dimensions
    casadi_assert_message(x.size1()==z.size1(), "Dimension error. Got x=" << x.dimString()
                          << " and z=" << z.dimString() << ".");
    casadi_assert_message(y.size2()==z.size2(), "Dimension error. Got y=" << y.dimString()
                          << " and z=" << z.dimString() << ".");
    casadi_assert_message(y.size1()==x.size2(), "Dimension error. Got y=" << y.dimString()
                          << " and x=" << x.dimString() << ".");

    // Direct access to the arrays
    int y_ncol = y.size2();
    const int* y_colind = y.colindPtr();
    const int* y_row = y.rowPtr();
    const std::vector<DataType> &y_data = y.data();
    const int* x_colind = x.colindPtr();
    const int* x_row = x.rowPtr();
    const std::vector<DataType> &x_data = x.data();
    const int* z_colind = z.colindPtr();
    const int* z_row = z.rowPtr();
    std::vector<DataType> &z_data = z.data();

    // loop over the cols of the first argument
    for (int i=0; i<y_ncol; ++i) {
      // loop over the non-zeros of the first argument
      for (int el=y_colind[i]; el<y_colind[i+1]; ++el) {
        int j = y_row[el];
        int el1 = z_colind[i];
        int el2 = x_colind[j];
        while (el1 < z_colind[i+1] && el2 < x_colind[j+1]) { // loop over matching non-zero elements
          int j1 = z_row[el1];
          int i2 = x_row[el2];
          if (j1==i2) {
            z_data[el1++] += y_data[el]*x_data[el2++];
          } else if (j1<i2) {
            el1++;
          } else {
            el2++;
          }
        }
      }
    }
  }

  template<typename DataType>
  void Matrix<DataType>::mul_no_alloc(const Matrix<DataType> &x,
                                      const std::vector<DataType> &y,
                                      std::vector<DataType>& z, bool transpose_x) {
    // Assert dimensions
    casadi_assert_message(z.size()==transpose_x ? x.size2() : x.size1(),
                          "Dimension error. Got transpose_x=" << transpose_x
                          << ", x=" << x.dimString() << " and z=" << z.size() << ".");
    casadi_assert_message(y.size()==transpose_x ? x.size1() : x.size2(),
                          "Dimension error. Got transpose_x=" << transpose_x
                          << ", x=" << x.dimString() << " and y=" << y.size() << ".");

    // Direct access to the arrays
    int x_ncol = x.size2();
    const int* x_colind = x.colindPtr();
    const int* x_row = x.rowPtr();
    const std::vector<DataType> &x_data = x.data();

    // loop over the columns of the matrix
    for (int i=0; i<x_ncol; ++i) {
      // loop over the non-zeros of the matrix
      for (int el=x_colind[i]; el<x_colind[i+1]; ++el) {
        int j = x_row[el];
        if (transpose_x) {
          z[i] += x_data[el] * y[j];
        } else {
          z[j] += x_data[el] * y[i];
        }
      }
    }
  }

  template<typename DataType>
  void Matrix<DataType>::mul_no_alloc_tn(const Matrix<DataType> &x_trans,
                                         const std::vector<DataType> &y,
                                         std::vector<DataType>& z) {
    mul_no_alloc(x_trans, y, z, true);
  }

  template<typename DataType>
  void Matrix<DataType>::mul_no_alloc_nn(const Matrix<DataType>& x,
                                         const std::vector<DataType> &y,
                                         std::vector<DataType> &z) {
    mul_no_alloc(x, y, z);
  }

  template<typename DataType>
  void Matrix<DataType>::mul_no_alloc_nt(const Matrix<DataType> &x,
                                         const Matrix<DataType>& y_trans,
                                         Matrix<DataType> &z) {
    // Assert dimensions
    casadi_assert_message(y_trans.size1()==z.size2(), "Dimension error. Got y_trans="
                          << y_trans.dimString() << " and z=" << z.dimString() << ".");
    casadi_assert_message(x.size1()==z.size1(), "Dimension error. Got x=" << x.dimString()
                          << " and z=" << z.dimString() << ".");
    casadi_assert_message(y_trans.size2()==x.size2(), "Dimension error. Got y_trans="
                          << y_trans.dimString() << " and x=" << x.dimString() << ".");

    // Direct access to the arrays
    const int* y_rowind = y_trans.colindPtr();
    int y_nrow = y_trans.size2();
    const int* y_col = y_trans.rowPtr();
    const std::vector<DataType> &y_trans_data = y_trans.data();
    const int* x_colind = x.colindPtr();
    const int* x_row = x.rowPtr();
    const std::vector<DataType> &x_data = x.data();
    const int* z_colind = z.colindPtr();
    const int* z_row = z.rowPtr();
    std::vector<DataType> &z_data = z.data();

    // loop over the rows of the first argument
    for (int i=0; i<y_nrow; ++i) {
      // loop over the non-zeros of the first argument
      for (int el=y_rowind[i]; el<y_rowind[i+1]; ++el) {
        int j = y_col[el];
        int el1 = x_colind[i];
        int el2 = z_colind[j];
        while (el1 < x_colind[i+1] && el2 < z_colind[j+1]) { // loop over matching non-zero elements
          int j1 = x_row[el1];
          int i2 = z_row[el2];
          if (j1==i2) {
            z_data[el2++] += y_trans_data[el] * x_data[el1++];
          } else if (j1<i2) {
            el1++;
          } else {
            el2++;
          }
        }
      }
    }
  }

  template<typename DataType>
  void Matrix<DataType>::mul_no_alloc_tn(const Matrix<DataType> &x_trans,
                                         const Matrix<DataType> &y,
                                         Matrix<DataType>& z) {
    // Assert dimensions
    casadi_assert_message(y.size2()==z.size2(), "Dimension error. Got y="
                          << y.dimString() << " and z=" << z.dimString() << ".");
    casadi_assert_message(x_trans.size2()==z.size1(), "Dimension error. Got x_trans="
                          << x_trans.dimString() << " and z=" << z.dimString() << ".");
    casadi_assert_message(y.size1()==x_trans.size1(), "Dimension error. Got y="
                          << y.dimString() << " and x_trans=" << x_trans.dimString() << ".");

    // Direct access to the arrays
    const int* y_colind = y.colindPtr();
    const int* y_row = y.rowPtr();
    const std::vector<DataType> &y_data = y.data();
    const int* x_rowind = x_trans.colindPtr();
    const int* x_col = x_trans.rowPtr();
    const std::vector<DataType> &x_trans_data = x_trans.data();
    int z_ncol = z.size2();
    const int* z_colind = z.colindPtr();
    const int* z_row = z.rowPtr();
    std::vector<DataType> &z_data = z.data();

    // loop over the cols of the resulting matrix
    for (int i=0; i<z_ncol; ++i) {
      // loop over the non-zeros of the resulting matrix
      for (int el=z_colind[i]; el<z_colind[i+1]; ++el) {
        int j = z_row[el];
        int el1 = y_colind[i];
        int el2 = x_rowind[j];
        while (el1 < y_colind[i+1] && el2 < x_rowind[j+1]) { // loop over non-zero elements
          int j1 = y_row[el1];
          int i2 = x_col[el2];
          if (j1==i2) {
            z_data[el] += y_data[el1++] * x_trans_data[el2++];
          } else if (j1<i2) {
            el1++;
          } else {
            el2++;
          }
        }
      }
    }
  }

  template<typename DataType>
  template<bool Fwd>
  void Matrix<DataType>::mul_sparsity(Matrix<DataType> &x,
                                      Matrix<DataType> &y,
                                      Matrix<DataType>& z,
                                      std::vector<DataType>& work) {

    // Assert dimensions
    casadi_assert_message(z.size1()==x.size1() && x.size2()==y.size1() && y.size2()==z.size2(),
                          "Dimension error. Got x=" << x.dimString() << ", y=" << y.dimString()
                          << " and z=" << z.dimString() << ".");

    // Make sure work vector large enough
    casadi_assert_message(work.size()>=z.size1(),
                          "Work vector too small: " << work.size() << " < " << z.size1());

    // Direct access to the arrays
    const int* y_colind = y.colindPtr();
    const int* y_row = y.rowPtr();
    const int* x_colind = x.colindPtr();
    const int* x_row = x.rowPtr();
    const int* z_colind = z.colindPtr();
    const int* z_row = z.rowPtr();

    // Convert data array to arrays of integers
    bvec_t *y_data = get_bvec_t(y.data());
    bvec_t *x_data = get_bvec_t(x.data());
    bvec_t *z_data = get_bvec_t(z.data());
    bvec_t *w = get_bvec_t(work);

    // Loop over the columns of y and z
    int ncol = z.size2();
    for (int cc=0; cc<ncol; ++cc) {
      // Get the dense column of z
      for (int kk=z_colind[cc]; kk<z_colind[cc+1]; ++kk) {
        w[z_row[kk]] = z_data[kk];
      }

      // Loop over the nonzeros of y
      for (int kk=y_colind[cc]; kk<y_colind[cc+1]; ++kk) {
        int rr = y_row[kk];

        // Loop over corresponding columns of x
        if (Fwd) {
          bvec_t yy = y_data[kk];
          for (int kk1=x_colind[rr]; kk1<x_colind[rr+1]; ++kk1) {
            w[x_row[kk1]] |= x_data[kk1] | yy;
          }
        } else {
          bvec_t yy = 0;
          for (int kk1=x_colind[rr]; kk1<x_colind[rr+1]; ++kk1) {
            yy |= w[x_row[kk1]];
            x_data[kk1] |= w[x_row[kk1]];
          }
          y_data[kk] |= yy;
        }
      }

      // Get the sparse column of z
      for (int kk=z_colind[cc]; kk<z_colind[cc+1]; ++kk) {
        z_data[kk] = w[z_row[kk]];
      }
    }
  }

  template<typename DataType>
  DataType Matrix<DataType>::quad_form(const std::vector<DataType>& x, const Matrix<DataType>& A) {
    // Assert dimensions
    casadi_assert_message(x.size()==A.size2() && x.size()==A.size1(),
                          "Dimension mismatch. Got x=" << x.size() << " and A=" << A.dimString());

    // Access the internal data of A
    const int* A_colind = A.colindPtr();
    const int* A_row = A.rowPtr();
    const std::vector<DataType> &A_data = A.data();

    // Return value
    DataType ret=0;

    // Loop over the cols of A
    for (int i=0; i<x.size(); ++i) {
      // Loop over the nonzeros of A
      for (int el=A_colind[i]; el<A_colind[i+1]; ++el) {
        // Get row
        int j = A_row[el];

        // Add contribution
        ret += x[i]*A_data[el]*x[j];
      }
    }

    return ret;
  }

  template<typename DataType>
  Matrix<DataType> Matrix<DataType>::T() const {
    // quick return if empty or scalar
    if ((size1()==0 && size2()==0) || isScalar()) return *this;

    // Create the new sparsity pattern and the mapping
    std::vector<int> mapping;
    Sparsity s = sparsity().transpose(mapping);

    // create the return matrix
    Matrix<DataType> ret(s);

    // Copy the content
    for (int i=0; i<mapping.size(); ++i)
      ret.at(i) = at(mapping[i]);

    return ret;
  }

  template<typename DataType>
  const DataType Matrix<DataType>::toScalar() const {
    // Make sure that the matrix is 1-by-1
    casadi_assert_message(isScalar(), "Can only convert 1-by-1 matrices to scalars");

    // return zero or the nonzero element
    if (nnz()==1)
      return data()[0];
    else
      return casadi_limits<DataType>::zero;
  }

  template<typename DataType>
  Matrix<DataType> Matrix<DataType>::binary(int op,
                                            const Matrix<DataType> &x,
                                            const Matrix<DataType> &y) {
    if (x.numel()==1)
      return scalar_matrix(op, x, y);
    else if (y.numel()==1)
      return matrix_scalar(op, x, y);
    else
      return matrix_matrix(op, x, y);
  }

  template<typename DataType>
  Matrix<DataType> Matrix<DataType>::scalar_matrix(int op,
                                                   const Matrix<DataType> &x,
                                                   const Matrix<DataType> &y) {
    // Return value
    Matrix<DataType> ret(y.sparsity());

    // Nonzeros
    std::vector<DataType>& ret_data = ret.data();
    const std::vector<DataType>& x_data = x.data();
    const DataType& x_val = x_data.empty() ? casadi_limits<DataType>::zero : x.front();
    const std::vector<DataType>& y_data = y.data();

    // Do the operation on all non-zero elements
    for (int el=0; el<y.nnz(); ++el) {
      casadi_math<DataType>::fun(op, x_val, y_data[el], ret_data[el]);
    }

    // Check the value of the structural zero-entries, if there are any
    if (!y.isDense() && !operation_checker<Function0Checker>(op)) {
      // Get the value for the structural zeros
      DataType fcn_0;
      casadi_math<DataType>::fun(op, x_val, casadi_limits<DataType>::zero, fcn_0);
      if (!casadi_limits<DataType>::isZero(fcn_0)) { // Remove this if?
        ret.makeDense(fcn_0);
      }
    }

    return ret;
  }

  template<typename DataType>
  Matrix<DataType> Matrix<DataType>::matrix_scalar(int op,
                                                   const Matrix<DataType> &x,
                                                   const Matrix<DataType> &y) {
    // Return value
    Matrix<DataType> ret(x.sparsity());

    // Nonzeros
    std::vector<DataType>& ret_data = ret.data();
    const std::vector<DataType>& x_data = x.data();
    const std::vector<DataType>& y_data = y.data();
    const DataType& y_val = y_data.empty() ? casadi_limits<DataType>::zero : y.front();

    // Do the operation on all non-zero elements
    for (int el=0; el<x.nnz(); ++el) {
      casadi_math<DataType>::fun(op, x_data[el], y_val, ret_data[el]);
    }

    // Check the value of the structural zero-entries, if there are any
    if (!x.isDense() && !operation_checker<F0XChecker>(op)) {
      // Get the value for the structural zeros
      DataType fcn_0;
      casadi_math<DataType>::fun(op, casadi_limits<DataType>::zero, y_val, fcn_0);
      if (!casadi_limits<DataType>::isZero(fcn_0)) { // Remove this if?
        ret.makeDense(fcn_0);
      }
    }

    return ret;
  }

  template<typename DataType>
  Matrix<DataType> Matrix<DataType>::matrix_matrix(int op,
                                                   const Matrix<DataType> &x,
                                                   const Matrix<DataType> &y) {

    if (!(x.size2() == y.size2() && x.size1() == y.size1())) {
      std::stringstream ss;
      casadi_math<DataType>::print(op, ss, "lhs", "rhs");
      casadi_error("matrix_matrix: dimension mismatch in element-wise matrix operation "
                   << ss.str() <<"." << std::endl << "Left argument has shape " << x.dimString()
                   << ", right has shape " << y.dimString() << ". They should be equal.");
    }

    // Get the sparsity pattern of the result
    // (ignoring structural zeros giving rise to nonzero result)
    const Sparsity& x_sp = x.sparsity();
    const Sparsity& y_sp = y.sparsity();
    Sparsity r_sp = x_sp.patternCombine(y_sp,
                                        operation_checker<F0XChecker>(op),
                                        operation_checker<Function0Checker>(op));

    // Return value
    Matrix<DataType> r(r_sp);

    // Perform the operations elementwise
    if (x_sp==y_sp) {
      // Matching sparsities
      casadi_math<DataType>::fun(op, getPtr(x.data()), getPtr(y.data()),
                                 getPtr(r.data()), r_sp.nnz());
    } else if (y_sp==r_sp) {
      // Project first argument
      Matrix<DataType> x_mod = x(r_sp);
      casadi_math<DataType>::fun(op, getPtr(x_mod.data()), getPtr(y.data()),
                                 getPtr(r.data()), r_sp.nnz());
    } else if (x_sp==r_sp) {
      // Project second argument
      Matrix<DataType> y_mod = y(r_sp);
      casadi_math<DataType>::fun(op, getPtr(x.data()),
                                 getPtr(y_mod.data()), getPtr(r.data()), r_sp.nnz());
    } else {
      // Project both arguments
      Matrix<DataType> x_mod = x(r_sp);
      Matrix<DataType> y_mod = y(r_sp);
      casadi_math<DataType>::fun(op, getPtr(x_mod.data()), getPtr(y_mod.data()),
                                 getPtr(r.data()), r_sp.nnz());
    }

    // Handle structural zeros giving rise to nonzero result, e.g. cos(0) == 1
    if (!r.isDense() && !operation_checker<F00Checker>(op)) {
      // Get the value for the structural zeros
      DataType fcn_0;
      casadi_math<DataType>::fun(op, casadi_limits<DataType>::zero,
                                 casadi_limits<DataType>::zero, fcn_0);
      r.makeDense(fcn_0);
    }

    return r;
  }

  template<typename DataType>
  Matrix<DataType> Matrix<DataType>::triplet(const std::vector<int>& row,
                                             const std::vector<int>& col,
                                             const std::vector<DataType>& d) {
    return triplet(row, col, d, *std::max_element(row.begin(), row.end()),
                   *std::max_element(col.begin(), col.end()));
  }

  template<typename DataType>
  Matrix<DataType> Matrix<DataType>::triplet(const std::vector<int>& row,
                                             const std::vector<int>& col,
                                             const std::vector<DataType>& d,
                                             const std::pair<int, int>& rc) {
    return triplet(row, col, d, rc.first, rc.second);
  }

  template<typename DataType>
  Matrix<DataType> Matrix<DataType>::triplet(const std::vector<int>& row,
                                             const std::vector<int>& col,
                                             const std::vector<DataType>& d,
                                             int nrow, int ncol) {
    casadi_assert_message(col.size()==row.size() && col.size()==d.size(),
                          "Argument error in Matrix<DataType>::triplet(row, col, d): "
                          "supplied lists must all be of equal length, but got: "
                          << row.size() << ", " << col.size()  << " and " << d.size());
    std::vector<int> mapping;
    Sparsity sp = Sparsity::triplet(nrow, ncol, row, col, mapping);
    std::vector<DataType> v(mapping.size());
    for (int k=0; k<v.size(); ++k) v[k] = d[mapping[k]];
    return Matrix<DataType>(sp, v, false);
  }

  template<typename DataType>
  Matrix<DataType> Matrix<DataType>::eye(int n) {
    return Matrix<DataType>::ones(Sparsity::diag(n));
  }

  template<typename DataType>
  Matrix<DataType> Matrix<DataType>::inf(const Sparsity& sp) {
    casadi_assert_message(std::numeric_limits<DataType>::has_infinity,
                          "Datatype cannot represent infinity");
    return Matrix<DataType>(sp, std::numeric_limits<DataType>::infinity(), false);
  }


  template<typename DataType>
  Matrix<DataType> Matrix<DataType>::inf(const std::pair<int, int>& rc) {
    return inf(rc.first, rc.second);
  }

  template<typename DataType>
  Matrix<DataType> Matrix<DataType>::inf(int nrow, int ncol) {
    return inf(Sparsity::dense(nrow, ncol));
  }

  template<typename DataType>
  Matrix<DataType> Matrix<DataType>::nan(const Sparsity& sp) {
    casadi_assert_message(std::numeric_limits<DataType>::has_quiet_NaN,
                          "Datatype cannot represent not-a-number");
    return Matrix<DataType>(sp, std::numeric_limits<DataType>::quiet_NaN(), false);
  }

  template<typename DataType>
  Matrix<DataType> Matrix<DataType>::nan(const std::pair<int, int>& rc) {
    return nan(rc.first, rc.second);
  }

  template<typename DataType>
  Matrix<DataType> Matrix<DataType>::nan(int nrow, int ncol) {
    return nan(Sparsity::dense(nrow, ncol));
  }

  template<typename DataType>
  void Matrix<DataType>::append(const Matrix<DataType>& y) {
    // Quick return if expr is empty
    if (size2()==0 && size1()==0) {
      *this=y;
      return;
    }

    // Quick return if empty
    if (y.size2()==0 && y.size1()==0) return;

    // Appending can be done efficiently if vectors
    if (isVector()) {
      // Append the sparsity pattern vertically
      sparsityRef().append(y.sparsity());

      // Add the non-zeros at the end
      data().insert(end(), y.begin(), y.end());
    } else {
      // Fall back on vertical concatenation
      *this = vertcat(*this, y);
    }
  }

  template<typename DataType>
  void Matrix<DataType>::appendColumns(const Matrix<DataType>& y) {

    // Quick return if expr is empty
    if (size2()==0 && size1()==0) {
      *this=y;
      return;
    }

    // Quick return if empty
    if (y.size2()==0 && y.size1()==0) return;

    // Append the sparsity pattern
    sparsityRef().appendColumns(y.sparsity());

    // Add the non-zeros at the end
    data().insert(end(), y.begin(), y.end());
  }

  template<typename DataType>
  NonZeroIterator<DataType>::NonZeroIterator(const Matrix<DataType> & m)
    : m_(m) {
    nz.i  = 0;
    nz.j  = 0;
    nz.k  = 0;
  }

  template<typename DataType>
  bool NonZeroIterator<DataType>::operator==(const NonZeroIterator<DataType>& rhs)
  { return (m_ == rhs.m_) && (nz.k==rhs.nz.k); }

  template<typename DataType>
  NonZero<DataType>& NonZeroIterator<DataType>::operator*() {
    return nz;
  }

  template<typename DataType>
  NonZeroIterator<DataType>& NonZeroIterator<DataType>::operator++() {
    nz.k ++;

    if (nz.k < m_.size()) {
      nz.j = m_.row()[nz.k];
      nz.el = m_.data()[nz.k];
      while (nz.k>=m_.colind(nz.i)) {nz.i++; }
    }
    return *this;
  }

  template<typename DataType>
  NonZeroIterator<DataType> NonZeroIterator<DataType>::begin() {
    NonZeroIterator<DataType> it = NonZeroIterator<DataType>(m_);
    return it;

  }
  template<typename DataType>
  NonZeroIterator<DataType> NonZeroIterator<DataType>::end() {
    NonZeroIterator<DataType> it = NonZeroIterator<DataType>(m_);
    it.nz.k = m_.size()-1;
    return it;
  }

  template<typename DataType>
  bool Matrix<DataType>::isRegular() const {
    return casadi::isRegular(data_);
  }

  template<typename DataType>
  bool Matrix<DataType>::isSmooth() const {
    return true;
  }

  template<typename DataType>
  long Matrix<DataType>::getElementHash() const {
    throw CasadiException("\"getElementHash\" not defined for instantiation");
  }

  template<typename DataType>
  bool Matrix<DataType>::isLeaf() const {
    throw CasadiException("\"isLeaf\" not defined for instantiation");
  }

  template<typename DataType>
  bool Matrix<DataType>::isCommutative() const {
    throw CasadiException("\"isCommutative\" not defined for instantiation");
  }

  template<typename DataType>
  bool Matrix<DataType>::isSymbolic() const {
    return false;
  }

  template<typename DataType>
  bool Matrix<DataType>::isSymbolicSparse() const {
    return false;
  }

  template<typename DataType>
  bool Matrix<DataType>::isInteger() const {
    // loop over non-zero elements
    for (int k=0; k<nnz(); ++k)
      if (!casadi_limits<DataType>::isInteger(at(k))) // if an element is not integer
        return false;

    // Integer if reached this point
    return true;
  }

  template<typename DataType>
  bool Matrix<DataType>::isConstant() const {
    // loop over non-zero elements
    for (int k=0; k<nnz(); ++k)
      if (!casadi_limits<DataType>::isConstant(at(k))) // if an element is not constant
        return false;

    // Constant if we reach this point
    return true;
  }

  template<typename DataType>
  bool Matrix<DataType>::isOne() const {
    if (!isDense()) {
      return false;
    }

    // loop over non-zero elements
    for (int el=0; el<nnz(); ++el)
      if (!casadi_limits<DataType>::isOne(at(el)))
        return false;

    return true;
  }

  template<typename DataType>
  bool Matrix<DataType>::isMinusOne() const {
    if (!isDense()) {
      return false;
    }

    // loop over non-zero elements
    for (int el=0; el<nnz(); ++el)
      if (!casadi_limits<DataType>::isMinusOne(at(el)))
        return false;

    return true;
  }

  template<typename DataType>
  bool Matrix<DataType>::isZero() const {

    // loop over (potentially) non-zero elements
    for (int el=0; el<nnz(); ++el)
      if (!casadi_limits<DataType>::isZero(at(el)))
        return false;

    return true;
  }

  template<typename DataType>
  bool Matrix<DataType>::isIdentity() const {

    // Make sure that the matrix is diagonal
    if (!sparsity().isDiagonal())
      return false;

    // Make sure that all entries are one
    for (typename Matrix<DataType>::const_iterator it=begin(); it!=end(); ++it) {
      if (!casadi_limits<DataType>::isOne(*it))
        return false;
    }

    return true;
  }

  template<typename DataType>
  bool Matrix<DataType>::zz_isEqual(const Matrix<DataType> &ex2, int depth) const {
    // Assert matching dimensions
    casadi_assert_message(shape() == ex2.shape(), "Dimension mismatch");

    // Project to union of patterns and call recursively if different sparsity
    if (sparsity() != ex2.sparsity()) {
      Sparsity sp = sparsity() + ex2.sparsity();
      return isEqual(setSparse(sp), ex2.setSparse(sp), depth);
    }

    // Check individual elements
    for (int k=0; k<nnz(); ++k) {
      if (!isEqual(at(k), ex2.at(k), depth)) return false;
    }

    // True if reched this point
    return true;
  }

  template<typename DataType>
  bool Matrix<DataType>::hasNonStructuralZeros() const {
    // Check if the structural nonzero is known to be zero
    for (int el=0; el<nnz(); ++el) {
      if (casadi_limits<DataType>::isZero(at(el)))
        return true;
    }

    // No known zeros amongst the structurally nonzero entries
    return false;
  }

  template<typename DataType>
  double Matrix<DataType>::getValue() const {
    return static_cast<double>(toScalar());
  }

  template<typename DataType>
  std::string Matrix<DataType>::getName() const {
    throw CasadiException("\"getName\" not defined for instantiation");
  }

  template<typename DataType>
  Matrix<DataType> Matrix<DataType>::getDep(int ch) const {
    throw CasadiException("\"getDep\" not defined for instantiation");
  }

  template<typename DataType>
  int Matrix<DataType>::getNdeps() const {
    throw CasadiException("\"getNdeps\" not defined for instantiation");
  }

  template<typename DataType>
  Matrix<DataType> Matrix<DataType>::setSparse(const Sparsity& sp, bool intersect) const {
    if (intersect) {
      return setSparse(sp.patternIntersection(sparsity()), false);
    } else {
      Matrix<DataType> ret(sp);
      ret.set(*this);
      return ret;
    }
  }

  template<typename DataType>
  void Matrix<DataType>::setMaxNumCallsInPrint(long num) {
    throw CasadiException("\"setMaxNumCallsInPrint\" not defined for instantiation");
  }

  template<typename DataType>
  long Matrix<DataType>::getMaxNumCallsInPrint() {
    throw CasadiException("\"getMaxNumCallsInPrint\" not defined for instantiation");
  }

  template<typename DataType>
  void Matrix<DataType>::setEqualityCheckingDepth(int eq_depth) {
    throw CasadiException("\"setEqualityCheckingDepth\" not defined for instantiation");
  }

  template<typename DataType>
  int Matrix<DataType>::getEqualityCheckingDepth() {
    throw CasadiException("\"getEqualityCheckingDepth\" not defined for instantiation");
  }

  template<typename DataType>
  Matrix<DataType> Matrix<DataType>::zz_det() const {
    int n = size2();
    casadi_assert_message(n == size1(), "matrix must be square");

    // Trivial return if scalar
    if (isScalar()) return toScalar();

    // Trivial case 2 x 2
    if (n==2) return elem(0, 0) * elem(1, 1) - elem(0, 1) * elem(1, 0);

    // Return expression
    Matrix<DataType> ret = 0;

    // Find out which is the best direction to expand along

    // Build up an IMatrix with ones on the non-zeros
    Matrix<int> sp = IMatrix::ones(sparsity());

    // Have a count of the nonzeros for each row
    Matrix<int> row_count = sp.zz_sumCols();

    // A blank row? determinant is structurally zero
    if (!row_count.isDense()) return 0;

    // Have a count of the nonzeros for each col
    Matrix<int> col_count = sp.zz_sumRows().T();

    // A blank col? determinant is structurally zero
    if (!row_count.isDense()) return 0;

    int min_row = std::distance(row_count.data().begin(),
                                std::min_element(row_count.data().begin(),
                                                 row_count.data().end()));
    int min_col = std::distance(col_count.data().begin(),
                                std::min_element(col_count.data().begin(),
                                                 col_count.data().end()));

    if (min_row <= min_col) {
      // Expand along row j
      int j = row_count.sparsity().row(min_row);

      Matrix<DataType> row = (*this)(j, Slice(0, n));

      std::vector< int > col_i = row.sparsity().getCol();

      for (int k=0; k<row.nnz(); ++k) {
        // Sum up the cofactors
        ret += row.at(k)*cofactor(*this, col_i.at(k), j);
      }
      return ret;
    } else {
      // Expand along col i
      int i = col_count.sparsity().row(min_col);

      Matrix<DataType> col = (*this)(Slice(0, n), i);

      const int* row_i = col.rowPtr();

      for (int k=0; k<col.nnz(); ++k) {
        // Sum up the cofactors
        ret += col.at(k)*cofactor(*this, i, row_i[k]);
      }
      return ret;
    }

  }

  template<typename DataType>
  Matrix<DataType> Matrix<DataType>::zz_sumAll() const {
    // Quick return if empty
    if (isEmpty()) return Matrix<DataType>(1, 1);
    // Sum non-zero elements
    DataType res=0;
    for (int k=0; k<nnz(); k++) {
      res += data()[k];
    }
    return res;
  }

  template<typename DataType>
  Matrix<DataType> Matrix<DataType>::zz_sumCols() const {
    return mul(*this, Matrix<DataType>::ones(size2(), 1));
  }

  template<typename DataType>
  Matrix<DataType> Matrix<DataType>::zz_sumRows() const {
    return mul(Matrix<DataType>::ones(1, size1()), *this);
  }

  template<typename DataType>
  Matrix<DataType> Matrix<DataType>::zz_getMinor(int i, int j) const {
    int n = size2();
    casadi_assert_message(n == size1(), "getMinor: matrix must be square");

    // Trivial return if scalar
    if (n==1) return 1;

    // Remove col i and row j
    Matrix<DataType> M = Matrix<DataType>(n-1, n-1);

    std::vector<int> col = sparsity().getCol();
    const int* row = sparsity().rowPtr();

    for (int k=0;k<nnz();++k) {
      int i1 = col[k];
      int j1 = row[k];

      if (i1 == i || j1 == j) continue;

      int i2 = (i1<i)?i1:i1-1;
      int j2 = (j1<j)?j1:j1-1;

      M(j2, i2) = (*this)(j1, i1);
    }
    return det(M);
  }

  template<typename DataType>
  Matrix<DataType> Matrix<DataType>::zz_cofactor(int i, int j) const {

    // Calculate the i, j minor
    Matrix<DataType> minor_ij = getMinor(*this, i, j);
    // Calculate the cofactor
    int sign_i = 1-2*((i+j) % 2);

    return sign_i * minor_ij;
  }

  template<typename DataType>
  Matrix<DataType> Matrix<DataType>::zz_adj() const {
    int n = size2();
    casadi_assert_message(n == size1(), "adj: matrix must be square");

    // Temporary placeholder
    Matrix<DataType> temp;

    // Cofactor matrix
    Matrix<DataType> C = Matrix<DataType>(n, n);
    for (int i=0; i<n; ++i)
      for (int j=0; j<n; ++j) {
        temp = cofactor(*this, i, j);
        if (!temp.isZero()) C(j, i) = temp;
      }

    return C.T();
  }

  template<typename DataType>
  Matrix<DataType> Matrix<DataType>::zz_inv() const {
    // laplace formula
    return adj(*this)/det(*this);
  }

  template<typename DataType>
  Matrix<DataType> Matrix<DataType>::zz_reshape(int nrow, int ncol) const {
    Sparsity sp = reshape(sparsity(), nrow, ncol);
    return Matrix<DataType>(sp, data(), false);
  }

  template<typename DataType>
  Matrix<DataType> Matrix<DataType>::zz_reshape(const Sparsity& sp) const {
    // quick return if already the right shape
    if (sp==sparsity()) return *this;

    // make sure that the patterns match
    casadi_assert(sp.isReshape(sparsity()));

    return Matrix<DataType>(sp, data(), false);
  }

  template<typename DataType>
  Matrix<DataType> Matrix<DataType>::zz_trace() const {
    casadi_assert_message(size2() == size1(), "trace: must be square");
    DataType res=0;
    for (int i=0; i< size2(); i ++) {
      res += elem(i, i);
    }
    return res;
  }

  template<typename DataType>
  Matrix<DataType> Matrix<DataType>::zz_vecNZ() const {
    return Matrix<DataType>(data());
  }

  template<typename DataType>
  Matrix<DataType>
  Matrix<DataType>::zz_blockcat(const std::vector< std::vector<Matrix<DataType> > > &v) {
    std::vector< Matrix<DataType> > ret;
    for (int i=0; i<v.size(); ++i)
      ret.push_back(horzcat(v[i]));
    return vertcat(ret);
  }

  template<typename DataType>
  Matrix<DataType> Matrix<DataType>::zz_horzcat(const std::vector<Matrix<DataType> > &v) {
    Matrix<DataType> ret;
    for (int i=0; i<v.size(); ++i)
      ret.appendColumns(v[i]);
    return ret;
  }

  template<typename DataType>
  std::vector<Matrix<DataType> >
  Matrix<DataType>::zz_horzsplit(const std::vector<int>& offset) const {
    // Split up the sparsity pattern
    std::vector<Sparsity> sp = horzsplit(sparsity(), offset);

    // Return object
    std::vector<Matrix<DataType> > ret;
    ret.reserve(sp.size());

    // Copy data
    typename std::vector<DataType>::const_iterator data_start=begin(), data_stop;
    for (std::vector<Sparsity>::const_iterator j=sp.begin(); j!=sp.end(); ++j) {
      data_stop = data_start + j->nnz();
      ret.push_back(Matrix<DataType>(*j, std::vector<DataType>(data_start,
                                                               data_stop), false));
      data_start = data_stop;
    }

    // Return the assembled matrix
    casadi_assert(data_stop==end());
    return ret;
  }

  template<typename DataType>
  Matrix<DataType> Matrix<DataType>::zz_vertcat(const std::vector<Matrix<DataType> > &v) {
    Matrix<DataType> ret;
    for (int i=0; i<v.size(); ++i)
      ret.appendColumns(v[i].T());
    return ret.T();
  }

  template<typename DataType>
  std::vector< Matrix<DataType> >
  Matrix<DataType>::zz_vertsplit(const std::vector<int>& offset) const {
    std::vector< Matrix<DataType> > ret = horzsplit(this->T(), offset);
    for (typename std::vector< Matrix<DataType> >::iterator it=ret.begin();
         it!=ret.end(); ++it) {
      *it = it->T();
    }
    return ret;
  }

  template<typename DataType>
  std::vector< Matrix<DataType> >
  Matrix<DataType>::zz_diagsplit(const std::vector<int>& offset1,
                                 const std::vector<int>& offset2) const {
    // Consistency check
    casadi_assert(offset1.size()>=1);
    casadi_assert(offset1.front()==0);
    casadi_assert(offset1.back()==size1());
    casadi_assert(isMonotone(offset1));

    // Consistency check
    casadi_assert(offset2.size()>=1);
    casadi_assert(offset2.front()==0);
    casadi_assert(offset2.back()==size2());
    casadi_assert(isMonotone(offset2));

    // Number of outputs
    int n = offset1.size()-1;

    // Return value
    std::vector< Matrix<DataType> > ret;

    // Caveat: this is a very silly implementation
    for (int i=0; i<n; ++i) {
      ret.push_back((*this)(Slice(offset1[i], offset1[i+1]), Slice(offset2[i], offset2[i+1])));
    }

    return ret;
  }

  template<typename DataType>
  Matrix<DataType> Matrix<DataType>::zz_inner_prod(const Matrix<DataType> &y) const {
    casadi_assert_message(shape()==y.shape(), "inner_prod: Dimension mismatch");
    return sumAll(*this*y);
  }

  template<typename DataType>
  Matrix<DataType> Matrix<DataType>::zz_outer_prod(const Matrix<DataType> &y) const {
    return mul(*this, y.T());
  }

  template<typename DataType>
  Matrix<DataType> Matrix<DataType>::zz_all() const {
    if (!isDense()) return false;
    DataType ret=1;
    for (int i=0;i<nnz();++i) {
      ret = ret && at(i)==1;
    }
    return ret;
  }

  template<typename DataType>
  Matrix<DataType> Matrix<DataType>::zz_any() const {
    if (!isDense()) return false;
    DataType ret=0;
    for (int i=0;i<nnz();++i) {
      ret = ret || at(i)==1;
    }
    return ret;
  }

  template<typename DataType>
  Matrix<DataType> Matrix<DataType>::zz_norm_1() const {
    return sumAll(fabs(*this));
  }

  template<typename DataType>
  Matrix<DataType> Matrix<DataType>::zz_norm_2() const {
    if (isVector()) {
      return norm_F(*this);
    } else {
      casadi_error("2-norms currently only supported for vectors. "
                   "Did you intend to calculate a Frobenius norms (norm_F)?");
    }
  }

  template<typename DataType>
  Matrix<DataType> Matrix<DataType>::zz_norm_F() const {
    return sqrt(sumAll((*this**this)));
  }

  template<typename DataType>
  Matrix<DataType> Matrix<DataType>::zz_norm_inf() const {
    // Get largest element by absolute value
    Matrix<DataType> s = 0;
    for (typename std::vector<DataType>::const_iterator i=begin(); i!=end(); ++i) {
      s = fmax(s, fabs(Matrix<DataType>(*i)));
    }
    return s;
  }

  template<typename DataType>
  void Matrix<DataType>::zz_qr(Matrix<DataType>& Q, Matrix<DataType> &R) const {
    // The following algorithm is taken from J. Demmel:
    // Applied Numerical Linear Algebra (algorithm 3.1.)
    casadi_assert_message(size1()>=size2(), "qr: fewer rows than columns");

    // compute Q and R column by column
    Q = R = Matrix<DataType>();
    for (int i=0; i<size2(); ++i) {
      // Initialize qi to be the i-th column of *this
      Matrix<DataType> ai = (*this)(ALL, i);
      Matrix<DataType> qi = ai;
      // The i-th column of R
      Matrix<DataType> ri = Matrix<DataType>(size2(), 1);

      // subtract the projection of qi in the previous directions from ai
      for (int j=0; j<i; ++j) {

        // Get the j-th column of Q
        Matrix<DataType> qj = Q(ALL, j);

        ri(j, 0) = mul(qi.T(), qj); // Modified Gram-Schmidt
        // ri[j] = inner_prod(qj, ai); // Classical Gram-Schmidt

        // Remove projection in direction j
        if (ri.hasNZ(j, 0))
          qi -= ri(j, 0) * qj;
      }

      // Normalize qi
      ri(i, 0) = norm_2(qi);
      qi /= ri(i, 0);

      // Update R and Q
      Q.appendColumns(qi);
      R.appendColumns(ri);
    }
  }

  template<typename DataType>
  Matrix<DataType> Matrix<DataType>::zz_nullspace() const {
    int n = size1();
    int m = size2();

    Matrix<DataType> X = *this;

    casadi_assert_message(m>=n, "nullspace(): expecting a flat matrix (more columns than rows), "
                          "but got " << dimString() << ".");

    Matrix<DataType> seed = DMatrix::eye(m)(Slice(0, m), Slice(n, m));

    std::vector< Matrix<DataType> > us;
    std::vector< Matrix<DataType> > betas;

    Matrix<DataType> beta;

    for (int i=0;i<n;++i) {
      Matrix<DataType> x = X(i, Slice(i, m));
      Matrix<DataType> u = Matrix<DataType>(x);
      Matrix<DataType> sigma = sqrt(sumCols(x*x));
      const Matrix<DataType>& x0 = x(0, 0);
      u(0, 0) = 1;

      Matrix<DataType> b = -copysign(sigma, x0);

      u(Slice(0), Slice(1, m-i))*= 1/(x0-b);
      beta = 1-x0/b;

      X(Slice(i, n), Slice(i, m)) -=
        beta*mul(mul(X(Slice(i, n), Slice(i, m)), u.T()), u);
      us.push_back(u);
      betas.push_back(beta);
    }

    for (int i=n-1;i>=0;--i) {
      seed(Slice(i, m), Slice(0, m-n)) -=
        betas[i]*mul(us[i].T(), mul(us[i], seed(Slice(i, m), Slice(0, m-n))));
    }

    return seed;

  }

  template<typename DataType>
  Matrix<DataType> Matrix<DataType>::zz_solve(const Matrix<DataType>& b) const {
    // check dimensions
    casadi_assert_message(size1() == b.size1(), "solve Ax=b: dimension mismatch: b has "
                          << b.size1() << " rows while A has " << size1() << ".");
    casadi_assert_message(size1() == size2(), "solve: A not square but " << dimString());

    if (isTril()) {
      // forward substitution if lower triangular
      Matrix<DataType> x = b;
      const int*  Arow = rowPtr();
      const int*  Acolind = colindPtr();
      const std::vector<DataType> & Adata = data();
      for (int i=0; i<size2(); ++i) { // loop over columns forwards
        for (int k=0; k<b.size2(); ++k) { // for every right hand side
          if (!x.hasNZ(i, k)) continue;
          x(i, k) /= (*this)(i, i);
          for (int kk=Acolind[i+1]-1; kk>=Acolind[i] && Arow[kk]>i; --kk) {
            int j = Arow[kk];
            x(j, k) -= Adata[kk]*x(i, k);
          }
        }
      }
      return x;
    } else if (isTriu()) {
      // backward substitution if upper triangular
      Matrix<DataType> x = b;
      const int*  Arow = rowPtr();
      const int*  Acolind = colindPtr();
      const std::vector<DataType> & Adata = data();
      for (int i=size2()-1; i>=0; --i) { // loop over columns backwards
        for (int k=0; k<b.size2(); ++k) { // for every right hand side
          if (!x.hasNZ(i, k)) continue;
          x(i, k) /= (*this)(i, i);
          for (int kk=Acolind[i]; kk<Acolind[i+1] && Arow[kk]<i; ++kk) {
            int j = Arow[kk];
            x(j, k) -= Adata[kk]*x(i, k);
          }
        }
      }
      return x;
    } else if (hasNonStructuralZeros()) {

      // If there are structurally nonzero entries that are known to be zero,
      // remove these and rerun the algorithm
      Matrix<DataType> A_sparse = *this;
      A_sparse.makeSparse();
      return solve(A_sparse, b);

    } else {

      // Make a BLT transformation of A
      std::vector<int> rowperm, colperm, rowblock, colblock, coarse_rowblock, coarse_colblock;
      sparsity().dulmageMendelsohn(rowperm, colperm, rowblock, colblock,
                                     coarse_rowblock, coarse_colblock);

      // Permute the right hand side
      Matrix<DataType> bperm = b(rowperm, ALL);

      // Permute the linear system
      Matrix<DataType> Aperm = (*this)(rowperm, colperm);

      // Solution
      Matrix<DataType> xperm;

      // Solve permuted system
      if (Aperm.isTril()) {

        // Forward substitution if lower triangular
        xperm = solve(Aperm, bperm);

      } else if (size2()<=3) {

        // Form inverse by minor expansion and multiply if very small (up to 3-by-3)
        xperm = mul(inv(Aperm), bperm);

      } else {

        // Make a QR factorization
        Matrix<DataType> Q, R;
        qr(Aperm, Q, R);

        // Solve the factorized system (note that solve will now be fast since it is triangular)
        xperm = solve(R, mul(Q.T(), bperm));
      }

      // get the inverted column permutation
      std::vector<int> inv_colperm(colperm.size());
      for (int k=0; k<colperm.size(); ++k)
        inv_colperm[colperm[k]] = k;

      // Permute back the solution and return
      Matrix<DataType> x = xperm(inv_colperm, ALL);
      return x;
    }
  }

  template<typename DataType>
  Matrix<DataType> Matrix<DataType>::zz_pinv() const {
    if (size2()>=size1()) {
      return solve(mul(*this, T()), *this).T();
    } else {
      return solve(mul(this->T(), *this), this->T());
    }
  }

  template<typename DataType>
  Matrix<DataType> Matrix<DataType>::zz_kron(const Matrix<DataType>& b) const {
    const Sparsity &a_sp = sparsity();
    Matrix<DataType> filler = Matrix<DataType>(b.shape());
    std::vector< std::vector< Matrix<DataType> > >
      blocks(size1(), std::vector< Matrix<DataType> >(size2(), filler));
    for (int i=0;i<size1();++i) {
      for (int j=0;j<size2();++j) {
        int k = a_sp.getNZ(i, j);
        if (k!=-1) {
          blocks[i][j] = (*this)[k]*b;
        }
      }
    }
    return blockcat(blocks);
  }

  template<typename DataType>
  Matrix<DataType> Matrix<DataType>::zz_repmat(const Sparsity& sp) const {
    casadi_assert_message(isScalar(),
                          "repmat(Matrix<DataType> x, Sparsity sp) only defined for scalar x");
    return Matrix<DataType>(sp, toScalar(), false);
  }

  template<typename DataType>
  Matrix<DataType> Matrix<DataType>::zz_repmat(int n, int m) const {
    if (n==1 &&  m==1) {
      // Quick return if possible
      return *this;
    } else if (isScalar()) {
      if (isDense()) {
        return Matrix<DataType>(Sparsity::dense(n, m), toScalar(), false);
      } else {
        return Matrix<DataType>(n, m);
      }
    } else {
      std::vector< Matrix<DataType> > v_hor(m, *this);
      std::vector< Matrix<DataType> > v_ver(n, horzcat(v_hor));
      return vertcat(v_ver);
    }
  }

  template<typename DataType>
  Matrix<DataType> Matrix<DataType>::zz_diag() const {
    // Nonzero mapping
    std::vector<int> mapping;
    // Get the sparsity
    Sparsity sp = sparsity().getDiag(mapping);

    Matrix<DataType> ret = Matrix<DataType>(sp);

    for (int k=0;k<mapping.size();k++) ret[k] = (*this)[mapping[k]];
    return ret;
  }

  /** \brief   Construct a matrix with given block on the diagonal */
  template<typename DataType>
  Matrix<DataType> Matrix<DataType>::zz_diagcat(const std::vector< Matrix<DataType> > &A) {
    std::vector<DataType> data;

    std::vector<Sparsity> sp;
    for (int i=0;i<A.size();++i) {
      data.insert(data.end(), A[i].data().begin(), A[i].data().end());
      sp.push_back(A[i].sparsity());
    }

    return Matrix<DataType>(diagcat(sp), data, false);
  }

  template<typename DataType>
  Matrix<DataType> Matrix<DataType>::zz_unite(const Matrix<DataType>& B) const {
    // Join the sparsity patterns
    std::vector<unsigned char> mapping;
    Sparsity sp = sparsity().patternUnion(B.sparsity(), mapping);

    // Create return matrix
    Matrix<DataType> ret(sp);

    // Copy sparsity
    int elA=0, elB=0;
    for (int k=0; k<mapping.size(); ++k) {
      if (mapping[k]==1) {
        ret.data()[k] = data()[elA++];
      } else if (mapping[k]==2) {
        ret.data()[k] = B.data()[elB++];
      } else {
        throw CasadiException("Pattern intersection not empty");
      }
    }

    casadi_assert(nnz()==elA);
    casadi_assert(B.nnz()==elB);

    return ret;
  }

  template<typename DataType>
  Matrix<DataType> Matrix<DataType>::zz_polyval(const Matrix<DataType>& x) const {
    casadi_assert_message(isDense(), "polynomial coefficients vector must be dense");
    casadi_assert_message(isVector() && nnz()>0, "polynomial coefficients must be a vector");
    Matrix<DataType> ret = (*this)[0];
    for (int i=1; i<nnz(); ++i) {
      ret = ret*x + (*this)[i];
    }
    return ret;
  }

  template<typename DataType>
  void Matrix<DataType>::zz_addMultiple(const std::vector<DataType>& v,
                                     std::vector<DataType>& res, bool trans_A) const {
    // Get dimension and sparsity
    int d1=size2(), d2=size1();
    const int* colind=this->colindPtr();
    const int* row=this->rowPtr();
    const std::vector<DataType>& data = this->data();

    // Assert consistent dimensions
    if (trans_A) {
      casadi_assert(v.size()==d1);
      casadi_assert(res.size()==d2);
    } else {
      casadi_assert(v.size()==d2);
      casadi_assert(res.size()==d1);
    }

    // Carry out multiplication
    for (int i=0; i<d1; ++i) { // loop over cols
      for (int el=colind[i]; el<colind[i+1]; ++el) { // loop over the non-zero elements
        int j=row[el];  // row
        // Add scalar product
        if (trans_A) {
          res[j] += v[i]*data[el];
        } else {
          res[i] += v[j]*data[el];
        }
      }
    }
  }

  template<typename DataType>
  Matrix<DataType> Matrix<DataType>::zz_project(const Sparsity& sp) const {
    // Check dimensions
    if (!(isEmpty() && sp.numel()==0)) {
      casadi_assert_message(size2()==sp.size2() && size1()==sp.size1(),
                            "Shape mismatch. Expecting " << dimString() << ", but got " <<
                            sp.dimString() << " instead.");
    }

    // Return value
    Matrix<DataType> ret(sp, 0, false);

    // Get the elements of the known matrix
    std::vector<int> known_ind = sparsity().find();

    // Find the corresponding nonzeros in the return matrix
    sp.getNZ(known_ind);

    // Set the element values
    const std::vector<DataType>& A_data = data();
    std::vector<DataType>& ret_data = ret.data();
    for (int k=0; k<known_ind.size(); ++k) {
      if (known_ind[k]!=-1) {
        ret_data[known_ind[k]] = A_data[k];
      }
    }
    return ret;
  }

  template<typename DataType>
  int Matrix<DataType>::zz_sprank() const {
    return sprank(sparsity());
  }

  template<typename DataType>
  int Matrix<DataType>::zz_norm_0_mul_nn(const Matrix<DataType> &A,
                                         std::vector<bool>& Bwork,
                                         std::vector<int>& Iwork) const {

    // Note: because the algorithm works with compressed row storage,
    // we have x=B and y=A

    casadi_assert_message(A.size1()==size2(), "Dimension error. Got " << dimString()
                          << " times " << A.dimString() << ".");

    int n_row = A.size2();
    int n_col = size1();

    casadi_assert_message(Bwork.size()>=n_col,
      "We need a bigger work vector (>=" << n_col << "), but got "<< Bwork.size() <<".");
    casadi_assert_message(Iwork.size()>=n_row+1+n_col,
      "We need a bigger work vector (>=" << n_row+1+n_col << "), but got "<< Iwork.size() <<".");

    const int* Aj = A.rowPtr();
    const int* Ap = A.colindPtr();

    const int* Bj = rowPtr();
    const int* Bp = colindPtr();

    int *Cp = getPtr(Iwork);
    int *mask = Cp+n_row+1;

    // Implementation borrowed from Scipy's sparsetools/csr.h

    // Pass 1

    // method that uses O(n) temp storage
    std::fill(mask, mask+n_col, -1);

    Cp[0] = 0;
    int nnz = 0;

    for (int i = 0; i < n_row; i++) {
      int row_nnz = 0;
      for (int jj = Ap[i]; jj < Ap[i+1]; jj++) {
        int j = Aj[jj];
        for (int kk = Bp[j]; kk < Bp[j+1]; kk++) {
          int k = Bj[kk];
          if (mask[k] != i) {
            mask[k] = i;
            row_nnz++;
          }
        }
      }
      int next_nnz = nnz + row_nnz;

      nnz = next_nnz;
      Cp[i+1] = nnz;
    }

    // Pass 2
    int *next = getPtr(Iwork) + n_row+1;
    std::fill(next, next+n_col, -1);

    std::vector<bool> & sums = Bwork;
    std::fill(sums.begin(), sums.end(), false);

    nnz = 0;

    Cp[0] = 0;

    for (int i = 0; i < n_row; i++) {
        int head   = -2;
        int length =  0;

        int jj_start = Ap[i];
        int jj_end   = Ap[i+1];
        for (int jj = jj_start; jj < jj_end; jj++) {
            int j = Aj[jj];

            int kk_start = Bp[j];
            int kk_end   = Bp[j+1];
            for (int kk = kk_start; kk < kk_end; kk++) {
                int k = Bj[kk];

                sums[k] = true;

                if (next[k] == -1) {
                    next[k] = head;
                    head  = k;
                    length++;
                }
            }
        }

        for (int jj = 0; jj < length; jj++) {

            if (sums[head]) {
                nnz++;
            }

            int temp = head;
            head = next[head];

            next[temp] = -1; //clear arrays
            sums[temp] =  0;
        }

        Cp[i+1] = nnz;
    }

    return nnz;

  }

  template<typename DataType>
  DataType Matrix<DataType>::zz_norm_inf_mul_nn(const Matrix<DataType> &A,
                                                std::vector<DataType>& Dwork,
                                                std::vector<int>& Iwork) const {
    const Matrix<DataType> &B = *this;
    // Note: because the algorithm works with compressed row storage,
    // we have x=B and y=A
    DataType res = 0;

    casadi_assert_message(A.size1()==B.size2(), "Dimension error. Got " << B.dimString()
                          << " times " << A.dimString() << ".");


    int n_row = A.size2();
    int n_col = B.size1();

    casadi_assert_message(Dwork.size()>=n_col,
                          "We need a bigger work vector (>="
                          << n_col << "), but got " << Dwork.size() <<".");
    casadi_assert_message(Iwork.size()>=n_row+1+n_col,
                          "We need a bigger work vector (>=" << n_row+1+n_col
                          << "), but got " << Iwork.size() <<".");

    const int* Aj = A.rowPtr();
    const int* Ap = A.colindPtr();
    const std::vector<DataType> &Ax = A.data();

    const int* Bj = B.rowPtr();
    const int* Bp = B.colindPtr();
    const std::vector<DataType> &Bx = B.data();

    int *Cp = getPtr(Iwork);
    int *mask = Cp + n_row+1;

    // Implementation borrowed from Scipy's sparsetools/csr.h

    // Pass 1

    // method that uses O(n) temp storage
    std::fill(mask, mask+n_col, -1);

    Cp[0] = 0;
    int nnz = 0;

    for (int i = 0; i < n_row; i++) {
      int row_nnz = 0;
      for (int jj = Ap[i]; jj < Ap[i+1]; jj++) {
        int j = Aj[jj];
        for (int kk = Bp[j]; kk < Bp[j+1]; kk++) {
          int k = Bj[kk];
          if (mask[k] != i) {
            mask[k] = i;
            row_nnz++;
          }
        }
      }
      int next_nnz = nnz + row_nnz;

      nnz = next_nnz;
      Cp[i+1] = nnz;
    }

    // Pass 2
    int *next = &Iwork[n_row+1];
    std::fill(next, next+n_col, -1);

    DataType* sums = &Dwork[0];
    std::fill(sums, sums+n_col, 0);

    nnz = 0;

    Cp[0] = 0;

    for (int i = 0; i < n_row; i++) {
      int head   = -2;
      int length =  0;
      int jj_start = Ap[i];
      int jj_end   = Ap[i+1];
      for (int jj = jj_start; jj < jj_end; jj++) {
        int j = Aj[jj];
        DataType v = Ax[jj];
        int kk_start = Bp[j];
        int kk_end   = Bp[j+1];
        for (int kk = kk_start; kk < kk_end; kk++) {
          int k = Bj[kk];
          sums[k] += v*Bx[kk];
          if (next[k] == -1) {
            next[k] = head;
            head  = k;
            length++;
          }
        }
      }

      for (int jj = 0; jj < length; jj++) {
        if (!casadi_limits<DataType>::isZero(sums[head])) {
          res = fmax(res, abs(sums[head]));
          nnz++;
        }
        int temp = head;
        head = next[head];
        next[temp] = -1; //clear arrays
        sums[temp] =  0;
      }

      Cp[i+1] = nnz;
    }

    return res;
  }

  template<typename DataType>
  void Matrix<DataType>::zz_expand(Matrix<DataType> &weights, Matrix<DataType>& terms) const {
    throw CasadiException("\"expand\" not defined for instantiation");
  }

  template<typename DataType>
  Matrix<DataType> Matrix<DataType>::zz_pw_const(const Matrix<DataType> &tval,
                                                 const Matrix<DataType> &val) const {
    throw CasadiException("\"pw_const\" not defined for instantiation");
    return Matrix<DataType>();
  }

  template<typename DataType>
  Matrix<DataType> Matrix<DataType>::zz_pw_lin(const Matrix<DataType> &tval,
                                               const Matrix<DataType> &val) const {
    throw CasadiException("\"pw_lin\" not defined for instantiation");
    return Matrix<DataType>();
  }

  template<typename DataType>
  Matrix<DataType> Matrix<DataType>::zz_if_else(const Matrix<DataType> &if_true,
                                                const Matrix<DataType> &if_false) const {
    throw CasadiException("\"if_else\" not defined for instantiation");
    return Matrix<DataType>();
  }

  template<typename DataType>
  Matrix<DataType> Matrix<DataType>::zz_heaviside() const {
    throw CasadiException("\"heaviside\" not defined for instantiation");
    return Matrix<DataType>();
  }

  template<typename DataType>
  Matrix<DataType> Matrix<DataType>::zz_rectangle() const {
    throw CasadiException("\"rectangle\" not defined for instantiation");
    return Matrix<DataType>();
  }

  template<typename DataType>
  Matrix<DataType> Matrix<DataType>::zz_triangle() const {
    throw CasadiException("\"triangle\" not defined for instantiation");
    return Matrix<DataType>();
  }

  template<typename DataType>
  Matrix<DataType> Matrix<DataType>::zz_ramp() const {
    throw CasadiException("\"ramp\" not defined for instantiation");
    return Matrix<DataType>();
  }

  template<typename DataType>
  Matrix<DataType> Matrix<DataType>::zz_gauss_quadrature(const Matrix<DataType> &x,
                                                         const Matrix<DataType> &a,
                                                         const Matrix<DataType> &b, int order,
                                                         const Matrix<DataType>& w) const {
    throw CasadiException("\"gauss_quadrature\" not defined for instantiation");
    return Matrix<DataType>();
  }

  template<typename DataType>
  Matrix<DataType> Matrix<DataType>::zz_simplify() const {
    throw CasadiException("\"simplify\" not defined for instantiation");
    return Matrix<DataType>();
  }

  template<typename DataType>
  Matrix<DataType> Matrix<DataType>::zz_substitute(const Matrix<DataType>& v,
                                                   const Matrix<DataType>& vdef) const {
    throw CasadiException("\"substitute\" not defined for instantiation");
    return Matrix<DataType>();
  }

  template<typename DataType>
  std::vector<Matrix<DataType> >
  Matrix<DataType>::zz_substitute(const std::vector<Matrix<DataType> >& ex,
                                  const std::vector<Matrix<DataType> >& v,
                                  const std::vector<Matrix<DataType> >& vdef) {
    throw CasadiException("\"substitute\" not defined for instantiation");
    return std::vector<Matrix<DataType> >();
  }

  template<typename DataType>
  void Matrix<DataType>::zz_substituteInPlace(const std::vector<Matrix<DataType> >& v,
                                              std::vector<Matrix<DataType> >& vdef,
                                              std::vector<Matrix<DataType> >& ex,
                                              bool reverse) {
    throw CasadiException("\"substituteInPlace\" not defined for instantiation");
  }

  template<typename DataType>
  Matrix<DataType> Matrix<DataType>::zz_spy() const {
    throw CasadiException("\"spy\" not defined for instantiation");
    return Matrix<DataType>();
  }

  template<typename DataType>
  bool Matrix<DataType>::zz_dependsOn(const Matrix<DataType> &arg) const {
    throw CasadiException("\"dependsOn\" not defined for instantiation");
    return false;
  }

  template<typename DataType>
  Matrix<DataType> Matrix<DataType>::zz_jacobian(const Matrix<DataType> &arg) const {
    throw CasadiException("\"jacobian\" not defined for instantiation");
    return Matrix<DataType>();
  }

  template<typename DataType>
  Matrix<DataType> Matrix<DataType>::zz_gradient(const Matrix<DataType> &arg) const {
    throw CasadiException("\"gradient\" not defined for instantiation");
    return Matrix<DataType>();
  }

  template<typename DataType>
  Matrix<DataType> Matrix<DataType>::zz_tangent(const Matrix<DataType> &arg) const {
    throw CasadiException("\"tangent\" not defined for instantiation");
    return Matrix<DataType>();
  }

  template<typename DataType>
  Matrix<DataType> Matrix<DataType>::zz_hessian(const Matrix<DataType> &arg) const {
    throw CasadiException("\"hessian\" not defined for instantiation");
    return Matrix<DataType>();
  }

  template<typename DataType>
  void Matrix<DataType>::zz_hessian(const Matrix<DataType> &arg, Matrix<DataType> &H,
                                    Matrix<DataType> &g) const {
    throw CasadiException("\"hessian\" not defined for instantiation");
  }

  template<typename DataType>
  Matrix<DataType> Matrix<DataType>::zz_jacobianTimesVector(const Matrix<DataType> &arg,
                                                            const Matrix<DataType> &v,
                                                            bool transpose_jacobian) const {
    throw CasadiException("\"jacobianTimesVector\" not defined for instantiation");
    return Matrix<DataType>();
  }

  template<typename DataType>
  Matrix<DataType> Matrix<DataType>::zz_taylor(const Matrix<DataType>& x,
                                               const Matrix<DataType>& a, int order) const {
    throw CasadiException("\"taylor\" not defined for instantiation");
    return Matrix<DataType>();
  }

  template<typename DataType>
  Matrix<DataType> Matrix<DataType>::zz_mtaylor(const Matrix<DataType>& x,
                                                const Matrix<DataType>& a, int order) const {
    throw CasadiException("\"mtaylor\" not defined for instantiation");
    return Matrix<DataType>();
  }

  template<typename DataType>
  Matrix<DataType> Matrix<DataType>::zz_mtaylor(const Matrix<DataType>& x,
                                                const Matrix<DataType>& a, int order,
                                                const std::vector<int>&order_contributions) const {
    throw CasadiException("\"mtaylor\" not defined for instantiation");
    return Matrix<DataType>();
  }

  template<typename DataType>
  int Matrix<DataType>::zz_countNodes() const {
    throw CasadiException("\"countNodes\" not defined for instantiation");
    return 0;
  }

  template<typename DataType>
  std::string
  Matrix<DataType>::zz_getOperatorRepresentation(const std::vector<std::string>& args) const {
    throw CasadiException("\"getOperatorRepresentation\" not defined for instantiation");
    return std::string();
  }

  template<typename DataType>
  std::vector<Matrix<DataType> > Matrix<DataType>::zz_getSymbols() const {
    throw CasadiException("\"getSymbols\" not defined for instantiation");
    return std::vector<Matrix<DataType> >();
  }

  template<typename DataType>
  std::vector<Matrix<DataType> >
  Matrix<DataType>::zz_getSymbols(const std::vector<Matrix<DataType> >& e) {
    throw CasadiException("\"getSymbols\" not defined for instantiation");
    return std::vector<Matrix<DataType> >();
  }

  template<typename DataType>
  void Matrix<DataType>::zz_extractShared(std::vector<Matrix<DataType> >& ex,
                                          std::vector<Matrix<DataType> >& v,
                                          std::vector<Matrix<DataType> >& vdef,
                                          const std::string& v_prefix,
                                          const std::string& v_suffix) {
    throw CasadiException("\"extractShared\" not defined for instantiation");
  }

  template<typename DataType>
  void Matrix<DataType>::zz_printCompact(std::ostream &stream) const {
    throw CasadiException("\"printCompact\" not defined for instantiation");
  }

  template<typename DataType>
  Matrix<DataType> Matrix<DataType>::zz_poly_coeff(const Matrix<DataType>&x) const {
    throw CasadiException("\"poly_coeff\" not defined for instantiation");
  }

  template<typename DataType>
  Matrix<DataType> Matrix<DataType>::zz_poly_roots() const {
    throw CasadiException("\"poly_roots\" not defined for instantiation");
  }

  template<typename DataType>
  Matrix<DataType> Matrix<DataType>::zz_eig_symbolic() const {
    throw CasadiException("\"eig_symbolic\" not defined for instantiation");
  }

  template<typename DataType>
  Matrix<DataType> Matrix<DataType>::zz_sparsify(double tol) const {
    Matrix<DataType> ret = *this;
    ret.makeSparse(tol);
    return ret;
  }

} // namespace casadi

/// \endcond

#endif // CASADI_MATRIX_IMPL_HPP
