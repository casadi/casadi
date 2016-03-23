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

#include "function/function.hpp"

#include "casadi_interrupt.hpp"

/// \cond INTERNAL

namespace casadi {
  // Implementations

  template<typename Scalar>
  int Matrix<Scalar>::stream_precision_ = 6;
  template<typename Scalar>
  int Matrix<Scalar>::stream_width_ = 0;
  template<typename Scalar>
  bool Matrix<Scalar>::stream_scientific_ = false;

  template<typename Scalar>
  bool Matrix<Scalar>::__nonzero__() const {
    if (numel()!=1) {casadi_error("Only scalar Matrix could have a truth value, but you "
                                  "provided a shape" << dim());}
    return nonzeros().at(0)!=0;
  }

  template<typename Scalar>
  void Matrix<Scalar>::get(Matrix<Scalar>& m, bool ind1,
                                const Slice& rr, const Slice& cc) const {
    // Both are scalar
    if (rr.is_scalar(size1()) && cc.is_scalar(size2())) {
      int k = sparsity().get_nz(rr.scalar(size1()), cc.scalar(size2()));
      if (k>=0) {
        m = nonzeros().at(k);
      } else {
        m = Matrix<Scalar>(1, 1);
      }
      return;
    }

    // Fall back on IM-IM
    get(m, ind1, rr.all(size1(), ind1), cc.all(size2(), ind1));
  }

  template<typename Scalar>
  void Matrix<Scalar>::get(Matrix<Scalar>& m, bool ind1,
                                const Slice& rr, const Matrix<int>& cc) const {
    // Fall back on IM-IM
    get(m, ind1, rr.all(size1(), ind1), cc);
  }

  template<typename Scalar>
  void Matrix<Scalar>::get(Matrix<Scalar>& m, bool ind1,
                                const Matrix<int>& rr, const Slice& cc) const {
    // Fall back on IM-IM
    get(m, ind1, rr, cc.all(size2(), ind1));
  }

  template<typename Scalar>
  void Matrix<Scalar>::get(Matrix<Scalar>& m, bool ind1,
                             const Matrix<int>& rr, const Matrix<int>& cc) const {
    // Scalar
    if (rr.is_scalar(true) && cc.is_scalar(true)) {
      return get(m, ind1, to_slice(rr, ind1), to_slice(cc, ind1));
    }

    // Make sure dense vectors
    casadi_assert_message(rr.is_dense() && rr.is_vector(),
                          "Marix::get: First index must be a dense vector");
    casadi_assert_message(cc.is_dense() && cc.is_vector(),
                          "Marix::get: Second index must be a dense vector");

    // Get the sparsity pattern - does bounds checking
    std::vector<int> mapping;
    Sparsity sp = sparsity().sub(rr.nonzeros(), cc.nonzeros(), mapping, ind1);

    // Copy nonzeros
    m = Matrix<Scalar>::zeros(sp);
    for (int k=0; k<mapping.size(); ++k) m->at(k) = nonzeros().at(mapping[k]);
  }

  template<typename Scalar>
  void Matrix<Scalar>::get(Matrix<Scalar>& m, bool ind1, const Slice& rr) const {
    // Scalar
    if (rr.is_scalar(numel())) {
      int r = rr.scalar(numel());
      int k = sparsity().get_nz(r % size1(), r / size1());
      if (k>=0) {
        m = nonzeros().at(k);
      } else {
        m = Matrix<Scalar>(1, 1);
      }
      return;
    }

    // Fall back on IM
    get(m, ind1, rr.all(numel(), ind1));
  }

  template<typename Scalar>
  void Matrix<Scalar>::get(Matrix<Scalar>& m, bool ind1, const Matrix<int>& rr) const {
    // Scalar
    if (rr.is_scalar(true)) {
      return get(m, ind1, to_slice(rr, ind1));
    }

    // If the indexed matrix is dense, use nonzero indexing
    if (is_dense()) {
      return get_nz(m, ind1, rr);
    }

    // Get the sparsity pattern - does bounds checking
    std::vector<int> mapping;
    Sparsity sp = sparsity().sub(rr.nonzeros(), rr.sparsity(), mapping, ind1);

    // If indexed matrix was a row/column vector, make sure that the result is too
    bool tr = (is_column() && rr.is_row()) || (is_row() && rr.is_column());

    // Copy nonzeros
    m = Matrix<Scalar>::zeros(tr ? sp.T() : sp);
    for (int k=0; k<mapping.size(); ++k) m->at(k) = nonzeros().at(mapping[k]);
  }

  template<typename Scalar>
  void Matrix<Scalar>::get(Matrix<Scalar>& m, bool ind1, const Sparsity& sp) const {
    casadi_assert_message(size()==sp.size(),
                          "get(Sparsity sp): shape mismatch. This matrix has shape "
                          << size() << ", but supplied sparsity index has shape "
                          << sp.size() << ".");
    m = project(*this, sp);
  }

  template<typename Scalar>
  void Matrix<Scalar>::set(const Matrix<Scalar>& m, bool ind1,
                             const Slice& rr, const Slice& cc) {
    // Both are scalar
    if (rr.is_scalar(size1()) && cc.is_scalar(size2()) && m.is_dense()) {
      int oldsize = sparsity_.nnz();
      int ind = sparsity_.add_nz(rr.scalar(size1()), cc.scalar(size2()));
      if (oldsize == sparsity_.nnz()) {
        nonzeros_.at(ind) = m.scalar();
      } else {
        nonzeros_.insert(nonzeros_.begin()+ind, m.scalar());
      }
      return;
    }

    // Fall back on (IM, IM)
    set(m, ind1, rr.all(size1(), ind1), cc.all(size2(), ind1));
  }

  template<typename Scalar>
  void Matrix<Scalar>::set(const Matrix<Scalar>& m, bool ind1,
                             const Slice& rr, const Matrix<int>& cc) {
    // Fall back on (IM, IM)
    set(m, ind1, rr.all(size1(), ind1), cc);
  }

  template<typename Scalar>
  void Matrix<Scalar>::set(const Matrix<Scalar>& m, bool ind1,
                                const Matrix<int>& rr, const Slice& cc) {
    // Fall back on (IM, IM)
    set(m, ind1, rr, cc.all(size2(), ind1));
  }

  template<typename Scalar>
  void Matrix<Scalar>::set(const Matrix<Scalar>& m, bool ind1,
                                const Matrix<int>& rr, const Matrix<int>& cc) {
    // Scalar
    if (rr.is_scalar(true) && cc.is_scalar(true) && m.is_dense()) {
      return set(m, ind1, to_slice(rr, ind1), to_slice(cc, ind1));
    }

    // Row vector rr (e.g. in MATLAB) is transposed to column vector
    if (rr.size1()==1 && rr.size2()>1) {
      return set(m, ind1, rr.T(), cc);
    }

    // Row vector cc (e.g. in MATLAB) is transposed to column vector
    if (cc.size1()==1 && cc.size2()>1) {
      return set(m, ind1, rr, cc.T());
    }

    // Make sure rr and cc are dense vectors
    casadi_assert_message(rr.is_dense() && rr.is_column(),
                          "Matrix::set: First index not dense vector");
    casadi_assert_message(cc.is_dense() && cc.is_column(),
                          "Matrix::set: Second index not dense vector");

    // Assert dimensions of assigning matrix
    if (rr.size1() != m.size1() || cc.size1() != m.size2()) {
      if (m.is_scalar()) {
        // m scalar means "set all"
        return set(repmat(m, rr.size1(), cc.size1()), ind1, rr, cc);
      } else if (rr.size1() == m.size2() && cc.size1() == m.size1()
                 && std::min(m.size1(), m.size2()) == 1) {
        // m is transposed if necessary
        return set(m.T(), ind1, rr, cc);
      } else {
        // Error otherwise
        casadi_error("Dimension mismatch." << "lhs is " << rr.size1() << "-by-"
                     << cc.size1() << ", while rhs is " << m.size());
      }
    }

    // Dimensions
    int sz1 = size1(), sz2 = size2();

    // Report out-of-bounds
    if (!inBounds(rr.nonzeros(), -sz1+ind1, sz1+ind1)) {
      casadi_error("set[., r, c] out of bounds. Your rr contains "
                   << *std::min_element(rr->begin(), rr->end()) << " up to "
                   << *std::max_element(rr->begin(), rr->end())
                   << ", which is outside the range [" << -sz1+ind1 << ","<< sz1+ind1 <<  ").");
    }
    if (!inBounds(cc.nonzeros(), -sz2+ind1, sz2+ind1)) {
      casadi_error("set [., r, c] out of bounds. Your cc contains "
                   << *std::min_element(cc->begin(), cc->end()) << " up to "
                   << *std::max_element(cc->begin(), cc->end())
                   << ", which is outside the range [" << -sz2+ind1 << ","<< sz2+ind1 <<  ").");
    }

    // If we are assigning with something sparse, first remove existing entries
    if (!m.is_dense()) {
      erase(rr.nonzeros(), cc.nonzeros(), ind1);
    }

    // Collect all assignments
    IM el = IM::zeros(m.sparsity());
    for (int j=0; j<el.size2(); ++j) { // Loop over columns of m
      int this_j = cc->at(j) - ind1; // Corresponding column in this
      if (this_j<0) this_j += sz2;
      for (int k=el.colind(j); k<el.colind(j+1); ++k) { // Loop over rows of m
        int i = m.row(k);
        int this_i = rr->at(i) - ind1; // Corresponding row in this
        if (this_i<0) this_i += sz1;
        el->at(k) = this_i + this_j*sz1;
      }
    }
    return set(m, false, el);
  }

  template<typename Scalar>
  void Matrix<Scalar>::set(const Matrix<Scalar>& m, bool ind1, const Slice& rr) {
    // Scalar
    if (rr.is_scalar(numel()) && m.is_dense()) {
      int r = rr.scalar(numel());
      int oldsize = sparsity_.nnz();
      int ind = sparsity_.add_nz(r % size1(), r / size1());
      if (oldsize == sparsity_.nnz()) {
        nonzeros_.at(ind) = m.scalar();
      } else {
        nonzeros_.insert(nonzeros_.begin()+ind, m.scalar());
      }
      return;
    }

    // Fall back on IM
    set(m, ind1, rr.all(numel(), ind1));
  }

  template<typename Scalar>
  void Matrix<Scalar>::set(const Matrix<Scalar>& m, bool ind1, const Matrix<int>& rr) {
    // Scalar
    if (rr.is_scalar(true) && m.is_dense()) {
      return set(m, ind1, to_slice(rr, ind1));
    }

    // Assert dimensions of assigning matrix
    if (rr.sparsity() != m.sparsity()) {
      if (rr.size() == m.size()) {
        // Remove submatrix to be replaced
        erase(rr.nonzeros(), ind1);

        // Find the intersection between rr's and m's sparsity patterns
        Sparsity sp = rr.sparsity() * m.sparsity();

        // Project both matrices to this sparsity
        return set(project(m, sp), ind1, Matrix<int>::project(rr, sp));
      } else if (m.is_scalar()) {
        // m scalar means "set all"
        if (m.is_dense()) {
          return set(Matrix<Scalar>(rr.sparsity(), m), ind1, rr);
        } else {
          return set(Matrix<Scalar>(rr.size()), ind1, rr);
        }
      } else if (rr.size1() == m.size2() && rr.size2() == m.size1()
                 && std::min(m.size1(), m.size2()) == 1) {
        // m is transposed if necessary
        return set(m.T(), ind1, rr);
      } else {
        // Error otherwise
        casadi_error("Dimension mismatch." << "lhs is " << rr.size()
                     << ", while rhs is " << m.size());
      }
    }

    // Dimensions of this
    int sz1 = size1(), sz2 = size2(), sz = nnz(), nel = numel(), rrsz = rr.nnz();

    // Quick return if nothing to set
    if (rrsz==0) return;

    // Check bounds
    if (!inBounds(rr.nonzeros(), -nel+ind1, nel+ind1)) {
      casadi_error("set[rr] out of bounds. Your rr contains "
                   << *std::min_element(rr->begin(), rr->end()) << " up to "
                   << *std::max_element(rr->begin(), rr->end())
                   << ", which is outside the range [" << -nel+ind1 << ","<< nel+ind1 <<  ").");
    }

    // Dense mode
    if (is_dense() && m.is_dense()) {
      return set_nz(m, ind1, rr);
    }

    // Construct new sparsity pattern
    std::vector<int> new_row=sparsity().get_row(), new_col=sparsity().get_col(), nz(rr.nonzeros());
    new_row.reserve(sz+rrsz);
    new_col.reserve(sz+rrsz);
    nz.reserve(rrsz);
    for (std::vector<int>::iterator i=nz.begin(); i!=nz.end(); ++i) {
      if (ind1) (*i)--;
      if (*i<0) *i += nel;
      new_row.push_back(*i % sz1);
      new_col.push_back(*i / sz1);
    }
    Sparsity sp = Sparsity::triplet(sz1, sz2, new_row, new_col);

    // If needed, update pattern
    if (sp != sparsity()) *this = project(*this, sp);

    // Find the nonzeros corresponding to rr
    sparsity().get_nz(nz);

    // Carry out the assignments
    for (int i=0; i<nz.size(); ++i) {
      nonzeros().at(nz[i]) = m->at(i);
    }
  }

  template<typename Scalar>
  void Matrix<Scalar>::set(const Matrix<Scalar>& m, bool ind1, const Sparsity& sp) {
    casadi_assert_message(size()==sp.size(),
                          "set(Sparsity sp): shape mismatch. This matrix has shape "
                          << size() << ", but supplied sparsity index has shape "
                          << sp.size() << ".");
    std::vector<int> ii = sp.find();
    if (m.is_scalar()) {
      (*this)(ii) = densify(m);
    } else {
      (*this)(ii) = densify(m(ii));
    }
  }

  template<typename Scalar>
  void Matrix<Scalar>::get_nz(Matrix<Scalar>& m, bool ind1, const Slice& kk) const {
    // Scalar
    if (kk.is_scalar(nnz())) {
      m = nonzeros().at(kk.scalar(nnz()));
      return;
    }

    // Fall back on IM
    get_nz(m, ind1, kk.all(nnz(), ind1));
  }

  template<typename Scalar>
  void Matrix<Scalar>::get_nz(Matrix<Scalar>& m, bool ind1, const Matrix<int>& kk) const {
    // Scalar
    if (kk.is_scalar(true)) {
      return get_nz(m, ind1, to_slice(kk, ind1));
    }

    // Get nonzeros of kk
    const std::vector<int>& k = kk.nonzeros();
    int sz = nnz();

    // Check bounds
    if (!inBounds(k, -sz+ind1, sz+ind1)) {
      casadi_error("get_nz[kk] out of bounds. Your kk contains "
                   << *std::min_element(k.begin(), k.end()) << " up to "
                   << *std::max_element(k.begin(), k.end())
                   << ", which is outside the range [" << -sz+ind1 << ","<< sz+ind1 <<  ").");
    }

    // If indexed matrix was a row/column vector, make sure that the result is too
    bool tr = (is_column() && kk.is_row()) || (is_row() && kk.is_column());

    // Copy nonzeros
    m = zeros(tr ? kk.sparsity().T() : kk.sparsity());
    for (int el=0; el<k.size(); ++el) {
      casadi_assert_message(!(ind1 && k[el]<=0), "Matlab is 1-based, but requested index " <<
                                                k[el] <<  ". Note that negative slices are" <<
                                                " disabled in the Matlab interface. " <<
                                                "Possibly you may want to use 'end'.");
      int k_el = k[el]-ind1;
      m->at(el) = nonzeros().at(k_el>=0 ? k_el : k_el+sz);
    }
  }

  template<typename Scalar>
  void Matrix<Scalar>::set_nz(const Matrix<Scalar>& m, bool ind1, const Slice& kk) {
    // Scalar
    if (kk.is_scalar(nnz())) {
      nonzeros().at(kk.scalar(nnz())) = m.scalar();
      return;
    }

    // Fallback on IM
    set_nz(m, ind1, kk.all(nnz(), ind1));
  }

  template<typename Scalar>
  void Matrix<Scalar>::set_nz(const Matrix<Scalar>& m, bool ind1, const Matrix<int>& kk) {
    // Scalar
    if (kk.is_scalar(true)) {
      return set_nz(m, ind1, to_slice(kk, ind1));
    }

    // Assert dimensions of assigning matrix
    if (kk.sparsity() != m.sparsity()) {
      if (m.is_scalar()) {
        // m scalar means "set all"
        if (!m.is_dense()) return; // Nothing to set
        return set_nz(Matrix<Scalar>(kk.sparsity(), m), ind1, kk);
      } else if (kk.size() == m.size()) {
        // Project sparsity if needed
        return set_nz(project(m, kk.sparsity()), ind1, kk);
      } else if (kk.size1() == m.size2() && kk.size2() == m.size1()
                 && std::min(m.size1(), m.size2()) == 1) {
        // m is transposed if necessary
        return set_nz(m.T(), ind1, kk);
      } else {
        // Error otherwise
        casadi_error("Dimension mismatch." << "lhs is " << kk.size()
                     << ", while rhs is " << m.size());
      }
    }

    // Get nonzeros
    const std::vector<int>& k = kk.nonzeros();
    int sz = nnz();

    // Check bounds
    if (!inBounds(k, -sz+ind1, sz+ind1)) {
      casadi_error("set_nz[kk] out of bounds. Your kk contains "
                   << *std::min_element(k.begin(), k.end()) << " up to "
                   << *std::max_element(k.begin(), k.end())
                   << ", which is outside the range [" << -sz+ind1 << ","<< sz+ind1 <<  ").");
    }

    // Set nonzeros, ignoring negative indices
    for (int el=0; el<k.size(); ++el) {
      casadi_assert_message(!(ind1 && k[el]<=0), "Matlab is 1-based, but requested index " <<
                                                k[el] <<  ". Note that negative slices are" <<
                                                " disabled in the Matlab interface. " <<
                                                "Possibly you may want to use 'end'.");
      int k_el = k[el]-ind1;
      nonzeros().at(k_el>=0 ? k_el : k_el+sz) = m->at(el);
    }
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::densify(const Matrix<Scalar>& x) {
    return densify(x, 0);
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::densify(const Matrix<Scalar>& x,
                                             const Matrix<Scalar>& val) {
    // Check argument
    casadi_assert(val.is_scalar());

    // Quick return if possible
    if (x.is_dense()) return x;

    // Get sparsity pattern
    int nrow = x.size1();
    int ncol = x.size2();
    const int* colind = x.colind();
    const int* row = x.row();
    auto it = x.nonzeros().cbegin();

    // New data vector
    std::vector<Scalar> d(nrow*ncol, val.scalar());

    // Copy nonzeros
    for (int cc=0; cc<ncol; ++cc) {
      for (int el=colind[cc]; el<colind[cc+1]; ++el) {
        d[cc*nrow + row[el]] = *it++;
      }
    }

    // Construct return matrix
    return Matrix<Scalar>(Sparsity::dense(x.size()), d);
  }

  template<typename Scalar>
  Matrix<Scalar>::Matrix() : sparsity_(Sparsity(0, 0)) {
  }

  template<typename Scalar>
  Matrix<Scalar>::Matrix(const Matrix<Scalar>& m) : sparsity_(m.sparsity_), nonzeros_(m.nonzeros_) {
  }

  template<typename Scalar>
  Matrix<Scalar>::Matrix(const std::vector<Scalar>& x) :
      sparsity_(Sparsity::dense(x.size(), 1)), nonzeros_(x) {
  }

  template<typename Scalar>
  Matrix<Scalar>& Matrix<Scalar>::operator=(const Matrix<Scalar>& m) {
    sparsity_ = m.sparsity_;
    nonzeros_ = m.nonzeros_;
    return *this;
  }

  template<typename Scalar>
  std::string Matrix<Scalar>::type_name() { return matrixName<Scalar>(); }

  template<typename Scalar>
  void Matrix<Scalar>::print_scalar(std::ostream &stream, bool trailing_newline) const {
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
      stream << scalar();
    }

    if (trailing_newline) stream << std::endl;
    stream << std::flush;
    stream.precision(precision);
    stream.width(width);
    stream.flags(flags);
  }

  template<typename Scalar>
  void Matrix<Scalar>::print_vector(std::ostream &stream, bool trailing_newline) const {
    casadi_assert_message(is_column(), "Not a vector");

    // Get components
    std::vector<std::string> nz, inter;
    print_split(nz, inter);

    // Print intermediate expressions
    for (int i=0; i<inter.size(); ++i)
      stream << "@" << (i+1) << "=" << inter[i] << ", ";
    inter.clear();

    // Access data structures
    const int* r = row();
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
        stream << nz.at(el++);
      } else {
        stream << "00";
      }
    }
    stream << "]";

    if (trailing_newline) stream << std::endl;
    stream << std::flush;
  }

  template<typename Scalar>
  void Matrix<Scalar>::print_dense(std::ostream &stream, bool trailing_newline) const {
    // Print as a single line
    bool oneliner=this->size1()<=1;

    // Get components
    std::vector<std::string> nz, inter;
    print_split(nz, inter);

    // Print intermediate expressions
    for (int i=0; i<inter.size(); ++i)
      stream << "@" << (i+1) << "=" << inter[i] << ", ";
    inter.clear();

    // Index counter for each column
    const int* cptr = this->colind();
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
          stream << nz.at(cind[cc]++);
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
  }

  template<typename Scalar>
  void Matrix<Scalar>::print_sparse(std::ostream &stream, bool trailing_newline) const {
    if (nnz()==0) {
      stream << "all zero sparse: " << size1() << "-by-" << size2();
    } else {
      // Print header
      stream << "sparse: " << size1() << "-by-" << size2() << ", " << nnz() << " nnz";

      // Get components
      std::vector<std::string> nz, inter;
      print_split(nz, inter);

      // Print intermediate expressions
      for (int i=0; i<inter.size(); ++i)
        stream << std::endl << " @" << (i+1) << "=" << inter[i] << ",";
      inter.clear();

      // Print nonzeros
      for (int cc=0; cc<size2(); ++cc) {
        for (int el=colind(cc); el<colind(cc+1); ++el) {
          int rr=row(el);
          stream << std::endl << " (" << rr << ", " << cc << ") -> " << nz.at(el);
          InterruptHandler::check();
        }
      }
    }
    if (trailing_newline) stream << std::endl;
    stream << std::flush;
  }

  template<typename Scalar>
  void Matrix<Scalar>::print_split(std::vector<std::string>& nz,
                                    std::vector<std::string>& inter) const {
    nz.resize(nnz());
    inter.resize(0);

    // Temporary
    std::stringstream ss;
    ss.precision(stream_precision_);
    ss.width(stream_width_);
    if (stream_scientific_) {
      ss.setf(std::ios::scientific);
    } else {
      ss.unsetf(std::ios::scientific);
    }

    // Print nonzeros
    for (int i=0; i<nz.size(); ++i) {
      ss.str(std::string());
      ss << nonzeros().at(i);
      nz[i] = ss.str();
    }
  }

  template<typename Scalar>
  void Matrix<Scalar>::print(std::ostream &stream, bool trailing_newline) const {
    if (is_empty()) {
      stream << "[]";
    } else if (numel()==1) {
      print_scalar(stream, false);
    } else if (is_column()) {
      print_vector(stream, false);
    } else if (std::max(size1(), size2())<=10 || static_cast<double>(nnz())/numel()>=0.5) {
      // if "small" or "dense"
      print_dense(stream, false);
    } else {
      print_sparse(stream, false);
    }
    if (trailing_newline) stream << std::endl;
  }

  template<typename Scalar>
  void Matrix<Scalar>::repr(std::ostream &stream, bool trailing_newline) const {
    stream << type_name() << "(";
    print(stream, false);
    stream << ")";
    if (trailing_newline) stream << std::endl;
    stream << std::flush;
  }

  template<typename Scalar>
  void Matrix<Scalar>::reserve(int nnz) {
    reserve(nnz, size2());
  }

  template<typename Scalar>
  void Matrix<Scalar>::reserve(int nnz, int ncol) {
    nonzeros().reserve(nnz);
  }

  template<typename Scalar>
  void Matrix<Scalar>::resize(int nrow, int ncol) {
    sparsity_.resize(nrow, ncol);
  }

  template<typename Scalar>
  void Matrix<Scalar>::clear() {
    sparsity_ = Sparsity(0, 0);
    nonzeros().clear();
  }

  template<typename Scalar>
  Matrix<Scalar>::Matrix(double val) :
      sparsity_(Sparsity::dense(1, 1)), nonzeros_(std::vector<Scalar>(1, val)) {
  }

  template<typename Scalar>
  Matrix<Scalar>::Matrix(const std::vector< std::vector<double> >& d) {
    // Get dimensions
    int nrow=d.size();
    int ncol=d.empty() ? 1 : d.front().size();

    // Assert consistency
    for (int rr=0; rr<nrow; ++rr) {
      casadi_assert_message(ncol==d[rr].size(),
        "Matrix<Scalar>::Matrix(const std::vector< std::vector<Scalar> >& d): "
        "shape mismatch" << std::endl
        << "Attempting to construct a matrix from a nested list." << std::endl
        << "I got convinced that the desired size is ("<< nrow << " x " << ncol
        << " ), but now I encounter a vector of size ("
        << d[rr].size() <<  " )" << std::endl);
    }

    // Form matrix
    sparsity_ = Sparsity::dense(nrow, ncol);
    nonzeros().resize(nrow*ncol);
    typename std::vector<Scalar>::iterator it=nonzeros_.begin();
    for (int cc=0; cc<ncol; ++cc) {
      for (int rr=0; rr<nrow; ++rr) {
        *it++ = d[rr][cc];
      }
    }
  }

  template<typename Scalar>
  Matrix<Scalar>::Matrix(const Sparsity& sp) : sparsity_(sp), nonzeros_(sp.nnz(), 1) {
  }

  template<typename Scalar>
  Matrix<Scalar>::Matrix(int nrow, int ncol) : sparsity_(nrow, ncol) {
  }

  template<typename Scalar>
  Matrix<Scalar>::Matrix(const std::pair<int, int>& rc) : sparsity_(rc) {
  }

  template<typename Scalar>
  Matrix<Scalar>::Matrix(const Sparsity& sp, const Scalar& val, bool dummy) :
      sparsity_(sp), nonzeros_(sp.nnz(), val) {
  }

  template<typename Scalar>
  Matrix<Scalar>::Matrix(const Sparsity& sp, const std::vector<Scalar>& d, bool dummy) :
      sparsity_(sp), nonzeros_(d) {
    casadi_assert_message(sp.nnz()==d.size(), "Size mismatch." << std::endl
                          << "You supplied a sparsity of " << sp.dim()
                          << ", but the supplied vector is of length " << d.size());
  }

  template<typename Scalar>
  Matrix<Scalar>::Matrix(const Sparsity& sp, const Matrix<Scalar>& d) {
    if (d.is_scalar()) {
      *this = Matrix<Scalar>(sp, d.scalar(), false);
    } else if (d.is_column() || d.size1()==1) {
      casadi_assert(sp.nnz()==d.numel());
      if (d.is_dense()) {
        *this = Matrix<Scalar>(sp, d.nonzeros(), false);
      } else {
        *this = Matrix<Scalar>(sp, densify(d).nonzeros(), false);
      }
    } else {
      casadi_error("Matrix(Sparsisty, Matrix): Only allowed for scalars and vectors");
    }
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::unary(int op, const Matrix<Scalar> &x) {
    // Return value
    Matrix<Scalar> ret = Matrix<Scalar>::zeros(x.sparsity());

    // Nonzeros
    std::vector<Scalar>& ret_data = ret.nonzeros();
    const std::vector<Scalar>& x_data = x.nonzeros();

    // Do the operation on all non-zero elements
    for (int el=0; el<x.nnz(); ++el) {
      casadi_math<Scalar>::fun(op, x_data[el], x_data[el], ret_data[el]);
    }

    // Check the value of the structural zero-entries, if there are any
    if (!x.is_dense() && !operation_checker<F0XChecker>(op)) {
      // Get the value for the structural zeros
      Scalar fcn_0;
      casadi_math<Scalar>::fun(op, 0, 0, fcn_0);
      if (!casadi_limits<Scalar>::is_zero(fcn_0)) { // Remove this if?
        ret = densify(ret, fcn_0);
      }
    }

    return ret;
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::operator-() const {
    return unary(OP_NEG, *this);
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::operator+() const {
    return *this;
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::mrdivide(const Matrix<Scalar>& a,
                                              const Matrix<Scalar>& b) {
    casadi_assert_message(a.is_scalar() || b.is_scalar(), "Not implemented");
    return a/b;
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::mldivide(const Matrix<Scalar>& a,
                                              const Matrix<Scalar>& b) {
    casadi_assert_message(a.is_scalar() || b.is_scalar(), "Not implemented");
    return b/a;
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::mpower(const Matrix<Scalar>& a,
                                            const Matrix<Scalar>& b) {
    casadi_assert_message(a.is_scalar() || b.is_scalar(), "Not implemented");
    return pow(a, b);
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::printme(const Matrix<Scalar>& y) const {
    return binary(OP_PRINTME, *this, y);
  }

  template<typename Scalar>
  void Matrix<Scalar>::erase(const std::vector<int>& rr, const std::vector<int>& cc, bool ind1) {
    // Erase from sparsity pattern
    std::vector<int> mapping = sparsity_.erase(rr, cc, ind1);

    // Update non-zero entries
    for (int k=0; k<mapping.size(); ++k)
      nonzeros()[k] = nonzeros()[mapping[k]];

    // Truncate nonzero vector
    nonzeros().resize(mapping.size());
  }

  template<typename Scalar>
  void Matrix<Scalar>::erase(const std::vector<int>& rr, bool ind1) {
    // Erase from sparsity pattern
    std::vector<int> mapping = sparsity_.erase(rr, ind1);

    // Update non-zero entries
    for (int k=0; k<mapping.size(); ++k)
      nonzeros()[k] = nonzeros()[mapping[k]];

    // Truncate nonzero vector
    nonzeros().resize(mapping.size());
  }

  template<typename Scalar>
  void Matrix<Scalar>::remove(const std::vector<int>& rr, const std::vector<int>& cc) {
    if (!inBounds(rr, size1())) {
      casadi_error("Remove(rr, cc) out of bounds. Your rr contains "
                   << *std::min_element(rr.begin(), rr.end()) << " up to "
                   << *std::max_element(rr.begin(), rr.end())
                   << ", which is outside of the matrix shape " << dim() << ".");
    }
    if (!inBounds(cc, size2())) {
      casadi_error("Remove(rr, cc) out of bounds. Your cc contains "
                   << *std::min_element(cc.begin(), cc.end()) << " up to "
                   << *std::max_element(cc.begin(), cc.end())
                   << ", which is outside of the matrix shape " << dim() << ".");
    }

    // Remove by performing a complementary slice
    std::vector<int> rrc = complement(rr, size1());
    std::vector<int> ccc = complement(cc, size2());

    Matrix<Scalar> ret = (*this)(rrc, ccc);

    operator=(ret);

  }

  template<typename Scalar>
  void Matrix<Scalar>::enlarge(int nrow, int ncol, const std::vector<int>& rr,
                                 const std::vector<int>& cc, bool ind1) {
    sparsity_.enlarge(nrow, ncol, rr, cc, ind1);
  }

  template<typename Scalar>
  void Matrix<Scalar>::sanity_check(bool complete) const {
    sparsity_.sanity_check(complete);
    if (nonzeros_.size()!=sparsity_.nnz()) {
      std::stringstream s;
      s << "Matrix is not sane. The following must hold:" << std::endl;
      s << "  nonzeros().size() = sparsity().nnz(), but got nonzeros().size()  = "
        << nonzeros_.size()
        << "   and sparsity().nnz() = "  << sparsity_.nnz() << std::endl;
      casadi_error(s.str());
    }
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::mtimes(const Matrix<Scalar> &x, const Matrix<Scalar> &y) {
    if (x.is_scalar() || y.is_scalar()) {
      // Use element-wise multiplication if at least one factor scalar
      return x*y;
    } else {
      Matrix<Scalar> z = Matrix<Scalar>::zeros(Sparsity::mtimes(x.sparsity(), y.sparsity()));
      return mac(x, y, z);
    }
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::mac(const Matrix<Scalar> &x,
                                         const Matrix<Scalar> &y,
                                         const Matrix<Scalar> &z) {
    if (x.is_scalar() || y.is_scalar()) {
      // Use element-wise multiplication if at least one factor scalar
      return z + x*y;
    }

    // Check matching dimensions
    casadi_assert_message(x.size2()==y.size1(),
                          "Matrix product with incompatible dimensions. Lhs is "
                          << x.dim() << " and rhs is " << y.dim() << ".");

    casadi_assert_message(y.size2()==z.size2(),
                          "Matrix addition with incompatible dimensions. Lhs is "
                          << mtimes(x, y).dim() << " and rhs is " << z.dim() << ".");

    casadi_assert_message(x.size1()==z.size1(),
                          "Matrix addition with incompatible dimensions. Lhs is "
                          << mtimes(x, y).dim() << " and rhs is " << z.dim() << ".");

    // Check if we can simplify the product
    if (x.is_identity()) {
      return y + z;
    } else if (y.is_identity()) {
      return x + z;
    } else if (x.is_zero() || y.is_zero()) {
      return z;
    } else {
      // Carry out the matrix product
      Matrix<Scalar> ret = z;
      std::vector<Scalar> work(x.size1());
      casadi_mtimes(x.ptr(), x.sparsity(), y.ptr(), y.sparsity(),
                    ret.ptr(), ret.sparsity(), get_ptr(work), false);
      return ret;
    }
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::
  _bilin(const Matrix<Scalar>& A, const Matrix<Scalar>& x,
         const Matrix<Scalar>& y) {
    return casadi_bilin(A.ptr(), A.sparsity(), x.ptr(), y.ptr());
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::
  _rank1(const Matrix<Scalar>& A, const Matrix<Scalar>& alpha,
         const Matrix<Scalar>& x, const Matrix<Scalar>& y) {
    Matrix<Scalar> ret = A;
    casadi_rank1(ret.ptr(), ret.sparsity(), *alpha.ptr(), x.ptr(), y.ptr());
    return ret;
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::T() const {
    // quick return if empty or scalar
    if ((size1()==0 && size2()==0) || is_scalar()) return *this;

    // Create the new sparsity pattern and the mapping
    std::vector<int> mapping;
    Sparsity s = sparsity().transpose(mapping);

    // create the return matrix
    Matrix<Scalar> ret = zeros(s);

    // Copy the content
    for (int i=0; i<mapping.size(); ++i)
      ret->at(i) = nonzeros().at(mapping[i]);

    return ret;
  }

  template<typename Scalar>
  const Scalar Matrix<Scalar>::scalar() const {
    // Make sure that the matrix is 1-by-1
    casadi_assert_message(is_scalar(), "Can only convert 1-by-1 matrices to scalars");

    // return zero or the nonzero element
    if (nnz()==1)
      return nonzeros()[0];
    else
      return casadi_limits<Scalar>::zero;
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::binary(int op,
                                            const Matrix<Scalar> &x,
                                            const Matrix<Scalar> &y) {
    if (x.numel()==1)
      return scalar_matrix(op, x, y);
    else if (y.numel()==1)
      return matrix_scalar(op, x, y);
    else
      return matrix_matrix(op, x, y);
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::scalar_matrix(int op,
                                                   const Matrix<Scalar> &x,
                                                   const Matrix<Scalar> &y) {
    // Return value
    Matrix<Scalar> ret = Matrix<Scalar>::zeros(y.sparsity());

    // Nonzeros
    std::vector<Scalar>& ret_data = ret.nonzeros();
    const std::vector<Scalar>& x_data = x.nonzeros();
    const Scalar& x_val = x_data.empty() ? casadi_limits<Scalar>::zero : x->front();
    const std::vector<Scalar>& y_data = y.nonzeros();

    // Do the operation on all non-zero elements
    for (int el=0; el<y.nnz(); ++el) {
      casadi_math<Scalar>::fun(op, x_val, y_data[el], ret_data[el]);
    }

    // Check the value of the structural zero-entries, if there are any
    if (!y.is_dense() && !operation_checker<FX0Checker>(op)) {
      // Get the value for the structural zeros
      Scalar fcn_0;
      casadi_math<Scalar>::fun(op, x_val, casadi_limits<Scalar>::zero, fcn_0);
      if (!casadi_limits<Scalar>::is_zero(fcn_0)) { // Remove this if?
        ret = densify(ret, fcn_0);
      }
    }

    return ret;
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::matrix_scalar(int op,
                                                   const Matrix<Scalar> &x,
                                                   const Matrix<Scalar> &y) {
    // Return value
    Matrix<Scalar> ret = Matrix<Scalar>::zeros(x.sparsity());

    // Nonzeros
    std::vector<Scalar>& ret_data = ret.nonzeros();
    const std::vector<Scalar>& x_data = x.nonzeros();
    const std::vector<Scalar>& y_data = y.nonzeros();
    const Scalar& y_val = y_data.empty() ? casadi_limits<Scalar>::zero : y->front();

    // Do the operation on all non-zero elements
    for (int el=0; el<x.nnz(); ++el) {
      casadi_math<Scalar>::fun(op, x_data[el], y_val, ret_data[el]);
    }

    // Check the value of the structural zero-entries, if there are any
    if (!x.is_dense() && !operation_checker<F0XChecker>(op)) {
      // Get the value for the structural zeros
      Scalar fcn_0;
      casadi_math<Scalar>::fun(op, casadi_limits<Scalar>::zero, y_val, fcn_0);
      if (!casadi_limits<Scalar>::is_zero(fcn_0)) { // Remove this if?
        ret = densify(ret, fcn_0);
      }
    }

    return ret;
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::matrix_matrix(int op,
                                                   const Matrix<Scalar> &x,
                                                   const Matrix<Scalar> &y) {

    if (!(x.size2() == y.size2() && x.size1() == y.size1())) {
      std::stringstream ss;
      casadi_math<Scalar>::print(op, ss, "lhs", "rhs");
      casadi_error("matrix_matrix: dimension mismatch in element-wise matrix operation "
                   << ss.str() <<"." << std::endl << "Left argument has shape " << x.dim()
                   << ", right has shape " << y.dim() << ". They should be equal.");
    }

    // Get the sparsity pattern of the result
    // (ignoring structural zeros giving rise to nonzero result)
    const Sparsity& x_sp = x.sparsity();
    const Sparsity& y_sp = y.sparsity();
    Sparsity r_sp = x_sp.combine(y_sp,
                                        operation_checker<F0XChecker>(op),
                                        operation_checker<FX0Checker>(op));

    // Return value
    Matrix<Scalar> r = zeros(r_sp);

    // Perform the operations elementwise
    if (x_sp==y_sp) {
      // Matching sparsities
      casadi_math<Scalar>::fun(op, get_ptr(x.nonzeros()), get_ptr(y.nonzeros()),
                                 get_ptr(r.nonzeros()), r_sp.nnz());
    } else if (y_sp==r_sp) {
      // Project first argument
      Matrix<Scalar> x_mod = x(r_sp);
      casadi_math<Scalar>::fun(op, get_ptr(x_mod.nonzeros()), get_ptr(y.nonzeros()),
                                 get_ptr(r.nonzeros()), r_sp.nnz());
    } else if (x_sp==r_sp) {
      // Project second argument
      Matrix<Scalar> y_mod = y(r_sp);
      casadi_math<Scalar>::fun(op, get_ptr(x.nonzeros()),
                                 get_ptr(y_mod.nonzeros()), get_ptr(r.nonzeros()), r_sp.nnz());
    } else {
      // Project both arguments
      Matrix<Scalar> x_mod = x(r_sp);
      Matrix<Scalar> y_mod = y(r_sp);
      casadi_math<Scalar>::fun(op, get_ptr(x_mod.nonzeros()), get_ptr(y_mod.nonzeros()),
                                 get_ptr(r.nonzeros()), r_sp.nnz());
    }

    // Handle structural zeros giving rise to nonzero result, e.g. cos(0) == 1
    if (!r.is_dense() && !operation_checker<F00Checker>(op)) {
      // Get the value for the structural zeros
      Scalar fcn_0;
      casadi_math<Scalar>::fun(op, casadi_limits<Scalar>::zero,
                                 casadi_limits<Scalar>::zero, fcn_0);
      r = densify(r, fcn_0);
    }

    return r;
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::triplet(const std::vector<int>& row,
                                             const std::vector<int>& col,
                                             const Matrix<Scalar>& d) {
    return triplet(row, col, d, *std::max_element(row.begin(), row.end()),
                   *std::max_element(col.begin(), col.end()));
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::triplet(const std::vector<int>& row,
                                             const std::vector<int>& col,
                                             const Matrix<Scalar>& d,
                                             const std::pair<int, int>& rc) {
    return triplet(row, col, d, rc.first, rc.second);
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::triplet(const std::vector<int>& row,
                                             const std::vector<int>& col,
                                             const Matrix<Scalar>& d,
                                             int nrow, int ncol) {
    casadi_assert_message(col.size()==row.size() && col.size()==d.nnz(),
                          "Argument error in Matrix<Scalar>::triplet(row, col, d): "
                          "supplied lists must all be of equal length, but got: "
                          << row.size() << ", " << col.size()  << " and " << d.nnz());
    std::vector<int> mapping;
    Sparsity sp = Sparsity::triplet(nrow, ncol, row, col, mapping);
    return Matrix<Scalar>(sp, d[mapping]);
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::eye(int n) {
    return Matrix<Scalar>::ones(Sparsity::diag(n));
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::inf(const Sparsity& sp) {
    casadi_assert_message(std::numeric_limits<Scalar>::has_infinity,
                          "Datatype cannot represent infinity");
    return Matrix<Scalar>(sp, std::numeric_limits<Scalar>::infinity(), false);
  }


  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::inf(const std::pair<int, int>& rc) {
    return inf(rc.first, rc.second);
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::inf(int nrow, int ncol) {
    return inf(Sparsity::dense(nrow, ncol));
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::nan(const Sparsity& sp) {
    casadi_assert_message(std::numeric_limits<Scalar>::has_quiet_NaN,
                          "Datatype cannot represent not-a-number");
    return Matrix<Scalar>(sp, std::numeric_limits<Scalar>::quiet_NaN(), false);
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::nan(const std::pair<int, int>& rc) {
    return nan(rc.first, rc.second);
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::nan(int nrow, int ncol) {
    return nan(Sparsity::dense(nrow, ncol));
  }

  template<typename Scalar>
  bool Matrix<Scalar>::is_regular() const {
    return casadi::is_regular(nonzeros_);
  }

  template<typename Scalar>
  bool Matrix<Scalar>::is_smooth() const {
    return true;
  }

  template<typename Scalar>
  size_t Matrix<Scalar>::element_hash() const {
    throw CasadiException("\"element_hash\" not defined for instantiation");
  }

  template<typename Scalar>
  bool Matrix<Scalar>::is_leaf() const {
    throw CasadiException("\"is_leaf\" not defined for instantiation");
  }

  template<typename Scalar>
  bool Matrix<Scalar>::is_commutative() const {
    throw CasadiException("\"is_commutative\" not defined for instantiation");
  }

  template<typename Scalar>
  bool Matrix<Scalar>::is_symbolic() const {
    return false;
  }

  template<typename Scalar>
  bool Matrix<Scalar>::is_valid_input() const {
    return false;
  }

  template<typename Scalar>
  bool Matrix<Scalar>::has_duplicates() {
    throw CasadiException("\"has_duplicates\" not defined for instantiation");
  }

  template<typename Scalar>
  void Matrix<Scalar>::resetInput() {
    throw CasadiException("\"resetInput\" not defined for instantiation");
  }

  template<typename Scalar>
  bool Matrix<Scalar>::is_integer() const {
    // Look for non-integers
    for (auto&& e : nonzeros()) if (!casadi_limits<Scalar>::is_integer(e)) return false;

    // Integer if reached this point
    return true;
  }

  template<typename Scalar>
  bool Matrix<Scalar>::is_constant() const {
    // Look for non-constants
    for (auto&& e : nonzeros()) if (!casadi_limits<Scalar>::is_constant(e)) return false;

    // Constant if we reach this point
    return true;
  }

  template<typename Scalar>
  bool Matrix<Scalar>::is_one() const {
    if (!is_dense()) return false;

    // Look for non-ones
    for (auto&& e : nonzeros()) if (!casadi_limits<Scalar>::is_one(e)) return false;

    return true;
  }

  template<typename Scalar>
  bool Matrix<Scalar>::is_minus_one() const {
    if (!is_dense()) return false;

    // Look for non-minus-ones
    for (auto&& e : nonzeros()) if (!casadi_limits<Scalar>::is_minus_one(e)) return false;

    return true;
  }

  template<typename Scalar>
  bool Matrix<Scalar>::is_zero() const {

    // Look for non-zeros
    for (auto&& e : nonzeros()) if (!casadi_limits<Scalar>::is_zero(e)) return false;

    return true;
  }

  template<typename Scalar>
  bool Matrix<Scalar>::is_identity() const {

    // Make sure that the matrix is diagonal
    if (!sparsity().is_diag()) return false;

    // Make sure that all entries are one
    for (auto&& e : nonzeros()) if (!casadi_limits<Scalar>::is_one(e)) return false;

    return true;
  }

  template<typename Scalar>
  bool Matrix<Scalar>::is_equal(const Matrix<Scalar> &x, const Matrix<Scalar> &y, int depth) {
    // Assert matching dimensions
    casadi_assert_message(x.size() == y.size(), "Dimension mismatch");

    // Project to union of patterns and call recursively if different sparsity
    if (x.sparsity() != y.sparsity()) {
      Sparsity sp = x.sparsity() + y.sparsity();
      return is_equal(project(x, sp), project(y, sp), depth);
    }

    // Check individual elements
    auto y_it = y.nonzeros().begin();
    for (auto&& e : x.nonzeros()) {
      if (!casadi_limits<Scalar>::is_equal(e, *y_it++, depth)) return false;
    }

    // True if reched this point
    return true;
  }

  template<typename Scalar>
  bool Matrix<Scalar>::has_zeros() const {
    // Check if the structural nonzero is known to be zero
    for (auto&& e : nonzeros()) if (casadi_limits<Scalar>::is_zero(e)) return true;

    // No known zeros amongst the structurally nonzero entries
    return false;
  }

  template<typename Scalar>
  std::string Matrix<Scalar>::name() const {
    throw CasadiException("\"name\" not defined for instantiation");
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::dep(int ch) const {
    throw CasadiException("\"dep\" not defined for instantiation");
  }

  template<typename Scalar>
  int Matrix<Scalar>::n_dep() const {
    throw CasadiException("\"n_dep\" not defined for instantiation");
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::project(const Matrix<Scalar>& x,
                                             const Sparsity& sp, bool intersect) {
    if (intersect) {
      return project(x, sp.intersect(x.sparsity()), false);
    } else {
      Matrix<Scalar> ret = Matrix<Scalar>::zeros(sp);
      std::vector<Scalar> w(x.size1());
      casadi_project(x.ptr(), x.sparsity(), ret.ptr(), sp, get_ptr(w));
      return ret;
    }
  }

  template<typename Scalar>
  void Matrix<Scalar>::setEqualityCheckingDepth(int eq_depth) {
    throw CasadiException("\"setEqualityCheckingDepth\" not defined for instantiation");
  }

  template<typename Scalar>
  int Matrix<Scalar>::getEqualityCheckingDepth() {
    throw CasadiException("\"getEqualityCheckingDepth\" not defined for instantiation");
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::det(const Matrix<Scalar>& x) {
    int n = x.size2();
    casadi_assert_message(n == x.size1(), "matrix must be square");

    // Trivial return if scalar
    if (x.is_scalar()) return x;

    // Trivial case 2 x 2
    if (n==2) return x(0, 0) * x(1, 1) - x(0, 1) * x(1, 0);

    // Return expression
    Matrix<Scalar> ret = 0;

    // Find out which is the best direction to expand along

    // Build up an IM with ones on the non-zeros
    Matrix<int> sp = IM::ones(x.sparsity());

    // Have a count of the nonzeros for each row
    Matrix<int> row_count = Matrix<int>::sum2(sp);

    // A blank row? determinant is structurally zero
    if (!row_count.is_dense()) return 0;

    // Have a count of the nonzeros for each col
    Matrix<int> col_count = Matrix<int>::sum1(sp).T();

    // A blank col? determinant is structurally zero
    if (!row_count.is_dense()) return 0;

    int min_row = std::distance(row_count.nonzeros().begin(),
                                std::min_element(row_count.nonzeros().begin(),
                                                 row_count.nonzeros().end()));
    int min_col = std::distance(col_count.nonzeros().begin(),
                                std::min_element(col_count.nonzeros().begin(),
                                                 col_count.nonzeros().end()));

    if (min_row <= min_col) {
      // Expand along row j
      int j = row_count.sparsity().row(min_row);

      Matrix<Scalar> row = x(j, Slice(0, n));

      std::vector< int > col_i = row.sparsity().get_col();

      for (int k=0; k<row.nnz(); ++k) {
        // Sum up the cofactors
        ret += row->at(k)*cofactor(x, col_i.at(k), j);
      }
      return ret;
    } else {
      // Expand along col i
      int i = col_count.sparsity().row(min_col);

      Matrix<Scalar> col = x(Slice(0, n), i);

      const int* row_i = col.row();

      for (int k=0; k<col.nnz(); ++k) {
        // Sum up the cofactors
        ret += col->at(k)*cofactor(x, i, row_i[k]);
      }
      return ret;
    }

  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::sum2(const Matrix<Scalar>& x) {
    return mtimes(x, Matrix<Scalar>::ones(x.size2(), 1));
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::sum1(const Matrix<Scalar>& x) {
    return mtimes(Matrix<Scalar>::ones(1, x.size1()), x);
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::getMinor(const Matrix<Scalar>& x,
                                              int i, int j) {
    int n = x.size2();
    casadi_assert_message(n == x.size1(), "getMinor: matrix must be square");

    // Trivial return if scalar
    if (n==1) return 1;

    // Remove col i and row j
    Matrix<Scalar> M = Matrix<Scalar>(n-1, n-1);

    std::vector<int> col = x.sparsity().get_col();
    const int* row = x.sparsity().row();

    for (int k=0; k<x.nnz(); ++k) {
      int i1 = col[k];
      int j1 = row[k];

      if (i1 == i || j1 == j) continue;

      int i2 = (i1<i)?i1:i1-1;
      int j2 = (j1<j)?j1:j1-1;

      M(j2, i2) = x(j1, i1);
    }
    return det(M);
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::cofactor(const Matrix<Scalar>& A, int i, int j) {

    // Calculate the i, j minor
    Matrix<Scalar> minor_ij = getMinor(A, i, j);
    // Calculate the cofactor
    int sign_i = 1-2*((i+j) % 2);

    return sign_i * minor_ij;
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::adj(const Matrix<Scalar>& x) {
    int n = x.size2();
    casadi_assert_message(n == x.size1(), "adj: matrix must be square");

    // Temporary placeholder
    Matrix<Scalar> temp;

    // Cofactor matrix
    Matrix<Scalar> C = Matrix<Scalar>(n, n);
    for (int i=0; i<n; ++i)
      for (int j=0; j<n; ++j) {
        temp = cofactor(x, i, j);
        if (!temp.is_zero()) C(j, i) = temp;
      }

    return C.T();
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::inv(const Matrix<Scalar>& x) {
    // laplace formula
    return adj(x)/det(x);
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::reshape(const Matrix<Scalar>& x, int nrow, int ncol) {
    Sparsity sp = Sparsity::reshape(x.sparsity(), nrow, ncol);
    return Matrix<Scalar>(sp, x.nonzeros(), false);
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::reshape(const Matrix<Scalar>& x, const Sparsity& sp) {
    // quick return if already the right shape
    if (sp==x.sparsity()) return x;

    // make sure that the patterns match
    casadi_assert(sp.isReshape(x.sparsity()));

    return Matrix<Scalar>(sp, x.nonzeros(), false);
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::trace(const Matrix<Scalar>& x) {
    casadi_assert_message(x.is_square(), "trace: must be square");
    Scalar res=0;
    const Scalar* d=x.ptr();
    int size2 = x.size2();
    const int *colind=x.colind(), *row=x.row();
    for (int c=0; c<size2; c++) {
      for (int k=colind[c]; k!=colind[c+1]; ++k) {
        if (row[k]==c) {
          res += d[k];
        }
      }
    }
    return res;
  }

  template<typename Scalar>
  Matrix<Scalar>
  Matrix<Scalar>::blockcat(const std::vector< std::vector<Matrix<Scalar> > > &v) {
    std::vector< Matrix<Scalar> > ret;
    for (int i=0; i<v.size(); ++i)
      ret.push_back(horzcat(v[i]));
    return vertcat(ret);
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::horzcat(const std::vector<Matrix<Scalar> > &v) {
    // Concatenate sparsity patterns
    std::vector<Sparsity> sp(v.size());
    for (int i=0; i<v.size(); ++i) sp[i] = v[i].sparsity();
    Matrix<Scalar> ret = zeros(Sparsity::horzcat(sp));

    // Copy nonzeros
    auto i=ret->begin();
    for (auto&& j : v) {
      std::copy(j->begin(), j->end(), i);
      i += j.nnz();
    }
    return ret;
  }

  template<typename Scalar>
  std::vector<Matrix<Scalar> >
  Matrix<Scalar>::horzsplit(const Matrix<Scalar>& x, const std::vector<int>& offset) {
    // Split up the sparsity pattern
    std::vector<Sparsity> sp = Sparsity::horzsplit(x.sparsity(), offset);

    // Return object
    std::vector<Matrix<Scalar> > ret;
    ret.reserve(sp.size());

    // Copy data
    auto i=x.nonzeros().begin();
    for (auto&& j : sp) {
      auto i_next = i + j.nnz();
      ret.push_back(Matrix<Scalar>(j, std::vector<Scalar>(i, i_next), false));
      i = i_next;
    }

    // Return the assembled matrix
    casadi_assert(i==x.nonzeros().end());
    return ret;
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::vertcat(const std::vector<Matrix<Scalar> > &v) {
    std::vector<Matrix<Scalar> > vT(v.size());
    for (int i=0; i<v.size(); ++i) vT[i] = v[i].T();
    return horzcat(vT).T();
  }

  template<typename Scalar>
  std::vector< Matrix<Scalar> >
  Matrix<Scalar>::vertsplit(const Matrix<Scalar>& x, const std::vector<int>& offset) {
    std::vector< Matrix<Scalar> > ret = horzsplit(x.T(), offset);
    for (auto&& e : ret) e = e.T();
    return ret;
  }

  template<typename Scalar>
  std::vector< Matrix<Scalar> >
  Matrix<Scalar>::diagsplit(const Matrix<Scalar>& x, const std::vector<int>& offset1,
                              const std::vector<int>& offset2) {
    // Consistency check
    casadi_assert(offset1.size()>=1);
    casadi_assert(offset1.front()==0);
    casadi_assert(offset1.back()==x.size1());
    casadi_assert(isMonotone(offset1));

    // Consistency check
    casadi_assert(offset2.size()>=1);
    casadi_assert(offset2.front()==0);
    casadi_assert(offset2.back()==x.size2());
    casadi_assert(isMonotone(offset2));

    // Number of outputs
    int n = offset1.size()-1;

    // Return value
    std::vector< Matrix<Scalar> > ret;

    // Caveat: this is a very silly implementation
    for (int i=0; i<n; ++i) {
      ret.push_back(x(Slice(offset1[i], offset1[i+1]), Slice(offset2[i], offset2[i+1])));
    }

    return ret;
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::dot(const Matrix<Scalar> &x,
                                     const Matrix<Scalar> &y) {
    casadi_assert_message(x.size()==y.size(), "dot: Dimension mismatch");
    if (x.sparsity()!=y.sparsity()) {
      Sparsity sp = x.sparsity() * y.sparsity();
      return dot(project(x, sp), project(y, sp));
    }
    return casadi_dot(x.nnz(), x.ptr(), y.ptr());
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::all(const Matrix<Scalar>& x) {
    if (!x.is_dense()) return false;
    Scalar ret=1;
    for (int i=0; i<x.nnz(); ++i) {
      ret = ret && x->at(i)==1;
    }
    return ret;
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::any(const Matrix<Scalar>& x) {
    if (!x.is_dense()) return false;
    Scalar ret=0;
    for (int i=0; i<x.nnz(); ++i) {
      ret = ret || x->at(i)==1;
    }
    return ret;
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::norm_1(const Matrix<Scalar>& x) {
    return dot(fabs(x), Matrix<Scalar>::ones(x.sparsity()));
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::norm_2(const Matrix<Scalar>& x) {
    if (x.is_vector()) {
      return norm_F(x);
    } else {
      casadi_error("2-norms currently only supported for vectors. "
                   "Did you intend to calculate a Frobenius norms (norm_F)?");
    }
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::norm_F(const Matrix<Scalar>& x) {
    return sqrt(sum_square(x));
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::norm_inf(const Matrix<Scalar>& x) {
    // Get largest element by absolute value
    Matrix<Scalar> s = 0;
    for (auto i=x.nonzeros().begin(); i!=x.nonzeros().end(); ++i) {
      s = fmax(s, fabs(Matrix<Scalar>(*i)));
    }
    return s;
  }

  template<typename Scalar>
  void Matrix<Scalar>::qr(const Matrix<Scalar>& A,
                            Matrix<Scalar>& Q, Matrix<Scalar> &R) {
    // The following algorithm is taken from J. Demmel:
    // Applied Numerical Linear Algebra (algorithm 3.1.)
    casadi_assert_message(A.size1()>=A.size2(), "qr: fewer rows than columns");

    // compute Q and R column by column
    Q = R = Matrix<Scalar>();
    for (int i=0; i<A.size2(); ++i) {
      // Initialize qi to be the i-th column of *this
      Matrix<Scalar> ai = A(Slice(), i);
      Matrix<Scalar> qi = ai;
      // The i-th column of R
      Matrix<Scalar> ri = Matrix<Scalar>(A.size2(), 1);

      // subtract the projection of qi in the previous directions from ai
      for (int j=0; j<i; ++j) {

        // Get the j-th column of Q
        Matrix<Scalar> qj = Q(Slice(), j);

        ri(j, 0) = mtimes(qi.T(), qj); // Modified Gram-Schmidt
        // ri[j] = dot(qj, ai); // Classical Gram-Schmidt

        // Remove projection in direction j
        if (ri.has_nz(j, 0))
          qi -= ri(j, 0) * qj;
      }

      // Normalize qi
      ri(i, 0) = norm_2(qi);
      qi /= ri(i, 0);

      // Update R and Q
      Q = Matrix<Scalar>::horzcat({Q, qi});
      R = Matrix<Scalar>::horzcat({R, ri});
    }
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::nullspace(const Matrix<Scalar>& A) {
    Matrix<Scalar> X = A;
    int n = X.size1();
    int m = X.size2();
    casadi_assert_message(m>=n, "nullspace(): expecting a flat matrix (more columns than rows), "
                          "but got " << X.dim() << ".");

    Matrix<Scalar> seed = DM::eye(m)(Slice(0, m), Slice(n, m));

    std::vector< Matrix<Scalar> > us;
    std::vector< Matrix<Scalar> > betas;

    Matrix<Scalar> beta;

    for (int i=0;i<n;++i) {
      Matrix<Scalar> x = X(i, Slice(i, m));
      Matrix<Scalar> u = Matrix<Scalar>(x);
      Matrix<Scalar> sigma = sqrt(sum2(x*x));
      const Matrix<Scalar>& x0 = x(0, 0);
      u(0, 0) = 1;

      Matrix<Scalar> b = -copysign(sigma, x0);

      u(Slice(0), Slice(1, m-i))*= 1/(x0-b);
      beta = 1-x0/b;

      X(Slice(i, n), Slice(i, m)) -=
        beta*mtimes(mtimes(X(Slice(i, n), Slice(i, m)), u.T()), u);
      us.push_back(u);
      betas.push_back(beta);
    }

    for (int i=n-1;i>=0;--i) {
      seed(Slice(i, m), Slice(0, m-n)) -=
        betas[i]*mtimes(us[i].T(), mtimes(us[i], seed(Slice(i, m), Slice(0, m-n))));
    }

    return seed;

  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::chol(const Matrix<Scalar>& A) {
    // Cholesky-Banachiewicz algorithm
    // Naive, dense implementation O(n^3)

    // check dimensions
    casadi_assert_message(A.size1() == A.size2(), "Cholesky decomposition requires square matrix."
                                              "Got " << A.dim() << " instead.");

    Matrix<Scalar> ret = Matrix<Scalar>(Sparsity::lower(A.size1()));

    for (int i=0; i<A.size1(); ++i) { // loop over rows
      for (int j=0; j<i; ++j) {
        // Loop over columns before diagonal
        Matrix<Scalar> sum=0;
        for (int k=0;k<j;++k) {
          sum += ret(i, k)*ret(j, k);
        }
        ret(i, j) = (A(i, j)-sum)/ret(j, j);
      }

      // Treat the diagonal element separately
      int j = i;
      Matrix<Scalar> sum = 0;
      for (int k=0;k<j;++k) {
        sum += pow(ret(j, k), 2);
      }
      ret(j, j) = sqrt(A(j, j)-sum);
    }
    return ret.T();
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::solve(const Matrix<Scalar>& a, const Matrix<Scalar>& b) {
    // check dimensions
    casadi_assert_message(a.size1() == b.size1(), "solve Ax=b: dimension mismatch: b has "
                          << b.size1() << " rows while A has " << a.size1() << ".");
    casadi_assert_message(a.size1() == a.size2(), "solve: A not square but " << a.dim());

    if (a.is_tril()) {
      // forward substitution if lower triangular
      Matrix<Scalar> x = b;
      const int*  Arow = a.row();
      const int*  Acolind = a.colind();
      const std::vector<Scalar> & Adata = a.nonzeros();
      for (int i=0; i<a.size2(); ++i) { // loop over columns forwards
        for (int k=0; k<b.size2(); ++k) { // for every right hand side
          if (!x.has_nz(i, k)) continue;
          x(i, k) /= a(i, i);
          for (int kk=Acolind[i+1]-1; kk>=Acolind[i] && Arow[kk]>i; --kk) {
            int j = Arow[kk];
            x(j, k) -= Adata[kk]*x(i, k);
          }
        }
      }
      return x;
    } else if (a.is_triu()) {
      // backward substitution if upper triangular
      Matrix<Scalar> x = b;
      const int*  Arow = a.row();
      const int*  Acolind = a.colind();
      const std::vector<Scalar> & Adata = a.nonzeros();
      for (int i=a.size2()-1; i>=0; --i) { // loop over columns backwards
        for (int k=0; k<b.size2(); ++k) { // for every right hand side
          if (!x.has_nz(i, k)) continue;
          x(i, k) /= a(i, i);
          for (int kk=Acolind[i]; kk<Acolind[i+1] && Arow[kk]<i; ++kk) {
            int j = Arow[kk];
            x(j, k) -= Adata[kk]*x(i, k);
          }
        }
      }
      return x;
    } else if (a.has_zeros()) {

      // If there are structurally nonzero entries that are known to be zero,
      // remove these and rerun the algorithm
      return solve(sparsify(a), b);

    } else {

      // Make a BLT transformation of A
      std::vector<int> rowperm, colperm, rowblock, colblock, coarse_rowblock, coarse_colblock;
      a.sparsity().btf(rowperm, colperm, rowblock, colblock,
                       coarse_rowblock, coarse_colblock);

      // Permute the right hand side
      Matrix<Scalar> bperm = b(rowperm, Slice());

      // Permute the linear system
      Matrix<Scalar> Aperm = a(rowperm, colperm);

      // Solution
      Matrix<Scalar> xperm;

      // Solve permuted system
      if (Aperm.is_tril()) {

        // Forward substitution if lower triangular
        xperm = solve(Aperm, bperm);

      } else if (a.size2()<=3) {

        // Form inverse by minor expansion and multiply if very small (up to 3-by-3)
        xperm = mtimes(inv(Aperm), bperm);

      } else {

        // Make a QR factorization
        Matrix<Scalar> Q, R;
        qr(Aperm, Q, R);

        // Solve the factorized system (note that solve will now be fast since it is triangular)
        xperm = solve(R, mtimes(Q.T(), bperm));
      }

      // get the inverted column permutation
      std::vector<int> inv_colperm(colperm.size());
      for (int k=0; k<colperm.size(); ++k)
        inv_colperm[colperm[k]] = k;

      // Permute back the solution and return
      Matrix<Scalar> x = xperm(inv_colperm, Slice());
      return x;
    }
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::
  solve(const Matrix<Scalar>& a, const Matrix<Scalar>& b,
           const std::string& lsolver, const Dict& dict) {
    throw CasadiException("\"solve\" not defined for instantiation");
    return Matrix<Scalar>();
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::pinv(const Matrix<Scalar>& A) {
    if (A.size2()>=A.size1()) {
      return solve(mtimes(A, A.T()), A).T();
    } else {
      return solve(mtimes(A.T(), A), A.T());
    }
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::
  pinv(const Matrix<Scalar>& A, const std::string& lsolver, const Dict& dict) {
    throw CasadiException("\"solve\" not defined for instantiation");
    return Matrix<Scalar>();
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::kron(const Matrix<Scalar>& a, const Matrix<Scalar>& b) {
    const Sparsity &a_sp = a.sparsity();
    Matrix<Scalar> filler = Matrix<Scalar>(b.size());
    std::vector< std::vector< Matrix<Scalar> > >
      blocks(a.size1(), std::vector< Matrix<Scalar> >(a.size2(), filler));
    for (int i=0; i<a.size1(); ++i) {
      for (int j=0; j<a.size2(); ++j) {
        int k = a_sp.get_nz(i, j);
        if (k!=-1) {
          blocks[i][j] = a[k]*b;
        }
      }
    }
    return blockcat(blocks);
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::diag(const Matrix<Scalar>& A) {
    // Nonzero mapping
    std::vector<int> mapping;
    // Get the sparsity
    Sparsity sp = A.sparsity().get_diag(mapping);

    Matrix<Scalar> ret = zeros(sp);

    for (int k=0; k<mapping.size(); k++) ret[k] = A[mapping[k]];
    return ret;
  }

  /** \brief   Construct a matrix with given block on the diagonal */
  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::diagcat(const std::vector< Matrix<Scalar> > &A) {
    std::vector<Scalar> data;

    std::vector<Sparsity> sp;
    for (int i=0;i<A.size();++i) {
      data.insert(data.end(), A[i].nonzeros().begin(), A[i].nonzeros().end());
      sp.push_back(A[i].sparsity());
    }

    return Matrix<Scalar>(Sparsity::diagcat(sp), data, false);
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::unite(const Matrix<Scalar>& A, const Matrix<Scalar>& B) {
    // Join the sparsity patterns
    std::vector<unsigned char> mapping;
    Sparsity sp = A.sparsity().unite(B.sparsity(), mapping);

    // Create return matrix
    Matrix<Scalar> ret = zeros(sp);

    // Copy sparsity
    int elA=0, elB=0;
    for (int k=0; k<mapping.size(); ++k) {
      if (mapping[k]==1) {
        ret.nonzeros()[k] = A.nonzeros()[elA++];
      } else if (mapping[k]==2) {
        ret.nonzeros()[k] = B.nonzeros()[elB++];
      } else {
        throw CasadiException("Pattern intersection not empty");
      }
    }

    casadi_assert(A.nnz()==elA);
    casadi_assert(B.nnz()==elB);

    return ret;
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::polyval(const Matrix<Scalar>& p, const Matrix<Scalar>& x) {
    casadi_assert_message(p.is_dense(), "polynomial coefficients vector must be dense");
    casadi_assert_message(p.is_vector() && p.nnz()>0, "polynomial coefficients must be a vector");
    Matrix<Scalar> ret = x;
    for (auto&& e : ret.nonzeros()) {
      e = casadi_polyval(p.ptr(), p.numel()-1, e);
    }
    return ret;
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::norm_inf_mul(const Matrix<Scalar>& x,
                                                  const Matrix<Scalar>& y) {
    casadi_assert_message(y.size1()==x.size2(), "Dimension error. Got " << x.dim()
                          << " times " << y.dim() << ".");

    // Allocate work vectors
    std::vector<Scalar> dwork(x.size1());
    std::vector<int> iwork(x.size1()+1+y.size2());

    // Call C runtime
    return casadi_norm_inf_mul(x.ptr(), x.sparsity(), y.ptr(), y.sparsity(),
                               get_ptr(dwork), get_ptr(iwork));
  }

  template<typename Scalar>
  void Matrix<Scalar>::expand(const Matrix<Scalar>& ex,
                                Matrix<Scalar> &weights, Matrix<Scalar>& terms) {
    throw CasadiException("\"expand\" not defined for instantiation");
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::pw_const(const Matrix<Scalar>& ex,
                                              const Matrix<Scalar>& tval,
                                              const Matrix<Scalar>& val) {
    throw CasadiException("\"pw_const\" not defined for instantiation");
    return Matrix<Scalar>();
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::pw_lin(const Matrix<Scalar>& ex,
                                            const Matrix<Scalar>& tval,
                                            const Matrix<Scalar>& val) {
    throw CasadiException("\"pw_lin\" not defined for instantiation");
    return Matrix<Scalar>();
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::if_else(const Matrix<Scalar> &cond,
                                             const Matrix<Scalar> &if_true,
                                             const Matrix<Scalar> &if_false,
                                             bool short_circuit) {
    return cond*if_true + !cond*if_false;
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::conditional(const Matrix<Scalar>& ind,
                                                 const std::vector<Matrix<Scalar> >& x,
                                                 const Matrix<Scalar>& x_default,
                                                 bool short_circuit) {
    Matrix<Scalar> ret = x_default;
    for (int k=0; k<x.size(); ++k) {
      ret = if_else(ind==k, x[k], ret, short_circuit);
    }
    return ret;
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::heaviside(const Matrix<Scalar>& x) {
    return (1+sign(x))/2;
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::rectangle(const Matrix<Scalar>& x) {
    return 0.5*(sign(x+0.5)-sign(x-0.5));
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::triangle(const Matrix<Scalar>& x) {
    return rectangle(x/2)*(1-fabs(x));
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::ramp(const Matrix<Scalar>& x) {
    return x*heaviside(x);
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::
  gauss_quadrature(const Matrix<Scalar> &f,
                   const Matrix<Scalar> &x, const Matrix<Scalar> &a,
                   const Matrix<Scalar> &b, int order) {
    return gauss_quadrature(f, x, a, b, order, Matrix<Scalar>());
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::gauss_quadrature(const Matrix<Scalar>& f,
                                                      const Matrix<Scalar>& x,
                                                      const Matrix<Scalar>& a,
                                                      const Matrix<Scalar>& b, int order,
                                                      const Matrix<Scalar>& w) {
    throw CasadiException("\"gauss_quadrature\" not defined for instantiation");
    return Matrix<Scalar>();
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::simplify(const Matrix<Scalar> &x) {
    return x;
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::substitute(const Matrix<Scalar>& ex,
                                                const Matrix<Scalar>& v,
                                                const Matrix<Scalar>& vdef) {
    throw CasadiException("\"substitute\" not defined for instantiation");
    return Matrix<Scalar>();
  }

  template<typename Scalar>
  std::vector<Matrix<Scalar> >
  Matrix<Scalar>::substitute(const std::vector<Matrix<Scalar> >& ex,
                               const std::vector<Matrix<Scalar> >& v,
                               const std::vector<Matrix<Scalar> >& vdef) {
    throw CasadiException("\"substitute\" not defined for instantiation");
    return std::vector<Matrix<Scalar> >();
  }

  template<typename Scalar>
  void Matrix<Scalar>::substitute_inplace(const std::vector<Matrix<Scalar> >& v,
                                           std::vector<Matrix<Scalar> >& vdef,
                                           std::vector<Matrix<Scalar> >& ex,
                                           bool reverse) {
    throw CasadiException("\"substitute_inplace\" not defined for instantiation");
  }

  template<typename Scalar>
  bool Matrix<Scalar>::depends_on(const Matrix<Scalar> &x, const Matrix<Scalar> &arg) {
    throw CasadiException("\"depends_on\" not defined for instantiation");
    return false;
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::jacobian(const Matrix<Scalar> &f,
                                              const Matrix<Scalar> &x,
                                              bool symmetric) {
    casadi_error("\"jacobian\" not defined for " + type_name());
    return Matrix<Scalar>();
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::gradient(const Matrix<Scalar> &f,
                                              const Matrix<Scalar> &x) {
    casadi_error("\"gradient\" not defined for " + type_name());
    return Matrix<Scalar>();
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::tangent(const Matrix<Scalar> &f,
                                                const Matrix<Scalar> &x) {
    casadi_error("\"tangent\" not defined for " + type_name());
    return Matrix<Scalar>();
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::hessian(const Matrix<Scalar> &f,
                                             const Matrix<Scalar> &x) {
    throw CasadiException("\"hessian\" not defined for instantiation");
    return Matrix<Scalar>();
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::hessian(const Matrix<Scalar> &f,
                                             const Matrix<Scalar> &x,
                                             Matrix<Scalar> &g) {
    casadi_error("\"hessian\" not defined for " + type_name());
    return Matrix<Scalar>();
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::jtimes(const Matrix<Scalar> &ex,
                                            const Matrix<Scalar> &arg,
                                            const Matrix<Scalar> &v,
                                            bool tr) {
    casadi_error("\"jtimes\" not defined for " + type_name());
    return Matrix<Scalar>();
  }

  template<typename Scalar>
  std::vector<bool>
  Matrix<Scalar>::nl_var(const Matrix<Scalar> &expr, const Matrix<Scalar> &var) {
    casadi_error("\"nl_var\" not defined for " + type_name());
    return std::vector<bool>();
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::taylor(const Matrix<Scalar>& f,
                                            const Matrix<Scalar>& x,
                                            const Matrix<Scalar>& a, int order) {
    throw CasadiException("\"taylor\" not defined for instantiation");
    return Matrix<Scalar>();
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::mtaylor(const Matrix<Scalar>& f,
                                             const Matrix<Scalar>& x,
                                             const Matrix<Scalar>& a, int order) {
    throw CasadiException("\"mtaylor\" not defined for instantiation");
    return Matrix<Scalar>();
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::mtaylor(const Matrix<Scalar>& f,
                                             const Matrix<Scalar>& x,
                                             const Matrix<Scalar>& a, int order,
                                             const std::vector<int>&order_contributions) {
    throw CasadiException("\"mtaylor\" not defined for instantiation");
    return Matrix<Scalar>();
  }

  template<typename Scalar>
  int Matrix<Scalar>::n_nodes(const Matrix<Scalar>& x) {
    throw CasadiException("\"n_nodes\" not defined for instantiation");
    return 0;
  }

  template<typename Scalar>
  std::string
  Matrix<Scalar>::print_operator(const Matrix<Scalar>& x,
                                   const std::vector<std::string>& args) {
    throw CasadiException("\"print_operator\" not defined for instantiation");
    return std::string();
  }

  template<typename Scalar>
  std::vector<Matrix<Scalar> > Matrix<Scalar>::symvar(const Matrix<Scalar>& x) {
    throw CasadiException("\"symvar\" not defined for instantiation");
    return std::vector<Matrix<Scalar> >();
  }

  template<typename Scalar>
  void Matrix<Scalar>::shared(std::vector<Matrix<Scalar> >& ex,
                                       std::vector<Matrix<Scalar> >& v,
                                       std::vector<Matrix<Scalar> >& vdef,
                                       const std::string& v_prefix,
                                       const std::string& v_suffix) {
    throw CasadiException("\"shared\" not defined for instantiation");
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::poly_coeff(const Matrix<Scalar>& f,
                                                const Matrix<Scalar>&x) {
    throw CasadiException("\"poly_coeff\" not defined for instantiation");
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::poly_roots(const Matrix<Scalar>& p) {
    throw CasadiException("\"poly_roots\" not defined for instantiation");
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::eig_symbolic(const Matrix<Scalar>& m) {
    throw CasadiException("\"eig_symbolic\" not defined for instantiation");
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::sparsify(const Matrix<Scalar>& x, double tol) {
    // Quick return if there are no entries to be removed
    bool remove_nothing = true;
    for (auto it=x.nonzeros().begin(); it!=x.nonzeros().end() && remove_nothing; ++it) {
      remove_nothing = !casadi_limits<Scalar>::isAlmostZero(*it, tol);
    }
    if (remove_nothing) return x;

    // Get the current sparsity pattern
    int size1 = x.size1();
    int size2 = x.size2();
    const int* colind = x.colind();
    const int* row = x.row();

    // Construct the new sparsity pattern
    std::vector<int> new_colind(1, 0), new_row;
    std::vector<Scalar> new_data;

    // Loop over the columns
    for (int cc=0; cc<size2; ++cc) {
      // Loop over existing nonzeros
      for (int el=colind[cc]; el<colind[cc+1]; ++el) {
        // If it is not known to be a zero
        if (!casadi_limits<Scalar>::isAlmostZero(x->at(el), tol)) {
          // Save the nonzero in its new location
          new_data.push_back(x->at(el));

          // Add to pattern
          new_row.push_back(row[el]);
        }
      }
      // Save the new column offset
      new_colind.push_back(new_row.size());
    }

    // Construct the sparsity pattern
    Sparsity sp(size1, size2, new_colind, new_row);

    // Construct matrix and return
    return Matrix<Scalar>(sp, new_data);
  }


  template<typename Scalar>
  std::vector<Matrix<Scalar> > Matrix<Scalar>::get_input(const Function& f) {
    throw CasadiException("\"get_input\" not defined for " + type_name());
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::jac(const Function& f, int iind, int oind,
                                         bool compact, bool symmetric) {
    throw CasadiException("\"jac\" not defined for " + type_name());
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::grad(const Function& f, int iind, int oind) {
    throw CasadiException("\"grad\" not defined for " + type_name());
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::tang(const Function& f, int iind, int oind) {
    throw CasadiException("\"tang\" not defined for " + type_name());
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::jac(const Function& f, const std::string& iname, int oind,
                                         bool compact, bool symmetric) {
    throw CasadiException("\"jac\" not defined for " + type_name());
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::jac(const Function& f, int iind,
                                         const std::string& oname,
                                         bool compact, bool symmetric) {
    throw CasadiException("\"jac\" not defined for " + type_name());
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::jac(const Function& f, const std::string& iname,
                                         const std::string& oname,
         bool compact, bool symmetric) {
    throw CasadiException("\"jac\" not defined for " + type_name());
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::grad(const Function& f, const std::string& iname, int oind) {
    throw CasadiException("\"grad\" not defined for " + type_name());
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::grad(const Function& f, int iind, const std::string& oname) {
    throw CasadiException("\"grad\" not defined for " + type_name());
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::grad(const Function& f, const std::string& iname,
                                          const std::string& oname) {
    throw CasadiException("\"grad\" not defined for " + type_name());
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::tang(const Function& f, const std::string& iname, int oind) {
    throw CasadiException("\"tang\" not defined for " + type_name());
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::tang(const Function& f, int iind, const std::string& oname) {
    throw CasadiException("'tang' not defined for " + type_name());
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::tang(const Function& f, const std::string& iname,
                                          const std::string& oname) {
    throw CasadiException("'tang' not defined for " + type_name());
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::hess(const Function& f,
                                          int iind, int oind) {
    throw CasadiException("'hess' not defined for " + type_name());
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::hess(const Function& f,
                                          const std::string& iname, int oind) {
    throw CasadiException("'hess' not defined for " + type_name());
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::hess(const Function& f, int iind,
                                          const std::string& oname) {
    throw CasadiException("'hess' not defined for " + type_name());
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::hess(const Function& f, const std::string& iname,
                                          const std::string& oname) {
    throw CasadiException("'hess' not defined for " + type_name());
  }

  template<typename Scalar>
  Matrix<Scalar>::operator double() const {
    casadi_assert(is_scalar());
    return static_cast<double>(scalar());
  }

  template<typename Scalar>
  template<typename A>
  std::vector<A> Matrix<Scalar>::get_nonzeros() const {
    std::vector<A> ret(nnz());
    auto r = ret.begin();
    for (auto&& e : nonzeros()) *r++ = static_cast<A>(e);
    return ret;
  }

  template<typename Scalar>
  template<typename A>
  Matrix<Scalar>::operator std::vector<A>() const {
    // Get sparsity pattern
    int size1 = this->size1(), size2 = this->size2();
    const int *colind = this->colind(), *row = this->row();
    // Copy the nonzeros
    auto it = nonzeros().begin();
    std::vector<A> ret(numel(), 0);
    for (int cc=0; cc<size2; ++cc) {
      for (int el=colind[cc]; el<colind[cc+1]; ++el) {
        ret[row[el] + cc*size1] = static_cast<A>(*it++);
      }
    }
    return ret;
  }

  template<typename Scalar>
  Matrix<Scalar>::operator int() const {
    casadi_assert(is_scalar());
    return static_cast<int>(scalar());
  }

  // Template specializations
  template<>
  CASADI_EXPORT Matrix<double> Matrix<double>::
  solve(const Matrix<double>& a, const Matrix<double>& b,
        const std::string& lsolver, const Dict& dict);

  template<>
  CASADI_EXPORT Matrix<double> Matrix<double>::
  pinv(const Matrix<double>& A, const std::string& lsolver, const Dict& dict);

  template<>
  CASADI_EXPORT bool Matrix<SXElem>::__nonzero__() const;

  // Specialize functions in GenericMatrix<SX> and SX
  template<> inline std::string matrixName<SXElem>() { return "SX"; }
  template<> SX GenericMatrix<SX>::sym(const std::string& name, const Sparsity& sp);
  template<> bool SX::is_regular() const;
  template<> bool SX::is_smooth() const;
  template<> bool SX::is_leaf() const;
  template<> bool SX::is_commutative() const;
  template<> bool SX::is_symbolic() const;
  template<> bool SX::is_valid_input() const;
  template<> bool SX::has_duplicates();
  template<> void SX::resetInput();
  template<> SX SX::dep(int ch) const;
  template<> int SX::n_dep() const;
  template<> std::string SX::name() const;
  template<> void SX::setEqualityCheckingDepth(int eq_depth);
  template<> int SX::getEqualityCheckingDepth();
  template<> size_t SX::element_hash() const;
  template<> void SX::expand(const SX& f, SX& weights, SX& terms);
  template<> SX SX::pw_const(const SX& t, const SX& tval, const SX& val);
  template<> SX SX::pw_lin(const SX& t, const SX &tval, const SX &val);
  template<> SX SX::if_else(const SX &cond, const SX &if_true, const SX &if_false,
                            bool short_circuit);
  template<> SX SX::gauss_quadrature(const SX& f, const SX &x, const SX &a,
                                     const SX &b, int order,
                                     const SX& w);
  template<> SX SX::simplify(const SX& x);
  template<> SX SX::substitute(const SX& ex, const SX& v, const SX& vdef);
  template<> std::vector<SX > SX::substitute(const std::vector<SX >& ex,
                                                const std::vector<SX >& v,
                                                const std::vector<SX >& vdef);
  template<> void SX::substitute_inplace(const std::vector<SX >& v,
                                        std::vector<SX >& vdef,
                                        std::vector<SX >& ex,
                                        bool reverse);
  template<> bool SX::depends_on(const SX &x, const SX &arg);
  template<> std::vector<SX > SX::symvar(const SX &x);
  template<> SX SX::jacobian(const SX &f, const SX &x, bool symmetric);
  template<> SX SX::gradient(const SX &f, const SX &x);
  template<> SX SX::tangent(const SX &f, const SX &x);
  template<> SX SX::hessian(const SX &f, const SX &x);
  template<> SX SX::hessian(const SX &f, const SX &x, SX &g);
  template<> SX SX::jtimes(const SX &ex, const SX &arg, const SX &v, bool tr);
  template<> std::vector<bool> SX::nl_var(const SX &expr, const SX &var);
  template<> SX SX::taylor(const SX& f, const SX& x, const SX& a, int order);
  template<> SX SX::mtaylor(const SX& f, const SX& x, const SX& a, int order);
  template<> SX SX::mtaylor(const SX& f, const SX& x, const SX& a, int order,
                            const std::vector<int>& order_contributions);
  template<> int SX::n_nodes(const SX& x);
  template<> std::string
  SX::print_operator(const SX& x, const std::vector<std::string>& args);
  template<> void SX::shared(std::vector<SX >& ex,
                                    std::vector<SX >& v,
                                    std::vector<SX >& vdef,
                                    const std::string& v_prefix,
                                    const std::string& v_suffix);
  template<> SX SX::poly_coeff(const SX& f, const SX& x);
  template<> SX SX::poly_roots(const SX& p);
  template<> SX SX::eig_symbolic(const SX& m);
  template<> void SX::print_split(std::vector<std::string>& nz,
                                 std::vector<std::string>& inter) const;

  template<> std::vector<SX> SX::get_input(const Function& f);

  template<> SX SX::jac(const Function& f, int iind, int oind,
                        bool compact, bool symmetric);
  template<> SX SX::jac(const Function& f, const std::string& iname, int oind,
                        bool compact, bool symmetric);
  template<> SX SX::jac(const Function& f, int iind, const std::string& oname,
                        bool compact, bool symmetric);
  template<> SX SX::jac(const Function& f, const std::string& iname, const std::string& oname,
                        bool compact, bool symmetric);

  template<> SX SX::grad(const Function& f, int iind, int oind);
  template<> SX SX::grad(const Function& f, const std::string& iname, int oind);
  template<> SX SX::grad(const Function& f, int iind, const std::string& oname);
  template<> SX SX::grad(const Function& f, const std::string& iname, const std::string& oname);

  template<> SX SX::tang(const Function& f, int iind, int oind);
  template<> SX SX::tang(const Function& f, const std::string& iname, int oind);
  template<> SX SX::tang(const Function& f, int iind, const std::string& oname);
  template<> SX SX::tang(const Function& f, const std::string& iname, const std::string& oname);

  template<> SX SX::hess(const Function& f, int iind, int oind);
  template<> SX SX::hess(const Function& f, const std::string& iname, int oind);
  template<> SX SX::hess(const Function& f, int iind, const std::string& oname);
  template<> SX SX::hess(const Function& f, const std::string& iname, const std::string& oname);

#ifndef CASADI_MATRIX_CPP
  // Templates instantiated in matrix.cpp
  extern template class Matrix<double>;
  extern template class Matrix<int>;
  extern template class Matrix<SXElem>;
#endif // CASADI_MATRIX_CPP

  template<typename MatType>
  MatType _jtimes(const MatType &ex, const MatType &arg, const MatType &v, bool tr) {
    Function f("tmp", {arg}, {ex});

    // Split up v
    std::vector<MatType> vv = horzsplit(v);

    // Make sure well-posed
    casadi_assert(vv.size() >= 1);
    casadi_assert(ex.is_column());
    casadi_assert(arg.is_column());
    if (tr) {
      casadi_assert(v.size1()==ex.size1());
    } else {
      casadi_assert(v.size1()==arg.size1());
    }

    // Assemble arguments and directional derivatives
    std::vector<MatType> argv = MatType::get_input(f);
    std::vector<MatType> resv = f(argv);
    std::vector<std::vector<MatType> > seed(vv.size()), sens;
    for (int dir=0; dir<vv.size(); ++dir) {
      seed[dir] = { vv[dir]};
    }

    // Calculate directional derivatives
    if (tr) {
      f.reverse(argv, resv, seed, sens);
    } else {
      f.forward(argv, resv, seed, sens);
    }

    // Get the results
    for (int dir=0; dir<vv.size(); ++dir) {
      vv[dir] = sens[dir].at(0);
    }
    return horzcat(vv);
  }

  template<typename MatType>
  std::vector<bool> _nl_var(const MatType &expr, const MatType &var) {
    // Create a function for calculating a forward-mode derivative
    MatType v = MatType::sym("v", var.sparsity());
    Function f("tmp", {var}, {jtimes(expr, var, v)});

    // Propagate sparsities backwards seeding all outputs
    std::vector<bvec_t> seed(f.nnz_out(0), 1);
    std::vector<bvec_t> sens(f.nnz_in(0), 0);
    f.rev({get_ptr(sens)}, {get_ptr(seed)});

    // Temporaries for evaluation
    std::vector<bool> ret(sens.size());
    std::copy(sens.begin(), sens.end(), ret.begin());
    return ret;
  }

} // namespace casadi

/// \endcond

#endif // CASADI_MATRIX_IMPL_HPP
