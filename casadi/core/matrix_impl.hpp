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

#include "dm.hpp"
#include "im.hpp"
#include "sx.hpp"

#include "sx_node.hpp"
#include "linsol.hpp"
#include "expm.hpp"
#include "serializing_stream.hpp"

using namespace std;

namespace casadi {
  template<typename Scalar>
  void Matrix<Scalar>::set_precision(casadi_int precision) { stream_precision_ = precision; }

  template<typename Scalar>
  void Matrix<Scalar>::set_width(casadi_int width) { stream_width_ = width; }

  template<typename Scalar>
  void Matrix<Scalar>::set_scientific(bool scientific) { stream_scientific_ = scientific; }

  template<typename Scalar>
  casadi_int Matrix<Scalar>::stream_precision_ = 6;
  template<typename Scalar>
  casadi_int Matrix<Scalar>::stream_width_ = 0;
  template<typename Scalar>
  bool Matrix<Scalar>::stream_scientific_ = false;

  template<typename Scalar>
  std::default_random_engine Matrix<Scalar>::rng_(
    // Seed with current time
    std::chrono::system_clock::now().time_since_epoch().count());

  template<typename Scalar>
  void Matrix<Scalar>::rng(casadi_int seed) {
    rng_.seed(seed);
  }

  template<typename Scalar>
  bool Matrix<Scalar>::__nonzero__() const {
    if (numel()!=1) {
      casadi_error("Only scalar Matrix could have a truth value, but you "
                   "provided a shape" + dim());
    }
    return nonzeros().at(0)!=0;
  }

  template<typename Scalar>
  void Matrix<Scalar>::get(Matrix<Scalar>& m, bool ind1,
                                const Slice& rr, const Slice& cc) const {
    // Both are scalar
    if (rr.is_scalar(size1()) && cc.is_scalar(size2())) {
      casadi_int k = sparsity().get_nz(rr.scalar(size1()), cc.scalar(size2()));
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
                                const Slice& rr, const Matrix<casadi_int>& cc) const {
    // Fall back on IM-IM
    get(m, ind1, rr.all(size1(), ind1), cc);
  }

  template<typename Scalar>
  void Matrix<Scalar>::get(Matrix<Scalar>& m, bool ind1,
                                const Matrix<casadi_int>& rr, const Slice& cc) const {
    // Fall back on IM-IM
    get(m, ind1, rr, cc.all(size2(), ind1));
  }

  template<typename Scalar>
  void Matrix<Scalar>::get(Matrix<Scalar>& m, bool ind1,
                             const Matrix<casadi_int>& rr, const Matrix<casadi_int>& cc) const {
    // Scalar
    if (rr.is_scalar(true) && cc.is_scalar(true)) {
      return get(m, ind1, to_slice(rr, ind1), to_slice(cc, ind1));
    }

    // Make sure dense vectors
    casadi_assert(rr.is_dense() && rr.is_vector(),
                          "Marix::get: First index must be a dense vector");
    casadi_assert(cc.is_dense() && cc.is_vector(),
                          "Marix::get: Second index must be a dense vector");

    // Get the sparsity pattern - does bounds checking
    std::vector<casadi_int> mapping;
    Sparsity sp = sparsity().sub(rr.nonzeros(), cc.nonzeros(), mapping, ind1);

    // Copy nonzeros
    m = Matrix<Scalar>::zeros(sp);
    for (casadi_int k=0; k<mapping.size(); ++k) m->at(k) = nonzeros().at(mapping[k]);
  }

  template<typename Scalar>
  void Matrix<Scalar>::get(Matrix<Scalar>& m, bool ind1, const Slice& rr) const {
    // Scalar
    if (rr.is_scalar(numel())) {
      casadi_int r = rr.scalar(numel());
      casadi_int k = sparsity().get_nz(r % size1(), r / size1());
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
  void Matrix<Scalar>::get(Matrix<Scalar>& m, bool ind1, const Matrix<casadi_int>& rr) const {
    // Scalar
    if (rr.is_scalar(true)) {
      return get(m, ind1, to_slice(rr, ind1));
    }

    // If the indexed matrix is dense, use nonzero indexing
    if (is_dense()) {
      return get_nz(m, ind1, rr);
    }

    // Get the sparsity pattern - does bounds checking
    std::vector<casadi_int> mapping;
    Sparsity sp = sparsity().sub(rr.nonzeros(), rr.sparsity(), mapping, ind1);

    // If indexed matrix was a row/column vector, make sure that the result is too
    bool tr = (is_column() && rr.is_row()) || (is_row() && rr.is_column());

    // Copy nonzeros
    m = Matrix<Scalar>::zeros(tr ? sp.T() : sp);
    for (casadi_int k=0; k<mapping.size(); ++k) m->at(k) = nonzeros().at(mapping[k]);
  }

  template<typename Scalar>
  void Matrix<Scalar>::get(Matrix<Scalar>& m, bool ind1, const Sparsity& sp) const {
    casadi_assert(size()==sp.size(),
                          "Shape mismatch. This matrix has shape "
                          + str(size()) + ", but supplied sparsity index has shape "
                          + str(sp.size()) + ".");
    m = project(*this, sp);
  }

  template<typename Scalar>
  void Matrix<Scalar>::set(const Matrix<Scalar>& m, bool ind1,
                             const Slice& rr, const Slice& cc) {
    // Both are scalar
    if (rr.is_scalar(size1()) && cc.is_scalar(size2()) && m.is_dense()) {
      casadi_int oldsize = sparsity_.nnz();
      casadi_int ind = sparsity_.add_nz(rr.scalar(size1()), cc.scalar(size2()));
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
                             const Slice& rr, const Matrix<casadi_int>& cc) {
    // Fall back on (IM, IM)
    set(m, ind1, rr.all(size1(), ind1), cc);
  }

  template<typename Scalar>
  void Matrix<Scalar>::set(const Matrix<Scalar>& m, bool ind1,
                                const Matrix<casadi_int>& rr, const Slice& cc) {
    // Fall back on (IM, IM)
    set(m, ind1, rr, cc.all(size2(), ind1));
  }

  template<typename Scalar>
  void Matrix<Scalar>::set(const Matrix<Scalar>& m, bool ind1,
                                const Matrix<casadi_int>& rr, const Matrix<casadi_int>& cc) {
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
    casadi_assert(rr.is_dense() && rr.is_column(),
                          "Matrix::set: First index not dense vector");
    casadi_assert(cc.is_dense() && cc.is_column(),
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
        casadi_error("Dimension mismatch. lhs is " + str(rr.size1()) + "-by-"
                     + str(cc.size1()) + ", while rhs is " + str(m.size()));
      }
    }

    // Dimensions
    casadi_int sz1 = size1(), sz2 = size2();

    // Report out-of-bounds
    casadi_assert_in_range(rr.nonzeros(), -sz1+ind1, sz1+ind1);
    casadi_assert_in_range(cc.nonzeros(), -sz2+ind1, sz2+ind1);

    // If we are assigning with something sparse, first remove existing entries
    if (!m.is_dense()) {
      erase(rr.nonzeros(), cc.nonzeros(), ind1);
    }

    // Collect all assignments
    IM el = IM::zeros(m.sparsity());
    for (casadi_int j=0; j<el.size2(); ++j) { // Loop over columns of m
      casadi_int this_j = cc->at(j) - ind1; // Corresponding column in this
      if (this_j<0) this_j += sz2;
      for (casadi_int k=el.colind(j); k<el.colind(j+1); ++k) { // Loop over rows of m
        casadi_int i = m.row(k);
        casadi_int this_i = rr->at(i) - ind1; // Corresponding row in this
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
      casadi_int r = rr.scalar(numel());
      casadi_int oldsize = sparsity_.nnz();
      casadi_int ind = sparsity_.add_nz(r % size1(), r / size1());
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
  void Matrix<Scalar>::set(const Matrix<Scalar>& m, bool ind1, const Matrix<casadi_int>& rr) {
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
        return set(project(m, sp), ind1, Matrix<casadi_int>::project(rr, sp));
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
        casadi_error("Dimension mismatch. lhs is " + str(rr.size())
                     + ", while rhs is " + str(m.size()));
      }
    }

    // Dimensions of this
    casadi_int sz1 = size1(), sz2 = size2(), sz = nnz(), nel = numel(), rrsz = rr.nnz();

    // Quick return if nothing to set
    if (rrsz==0) return;

    // Check bounds
    casadi_assert_in_range(rr.nonzeros(), -nel+ind1, nel+ind1);

    // Dense mode
    if (is_dense() && m.is_dense()) {
      return set_nz(m, ind1, rr);
    }

    // Construct new sparsity pattern
    std::vector<casadi_int> new_row =
      sparsity().get_row(), new_col=sparsity().get_col(), nz(rr.nonzeros());
    new_row.reserve(sz+rrsz);
    new_col.reserve(sz+rrsz);
    nz.reserve(rrsz);
    for (std::vector<casadi_int>::iterator i=nz.begin(); i!=nz.end(); ++i) {
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
    for (casadi_int i=0; i<nz.size(); ++i) {
      nonzeros().at(nz[i]) = m->at(i);
    }
  }

  template<typename Scalar>
  void Matrix<Scalar>::set(const Matrix<Scalar>& m, bool ind1, const Sparsity& sp) {
    casadi_assert(size()==sp.size(),
                          "set(Sparsity sp): shape mismatch. This matrix has shape "
                          + str(size()) + ", but supplied sparsity index has shape "
                          + str(sp.size()) + ".");
    std::vector<casadi_int> ii = sp.find();
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
  void Matrix<Scalar>::get_nz(Matrix<Scalar>& m, bool ind1, const Matrix<casadi_int>& kk) const {
    // Scalar
    if (kk.is_scalar(true)) {
      return get_nz(m, ind1, to_slice(kk, ind1));
    }

    // Get nonzeros of kk
    const std::vector<casadi_int>& k = kk.nonzeros();
    casadi_int sz = nnz();

    // Check bounds
    casadi_assert_in_range(k, -sz+ind1, sz+ind1);

    // If indexed matrix was a row/column vector, make sure that the result is too
    bool tr = (is_column() && kk.is_row()) || (is_row() && kk.is_column());

    // Copy nonzeros
    m = zeros(tr ? kk.sparsity().T() : kk.sparsity());
    for (casadi_int el=0; el<k.size(); ++el) {
      casadi_assert(!(ind1 && k[el]<=0), "Matlab is 1-based, but requested index "
                                                + str(k[el]) + ". Note that negative slices are"
                                                " disabled in the Matlab interface. "
                                                "Possibly you may want to use 'end'.");
      casadi_int k_el = k[el]-ind1;
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
  void Matrix<Scalar>::set_nz(const Matrix<Scalar>& m, bool ind1, const Matrix<casadi_int>& kk) {
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
        casadi_error("Dimension mismatch. lhs is " + str(kk.size())
                     + ", while rhs is " + str(m.size()));
      }
    }

    // Get nonzeros
    const std::vector<casadi_int>& k = kk.nonzeros();
    casadi_int sz = nnz();

    // Check bounds
    casadi_assert_in_range(k, -sz+ind1, sz+ind1);

    // Set nonzeros, ignoring negative indices
    for (casadi_int el=0; el<k.size(); ++el) {
      casadi_assert(!(ind1 && k[el]<=0),
        "Matlab is 1-based, but requested index " + str(k[el])
        +  ". Note that negative slices are disabled in the Matlab interface. "
           "Possibly you may want to use 'end'.");
      casadi_int k_el = k[el]-ind1;
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
    casadi_assert_dev(val.is_scalar());

    // Quick return if possible
    if (x.is_dense()) return x;

    // Get sparsity pattern
    casadi_int nrow = x.size1();
    casadi_int ncol = x.size2();
    const casadi_int* colind = x.colind();
    const casadi_int* row = x.row();
    auto it = x.nonzeros().cbegin();

    // New data vector
    std::vector<Scalar> d(nrow*ncol, val.scalar());

    // Copy nonzeros
    for (casadi_int cc=0; cc<ncol; ++cc) {
      for (casadi_int el=colind[cc]; el<colind[cc+1]; ++el) {
        d[cc*nrow + row[el]] = *it++;
      }
    }

    // Construct return matrix
    return Matrix<Scalar>(Sparsity::dense(x.size()), d);
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::cumsum(const Matrix<Scalar> &x, casadi_int axis) {
    if (axis==-1) axis = x.is_row();
    Matrix<Scalar> ret = x;
    if (axis==0) {
      for (casadi_int i=1;i<x.size1();++i)
        ret(i, Slice()) += ret(i-1, Slice());
    } else {
      for (casadi_int i=1;i<x.size2();++i)
        ret(Slice(), i) += ret(Slice(), i-1);
    }
    return ret;
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::einstein(
      const Matrix<Scalar>& A, const Matrix<Scalar>& B, const Matrix<Scalar>& C,
      const std::vector<casadi_int>& dim_a, const std::vector<casadi_int>& dim_b,
      const std::vector<casadi_int>& dim_c,
      const std::vector<casadi_int>& a, const std::vector<casadi_int>& b,
      const std::vector<casadi_int>& c) {
    std::vector<casadi_int> iter_dims;
    std::vector<casadi_int> strides_a;
    std::vector<casadi_int> strides_b;
    std::vector<casadi_int> strides_c;
    casadi_int n_iter = einstein_process(A, B, C, dim_a, dim_b, dim_c, a, b, c,
          iter_dims, strides_a, strides_b, strides_c);

    const std::vector<Scalar>& Av = A.nonzeros();
    const std::vector<Scalar>& Bv = B.nonzeros();

    Matrix<Scalar> ret = C;
    std::vector<Scalar>& Cv = ret.nonzeros();

    einstein_eval(n_iter, iter_dims, strides_a, strides_b, strides_c,
      get_ptr(Av), get_ptr(Bv), get_ptr(Cv));
    return ret;
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::einstein(const Matrix<Scalar>& A, const Matrix<Scalar>& B,
      const std::vector<casadi_int>& dim_a, const std::vector<casadi_int>& dim_b,
      const std::vector<casadi_int>& dim_c,
      const std::vector<casadi_int>& a, const std::vector<casadi_int>& b,
      const std::vector<casadi_int>& c) {
    return Matrix<Scalar>::einstein(A, B, Matrix<Scalar>::zeros(product(dim_c), 1),
      dim_a, dim_b, dim_c, a, b, c);
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
  void Matrix<Scalar>::print_scalar(std::ostream &stream) const {
    casadi_assert(numel()==1, "Not a scalar");

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
    stream << std::flush;
    stream.precision(precision);
    stream.width(width);
    stream.flags(flags);
  }

  template<typename Scalar>
  void Matrix<Scalar>::print_vector(std::ostream &stream, bool truncate) const {
    casadi_assert(is_column(), "Not a vector");

    // Get components
    std::vector<std::string> nz, inter;
    print_split(nz, inter);

    // Print intermediate expressions
    for (casadi_int i=0; i<inter.size(); ++i)
      stream << "@" << (i+1) << "=" << inter[i] << ", ";
    inter.clear();

    // Access data structures
    const casadi_int* row = this->row();
    casadi_int nnz = this->nnz();
    casadi_int size1 = this->size1();

    // No need to truncate if less than 1000 entries
    const casadi_int max_numel = 1000;
    if (truncate && size1<=max_numel) truncate=false;

    // Nonzero
    casadi_int el=0;

    // Loop over rows
    stream << "[";
    for (casadi_int rr=0; rr<size1; ++rr) {
      // String representation
      std::string s = el<nnz && rr==row[el] ? nz.at(el++) : "00";

      // Truncate?
      if (truncate && rr>=3 && rr<size1-3) {
        // Do not print
        if (rr==3) stream << ", ...";
      } else {
        // Print
        if (rr!=0) stream << ", ";
        stream << s;
      }
    }
    stream << "]" << std::flush;
  }

  template<typename Scalar>
  void Matrix<Scalar>::print_dense(std::ostream &stream, bool truncate) const {
    print_dense(stream, sparsity(), ptr(), truncate);
  }

  template<typename Scalar>
  void Matrix<Scalar>::print_sparse(std::ostream &stream, bool truncate) const {
    print_sparse(stream, sparsity(), ptr(), truncate);
  }

  template<typename Scalar>
  void Matrix<Scalar>::print_split(std::vector<std::string>& nz,
                                    std::vector<std::string>& inter) const {

    print_split(nnz(), ptr(), nz, inter);
  }

  template<typename Scalar>
  void Matrix<Scalar>::print_scalar(std::ostream &stream, const Scalar& e) {
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
    stream << e;
    stream << std::flush;

    stream.precision(precision);
    stream.width(width);
    stream.flags(flags);
  }

  template<typename Scalar>
  void Matrix<Scalar>::print_split(casadi_int nnz, const Scalar* nonzeros,
                                    std::vector<std::string>& nz,
                                    std::vector<std::string>& inter) {
    nz.resize(nnz);
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
    for (casadi_int i=0; i<nz.size(); ++i) {
      ss.str(std::string());
      ss << nonzeros[i];
      nz[i] = ss.str();
    }
  }

  template<typename Scalar>
  void Matrix<Scalar>::print_sparse(std::ostream &stream, const Sparsity& sp,
      const Scalar* nonzeros, bool truncate) {
    // Access data structures
    casadi_int size1 = sp.size1();
    casadi_int size2 = sp.size2();
    const casadi_int* colind = sp.colind();
    const casadi_int* row = sp.row();
    casadi_int nnz = sp.nnz();

    // Quick return if all zero sparse
    if (nnz==0) {
      stream << "all zero sparse: " << size1 << "-by-" << size2 << std::flush;
      return;
    }

    // Print header
    stream << "sparse: " << size1 << "-by-" << size2 << ", " << nnz << " nnz";

    // Get components
    std::vector<std::string> nz, inter;
    print_split(nnz, nonzeros, nz, inter);

    // Print intermediate expressions
    for (casadi_int i=0; i<inter.size(); ++i)
      stream << std::endl << " @" << (i+1) << "=" << inter[i] << ",";
    inter.clear();

    // No need to truncate if less than 1000 nonzeros
    const casadi_int max_nnz = 1000;
    if (truncate && nnz<=max_nnz) truncate=false;

    // Print nonzeros
    for (casadi_int cc=0; cc<size2; ++cc) {
      for (casadi_int el=colind[cc]; el<colind[cc+1]; ++el) {
        if (truncate && el>=3 && el<nnz-3) {
          if (el==3) stream << std::endl << " ...";
        } else {
          stream << std::endl << " (" << row[el] << ", " << cc << ") -> " << nz.at(el);
          InterruptHandler::check();
        }
      }
    }
    stream << std::flush;
  }

  template<typename Scalar>
  void Matrix<Scalar>::print_dense(std::ostream &stream, const Sparsity& sp,
      const Scalar* nonzeros, bool truncate) {
    // Get components
    std::vector<std::string> nz, inter;
    print_split(sp.nnz(), nonzeros, nz, inter);

    // Print intermediate expressions
    for (casadi_int i=0; i<inter.size(); ++i)
      stream << "@" << (i+1) << "=" << inter[i] << ", ";
    inter.clear();

    // Access data structures
    casadi_int size1 = sp.size1();
    casadi_int size2 = sp.size2();
    const casadi_int* colind = sp.colind();
    const casadi_int* row = sp.row();

    // No need to truncate if less than 1000 entries
    const casadi_int max_numel = 1000;
    if (truncate && size1*size2<=max_numel) truncate=false;

    // Truncate rows and/or columns
    bool truncate_rows = truncate && size1>=7;
    bool truncate_columns = truncate && size2>=7;

    // Index counter for each column
    std::vector<casadi_int> ind(colind, colind+size2+1);

    // Print as a single line?
    bool oneliner=size1<=1;

    // Loop over rows
    for (casadi_int rr=0; rr<size1; ++rr) {
      // Print row?
      bool print_row = !(truncate_rows && rr>=3 && rr<size1-3);

      // Beginning of row
      if (rr==0) {
        if (!oneliner) stream << std::endl;
        stream << "[[";
      } else if (print_row) {
        stream << " [";
      }

      // Loop over columns
      for (casadi_int cc=0; cc<size2; ++cc) {
        // String representation of element
        std::string s = ind[cc]<colind[cc+1] && row[ind[cc]]==rr
          ? nz.at(ind[cc]++) : "00";

        // Skip whole row?
        if (!print_row) continue;

        // Print column?
        bool print_column = !(truncate_columns && cc>=3 && cc<size2-3);

        // Print element
        if (print_column) {
          if (cc!=0) stream << ", ";
          stream << s;
        } else if (cc==3) {
          stream << ", ...";
        }
      }

      // End of row
      if (rr<size1-1) {
        if (print_row) {
          stream << "], ";
          if (!oneliner) stream << std::endl;
        } else if (rr==3) {
          stream << " ...," << std::endl;
        }
      } else {
        stream << "]]";
      }
    }
    stream << std::flush;
  }

  template<typename Scalar>
  void Matrix<Scalar>::to_file(const std::string& filename, const std::string& format) const {
    to_file(filename, sparsity(), ptr(), format);
  }

  template<typename Scalar>
  void Matrix<Scalar>::disp(std::ostream& stream, bool more) const {
    if (is_empty()) {
      stream << "[]";
    } else if (numel()==1) {
      print_scalar(stream);
    } else if (is_column()) {
      print_vector(stream);
    } else if (std::max(size1(), size2())<=10 ||
        static_cast<double>(nnz())/static_cast<double>(numel())>=0.5) {
      // if "small" or "dense"
      print_dense(stream);
    } else {
      print_sparse(stream);
    }
  }

  template<typename Scalar>
  std::string Matrix<Scalar>::get_str(bool more) const {
    std::stringstream ss;
    disp(ss, more);
    return ss.str();
  }

  template<typename Scalar>
  void Matrix<Scalar>::reserve(casadi_int nnz) {
    reserve(nnz, size2());
  }

  template<typename Scalar>
  void Matrix<Scalar>::reserve(casadi_int nnz, casadi_int ncol) {
    nonzeros().reserve(nnz);
  }

  template<typename Scalar>
  void Matrix<Scalar>::resize(casadi_int nrow, casadi_int ncol) {
    sparsity_.resize(nrow, ncol);
  }

  template<typename Scalar>
  void Matrix<Scalar>::clear() {
    sparsity_ = Sparsity(0, 0);
    nonzeros().clear();
  }

  template<typename Scalar>
  Matrix<Scalar>::Matrix(double val) :
      sparsity_(
        Sparsity::dense(1, 1)),
        nonzeros_(std::vector<Scalar>(1, static_cast<Scalar>(val))) {
  }

  template<typename Scalar>
  Matrix<Scalar>::Matrix(const std::vector< std::vector<double> >& d) {
    // Get dimensions
    casadi_int nrow=d.size();
    casadi_int ncol=d.empty() ? 1 : d.front().size();

    // Assert consistency
    for (casadi_int rr=0; rr<nrow; ++rr) {
      casadi_assert(ncol==d[rr].size(),
        "Shape mismatch.\n"
        "Attempting to construct a matrix from a nested list.\n"
        "I got convinced that the desired size is (" + str(nrow) + " x " + str(ncol)
        + " ), but now I encounter a vector of size (" + str(d[rr].size()) +  " )");
    }

    // Form matrix
    sparsity_ = Sparsity::dense(nrow, ncol);
    nonzeros().resize(nrow*ncol);
    typename std::vector<Scalar>::iterator it=nonzeros_.begin();
    for (casadi_int cc=0; cc<ncol; ++cc) {
      for (casadi_int rr=0; rr<nrow; ++rr) {
        *it++ = static_cast<Scalar>(d[rr][cc]);
      }
    }
  }

  template<typename Scalar>
  Matrix<Scalar>::Matrix(const Sparsity& sp) : sparsity_(sp), nonzeros_(sp.nnz(), 1) {
  }

  template<typename Scalar>
  Matrix<Scalar>::Matrix(casadi_int nrow, casadi_int ncol) : sparsity_(nrow, ncol) {
  }

  template<typename Scalar>
  Matrix<Scalar>::Matrix(const std::pair<casadi_int, casadi_int>& rc) : sparsity_(rc) {
  }

  template<typename Scalar>
  Matrix<Scalar>::Matrix(const Sparsity& sp, const Scalar& val, bool dummy) :
      sparsity_(sp), nonzeros_(sp.nnz(), val) {
  }

  template<typename Scalar>
  Matrix<Scalar>::Matrix(const Sparsity& sp, const std::vector<Scalar>& d, bool dummy) :
      sparsity_(sp), nonzeros_(d) {
    casadi_assert(sp.nnz()==d.size(), "Size mismatch.\n"
                          "You supplied a sparsity of " + sp.dim()
                          + ", but the supplied vector is of length " + str(d.size()));
  }

  template<typename Scalar>
  Matrix<Scalar>::Matrix(const Sparsity& sp, const Matrix<Scalar>& d) {
    if (d.is_scalar()) {
      *this = Matrix<Scalar>(sp, d.scalar(), false);
    } else if (sp.nnz()==0) {
      casadi_assert(d.nnz()==0,
        "You passed nonzeros (" + d.dim(true) +
        ") to the constructor of a fully sparse matrix (" + sp.dim(true) + ").");
      *this = Matrix<Scalar>(sp);
    } else if (d.is_column() || d.size1()==1) {
      casadi_assert_dev(sp.nnz()==d.numel());
      if (d.is_dense()) {
        *this = Matrix<Scalar>(sp, d.nonzeros(), false);
      } else {
        *this = Matrix<Scalar>(sp, densify(d).nonzeros(), false);
      }
    } else {
      casadi_error("Matrix(Sparsity, Matrix): Only allowed for scalars and vectors");
    }
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::unary(casadi_int op, const Matrix<Scalar> &x) {
    // Return value
    Matrix<Scalar> ret = Matrix<Scalar>::zeros(x.sparsity());

    // Nonzeros
    std::vector<Scalar>& ret_data = ret.nonzeros();
    const std::vector<Scalar>& x_data = x.nonzeros();

    // Do the operation on all non-zero elements
    for (casadi_int el=0; el<x.nnz(); ++el) {
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
  Matrix<Scalar> Matrix<Scalar>::mrdivide(const Matrix<Scalar>& b,
                                              const Matrix<Scalar>& a) {
    if (a.is_scalar() || b.is_scalar()) return b/a;
    return solve(a.T(), b.T()).T();
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::mldivide(const Matrix<Scalar>& a,
                                              const Matrix<Scalar>& b) {
    if (a.is_scalar() || b.is_scalar()) return b/a;
    return solve(a, b);
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::printme(const Matrix<Scalar>& y) const {
    return binary(OP_PRINTME, *this, y);
  }

  template<typename Scalar>
  void Matrix<Scalar>::erase(const std::vector<casadi_int>& rr,
      const std::vector<casadi_int>& cc, bool ind1) {
    // Erase from sparsity pattern
    std::vector<casadi_int> mapping = sparsity_.erase(rr, cc, ind1);

    // Update non-zero entries
    for (casadi_int k=0; k<mapping.size(); ++k)
      nonzeros()[k] = nonzeros()[mapping[k]];

    // Truncate nonzero vector
    nonzeros().resize(mapping.size());
  }

  template<typename Scalar>
  void Matrix<Scalar>::erase(const std::vector<casadi_int>& rr, bool ind1) {
    // Erase from sparsity pattern
    std::vector<casadi_int> mapping = sparsity_.erase(rr, ind1);

    // Update non-zero entries
    for (casadi_int k=0; k<mapping.size(); ++k)
      nonzeros()[k] = nonzeros()[mapping[k]];

    // Truncate nonzero vector
    nonzeros().resize(mapping.size());
  }

  template<typename Scalar>
  void Matrix<Scalar>::remove(const std::vector<casadi_int>& rr,
      const std::vector<casadi_int>& cc) {
    casadi_assert_bounded(rr, size1());
    casadi_assert_bounded(cc, size2());

    // Remove by performing a complementary slice
    std::vector<casadi_int> rrc = complement(rr, size1());
    std::vector<casadi_int> ccc = complement(cc, size2());

    Matrix<Scalar> ret = (*this)(rrc, ccc); // NOLINT(cppcoreguidelines-slicing)

    operator=(ret);

  }

  template<typename Scalar>
  void Matrix<Scalar>::enlarge(casadi_int nrow, casadi_int ncol, const std::vector<casadi_int>& rr,
                                 const std::vector<casadi_int>& cc, bool ind1) {
    sparsity_.enlarge(nrow, ncol, rr, cc, ind1);
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
    casadi_assert(x.size2()==y.size1(),
                          "Matrix product with incompatible dimensions. Lhs is "
                          + x.dim() + " and rhs is " + y.dim() + ".");

    casadi_assert(y.size2()==z.size2(),
                          "Matrix addition with incompatible dimensions. Lhs is "
                          + mtimes(x, y).dim() + " and rhs is " + z.dim() + ".");

    casadi_assert(x.size1()==z.size1(),
                          "Matrix addition with incompatible dimensions. Lhs is "
                          + mtimes(x, y).dim() + " and rhs is " + z.dim() + ".");

    // Check if we can simplify the product
    if (x.is_eye()) {
      return y + z;
    } else if (y.is_eye()) {
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
    std::vector<casadi_int> mapping;
    Sparsity s = sparsity().transpose(mapping);

    // create the return matrix
    Matrix<Scalar> ret = zeros(s);

    // Copy the content
    for (casadi_int i=0; i<mapping.size(); ++i)
      ret->at(i) = nonzeros().at(mapping[i]);

    return ret;
  }

  template<typename Scalar>
  const Scalar Matrix<Scalar>::scalar() const {
    // Make sure that the matrix is 1-by-1
    casadi_assert(is_scalar(), "Can only convert 1-by-1 matrices to scalars");

    // return zero or the nonzero element
    if (nnz()==1)
      return nonzeros()[0];
    else
      return casadi_limits<Scalar>::zero;
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::binary(casadi_int op,
                                            const Matrix<Scalar> &x,
                                            const Matrix<Scalar> &y) {
    if (x.is_scalar()) {
      return scalar_matrix(op, x, y);
    } else if (y.is_scalar()) {
      return matrix_scalar(op, x, y);
    } else {
      return matrix_matrix(op, x, y);
    }
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::
  scalar_matrix(casadi_int op, const Matrix<Scalar> &x, const Matrix<Scalar> &y) {
    if ( (operation_checker<FX0Checker>(op) && y.nnz()==0) ||
         (operation_checker<F0XChecker>(op) && x.nnz()==0))
            return Matrix<Scalar>::zeros(Sparsity(y.size()));

    // Return value
    Matrix<Scalar> ret = Matrix<Scalar>::zeros(y.sparsity());

    // Nonzeros
    std::vector<Scalar>& ret_data = ret.nonzeros();
    const std::vector<Scalar>& x_data = x.nonzeros();
    const Scalar& x_val = x_data.empty() ? casadi_limits<Scalar>::zero : x->front();
    const std::vector<Scalar>& y_data = y.nonzeros();

    // Do the operation on all non-zero elements
    for (casadi_int el=0; el<y.nnz(); ++el) {
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
  Matrix<Scalar> Matrix<Scalar>::
  matrix_scalar(casadi_int op, const Matrix<Scalar> &x, const Matrix<Scalar> &y) {

    if ( (operation_checker<FX0Checker>(op) && y.nnz()==0) ||
         (operation_checker<F0XChecker>(op) && x.nnz()==0))
            return Matrix<Scalar>::zeros(Sparsity(x.size()));

    // Return value
    Matrix<Scalar> ret = Matrix<Scalar>::zeros(x.sparsity());

    // Nonzeros
    std::vector<Scalar>& ret_data = ret.nonzeros();
    const std::vector<Scalar>& x_data = x.nonzeros();
    const std::vector<Scalar>& y_data = y.nonzeros();
    const Scalar& y_val = y_data.empty() ? casadi_limits<Scalar>::zero : y->front();

    // Do the operation on all non-zero elements
    for (casadi_int el=0; el<x.nnz(); ++el) {
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
  Matrix<Scalar> Matrix<Scalar>::
  matrix_matrix(casadi_int op, const Matrix<Scalar> &x, const Matrix<Scalar> &y) {
    // Check, correct dimensions
    if (x.size() != y.size()) {
      // x and y are horizontal multiples of each other?
      if (!x.is_empty() && !y.is_empty()) {
        if (x.size1() == y.size1() && x.size2() % y.size2() == 0) {
          return matrix_matrix(op, x, repmat(y, 1, x.size2() / y.size2()));
        } else if (y.size1() == x.size1() && y.size2() % x.size2() == 0) {
          return matrix_matrix(op, repmat(x, 1, y.size2() / x.size2()), y);
        }
      }
      // Dimension mismatch
      casadi_error("Dimension mismatch for " + casadi_math<Scalar>::print(op, "x", "y") +
                   ", x is " + x.dim() + ", while y is " + y.dim());
    }

    // Get the sparsity pattern of the result
    // (ignoring structural zeros giving rise to nonzero result)
    const Sparsity& x_sp = x.sparsity();
    const Sparsity& y_sp = y.sparsity();
    Sparsity r_sp = x_sp.combine(y_sp, operation_checker<F0XChecker>(op),
                                 operation_checker<FX0Checker>(op));

    // Return value
    Matrix<Scalar> r = zeros(r_sp);

    // Perform the operations elementwise
    if (x_sp==y_sp) {
      // Matching sparsities
      casadi_math<Scalar>::fun(op, x.ptr(), y.ptr(), r.ptr(), r_sp.nnz());
    } else if (y_sp==r_sp) {
      // Project first argument
      Matrix<Scalar> x_mod = x(r_sp);
      casadi_math<Scalar>::fun(op, x_mod.ptr(), y.ptr(), r.ptr(), r_sp.nnz());
    } else if (x_sp==r_sp) {
      // Project second argument
      Matrix<Scalar> y_mod = y(r_sp);
      casadi_math<Scalar>::fun(op, x.ptr(), y_mod.ptr(), r.ptr(), r_sp.nnz());
    } else {
      // Project both arguments
      Matrix<Scalar> x_mod = x(r_sp);
      Matrix<Scalar> y_mod = y(r_sp);
      casadi_math<Scalar>::fun(op, x_mod.ptr(), y_mod.ptr(), r.ptr(), r_sp.nnz());
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
  Matrix<Scalar> Matrix<Scalar>::triplet(const std::vector<casadi_int>& row,
                                             const std::vector<casadi_int>& col,
                                             const Matrix<Scalar>& d) {
    return triplet(row, col, d, *std::max_element(row.begin(), row.end()),
                   *std::max_element(col.begin(), col.end()));
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::triplet(const std::vector<casadi_int>& row,
                                             const std::vector<casadi_int>& col,
                                             const Matrix<Scalar>& d,
                                             const std::pair<casadi_int, casadi_int>& rc) {
    return triplet(row, col, d, rc.first, rc.second);
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::triplet(const std::vector<casadi_int>& row,
                                             const std::vector<casadi_int>& col,
                                             const Matrix<Scalar>& d,
                                             casadi_int nrow, casadi_int ncol) {
    casadi_assert(col.size()==row.size() && col.size()==d.nnz(),
                          "Argument error in Matrix<Scalar>::triplet(row, col, d): "
                          "supplied lists must all be of equal length, but got: "
                          + str(row.size()) + ", " + str(col.size()) + " and " + str(d.nnz()));
    std::vector<casadi_int> mapping;
    Sparsity sp = Sparsity::triplet(nrow, ncol, row, col, mapping, false);
    return Matrix<Scalar>(sp, d.nz(mapping));
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::eye(casadi_int n) {
    return Matrix<Scalar>::ones(Sparsity::diag(n));
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::inf(const Sparsity& sp) {
    casadi_assert(std::numeric_limits<Scalar>::has_infinity,
                          "Datatype cannot represent infinity");
    return Matrix<Scalar>(sp, std::numeric_limits<Scalar>::infinity(), false);
  }


  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::inf(const std::pair<casadi_int, casadi_int>& rc) {
    return inf(rc.first, rc.second);
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::inf(casadi_int nrow, casadi_int ncol) {
    return inf(Sparsity::dense(nrow, ncol));
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::nan(const Sparsity& sp) {
    casadi_assert(std::numeric_limits<Scalar>::has_quiet_NaN,
                          "Datatype cannot represent not-a-number");
    return Matrix<Scalar>(sp, std::numeric_limits<Scalar>::quiet_NaN(), false);
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::nan(const std::pair<casadi_int, casadi_int>& rc) {
    return nan(rc.first, rc.second);
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::nan(casadi_int nrow, casadi_int ncol) {
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
  casadi_int Matrix<Scalar>::element_hash() const {
    casadi_error("'element_hash' not defined for " + type_name());
  }

  template<typename Scalar>
  bool Matrix<Scalar>::is_leaf() const {
    casadi_error("'is_leaf' not defined for " + type_name());
  }

  template<typename Scalar>
  bool Matrix<Scalar>::is_commutative() const {
    casadi_error("'is_commutative' not defined for " + type_name());
  }

  template<typename Scalar>
  bool Matrix<Scalar>::is_symbolic() const {
    return false;
  }

  template<typename Scalar>
  casadi_int Matrix<Scalar>::op() const {
    casadi_error("'op' not defined for " + type_name());
  }

  template<typename Scalar>
  bool Matrix<Scalar>::is_op(casadi_int k) const {
    casadi_error("'is_op' not defined for " + type_name());
  }

  template<typename Scalar>
  void Matrix<Scalar>::export_code(const std::string& lang,
       std::ostream &stream, const Dict& options) const {
    casadi_error("'export_code' not defined for " + type_name());
  }

  template<typename Scalar>
  bool Matrix<Scalar>::is_valid_input() const {
    return false;
  }

  template<typename Scalar>
  bool Matrix<Scalar>::has_duplicates() const {
    casadi_error("'has_duplicates' not defined for " + type_name());
  }

  template<typename Scalar>
  void Matrix<Scalar>::reset_input() const {
    casadi_error("'reset_input' not defined for " + type_name());
  }

  template<typename Scalar>
  Matrix<double> Matrix<Scalar>::from_file(const std::string& filename,
      const std::string& format_hint) {
    casadi_error("'from_file' not defined for " + type_name());
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
  bool Matrix<Scalar>::is_eye() const {

    // Make sure that the matrix is diagonal
    if (!sparsity().is_diag()) return false;

    // Make sure that all entries are one
    for (auto&& e : nonzeros()) if (!casadi_limits<Scalar>::is_one(e)) return false;

    return true;
  }

  template<typename Scalar>
  bool Matrix<Scalar>::is_equal(const Matrix<Scalar> &x, const Matrix<Scalar> &y,
      casadi_int depth) {
    // Assert matching dimensions
    casadi_assert(x.size() == y.size(), "Dimension mismatch");

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

  // To avoid overloaded function name conflicts
  template<typename Scalar>
  inline Matrix<Scalar> mmin_nonstatic(const Matrix<Scalar> &x) {
    if (x.is_empty()) return Matrix<Scalar>();
    return casadi_mmin(x.ptr(), x.nnz(), x.is_dense());
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::mmin(const Matrix<Scalar> &x) {
    return mmin_nonstatic(x);
  }

  // To avoid overloaded function name conflicts
  template<typename Scalar>
  inline Matrix<Scalar> mmax_nonstatic(const Matrix<Scalar> &x) {
    if (x.is_empty()) return Matrix<Scalar>();
    return casadi_mmax(x.ptr(), x.nnz(), x.is_dense());
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::mmax(const Matrix<Scalar> &x) {
    return mmax_nonstatic(x);
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
    casadi_error("'name' not defined for " + type_name());
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::dep(casadi_int ch) const {
    casadi_error("'dep' not defined for " + type_name());
  }

  template<typename Scalar>
  casadi_int Matrix<Scalar>::n_dep() const {
    casadi_error("'n_dep' not defined for " + type_name());
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::project(const Matrix<Scalar>& x,
                                         const Sparsity& sp, bool intersect) {
    if (intersect) {
      return project(x, sp.intersect(x.sparsity()), false);
    } else {
      casadi_assert(sp.size()==x.size(), "Dimension mismatch");
      Matrix<Scalar> ret = Matrix<Scalar>::zeros(sp);
      std::vector<Scalar> w(x.size1());
      casadi_project(x.ptr(), x.sparsity(), ret.ptr(), sp, get_ptr(w));
      return ret;
    }
  }

  template<typename Scalar>
  void Matrix<Scalar>::set_max_depth(casadi_int eq_depth) {
    casadi_error("'set_max_depth' not defined for " + type_name());
  }

  template<typename Scalar>
  casadi_int Matrix<Scalar>::get_max_depth() {
    casadi_error("'get_max_depth' not defined for " + type_name());
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::det(const Matrix<Scalar>& x) {
    casadi_int n = x.size2();
    casadi_assert(n == x.size1(), "matrix must be square");

    // Trivial return if scalar
    if (x.is_scalar()) return x;

    // Trivial case 2 x 2
    if (n==2) return x(0, 0) * x(1, 1) - x(0, 1) * x(1, 0);

    // Return expression
    Matrix<Scalar> ret = 0;

    // Find out which is the best direction to expand along

    // Build up an IM with ones on the non-zeros
    Matrix<casadi_int> sp = IM::ones(x.sparsity());

    // Have a count of the nonzeros for each row
    Matrix<casadi_int> row_count = Matrix<casadi_int>::sum2(sp);

    // A blank row? determinant is structurally zero
    if (!row_count.is_dense()) return 0;

    // Have a count of the nonzeros for each col
    Matrix<casadi_int> col_count = Matrix<casadi_int>::sum1(sp).T();

    // A blank col? determinant is structurally zero
    if (!row_count.is_dense()) return 0;

    casadi_int min_row = std::distance(row_count.nonzeros().begin(),
                                std::min_element(row_count.nonzeros().begin(),
                                                 row_count.nonzeros().end()));
    casadi_int min_col = std::distance(col_count.nonzeros().begin(),
                                std::min_element(col_count.nonzeros().begin(),
                                                 col_count.nonzeros().end()));

    if (min_row <= min_col) {
      // Expand along row j
      casadi_int j = row_count.sparsity().row(min_row);

      Matrix<Scalar> row = x(j, Slice(0, n));

      std::vector< casadi_int > col_i = row.sparsity().get_col();

      for (casadi_int k=0; k<row.nnz(); ++k) {
        // Sum up the cofactors
        ret += row->at(k)*cofactor(x, col_i.at(k), j);
      }
      return ret;
    } else {
      // Expand along col i
      casadi_int i = col_count.sparsity().row(min_col);

      Matrix<Scalar> col = x(Slice(0, n), i);

      const casadi_int* row_i = col.row();

      for (casadi_int k=0; k<col.nnz(); ++k) {
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
  Matrix<Scalar> Matrix<Scalar>::minor(const Matrix<Scalar>& x,
                                              casadi_int i, casadi_int j) {
    casadi_int n = x.size2();
    casadi_assert(n == x.size1(), "minor: matrix must be square");

    // Trivial return if scalar
    if (n==1) return 1;

    // Remove col i and row j
    Matrix<Scalar> M = Matrix<Scalar>(n-1, n-1);

    std::vector<casadi_int> col = x.sparsity().get_col();
    const casadi_int* row = x.sparsity().row();

    for (casadi_int k=0; k<x.nnz(); ++k) {
      casadi_int i1 = col[k];
      casadi_int j1 = row[k];

      if (i1 == i || j1 == j) continue;

      casadi_int i2 = (i1<i)?i1:i1-1;
      casadi_int j2 = (j1<j)?j1:j1-1;

      M(j2, i2) = x(j1, i1);
    }
    return det(M);
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::cofactor(const Matrix<Scalar>& A, casadi_int i, casadi_int j) {

    // Calculate the i, j minor
    Matrix<Scalar> minor_ij = minor(A, i, j);
    // Calculate the cofactor
    casadi_int sign_i = 1-2*((i+j) % 2);

    return sign_i * minor_ij;
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::adj(const Matrix<Scalar>& x) {
    casadi_int n = x.size2();
    casadi_assert(n == x.size1(), "adj: matrix must be square");

    // Temporary placeholder
    Matrix<Scalar> temp;

    // Cofactor matrix
    Matrix<Scalar> C = Matrix<Scalar>(n, n);
    for (casadi_int i=0; i<n; ++i)
      for (casadi_int j=0; j<n; ++j) {
        temp = cofactor(x, i, j);
        if (!temp.is_zero()) C(j, i) = temp;
      }

    return C.T();
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::inv_minor(const Matrix<Scalar>& x) {
    // laplace formula
    return adj(x)/det(x);
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::reshape(const Matrix<Scalar>& x,
      casadi_int nrow, casadi_int ncol) {
    Sparsity sp = Sparsity::reshape(x.sparsity(), nrow, ncol);
    return Matrix<Scalar>(sp, x.nonzeros(), false);
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::reshape(const Matrix<Scalar>& x, const Sparsity& sp) {
    // quick return if already the right shape
    if (sp==x.sparsity()) return x;

    // make sure that the patterns match
    casadi_assert_dev(sp.is_reshape(x.sparsity()));

    return Matrix<Scalar>(sp, x.nonzeros(), false);
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::trace(const Matrix<Scalar>& x) {
    casadi_assert(x.is_square(), "trace: must be square");
    Scalar res=0;
    const Scalar* d=x.ptr();
    casadi_int size2 = x.size2();
    const casadi_int *colind=x.colind(), *row=x.row();
    for (casadi_int c=0; c<size2; c++) {
      for (casadi_int k=colind[c]; k!=colind[c+1]; ++k) {
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
    for (casadi_int i=0; i<v.size(); ++i)
      ret.push_back(horzcat(v[i]));
    return vertcat(ret);
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::horzcat(const std::vector<Matrix<Scalar> > &v) {
    // Concatenate sparsity patterns
    std::vector<Sparsity> sp(v.size());
    for (casadi_int i=0; i<v.size(); ++i) sp[i] = v[i].sparsity();
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
  Matrix<Scalar>::horzsplit(const Matrix<Scalar>& x, const std::vector<casadi_int>& offset) {
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
    casadi_assert_dev(i==x.nonzeros().end());
    return ret;
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::vertcat(const std::vector<Matrix<Scalar> > &v) {
    std::vector<Matrix<Scalar> > vT(v.size());
    for (casadi_int i=0; i<v.size(); ++i) vT[i] = v[i].T();
    return horzcat(vT).T();
  }

  template<typename Scalar>
  std::vector< Matrix<Scalar> >
  Matrix<Scalar>::vertsplit(const Matrix<Scalar>& x, const std::vector<casadi_int>& offset) {
    std::vector< Matrix<Scalar> > ret = horzsplit(x.T(), offset);
    for (auto&& e : ret) e = e.T();
    return ret;
  }

  template<typename Scalar>
  std::vector< Matrix<Scalar> >
  Matrix<Scalar>::diagsplit(const Matrix<Scalar>& x, const std::vector<casadi_int>& offset1,
                              const std::vector<casadi_int>& offset2) {
    // Consistency check
    casadi_assert_dev(!offset1.empty());
    casadi_assert_dev(offset1.front()==0);
    casadi_assert_dev(offset1.back()==x.size1());
    casadi_assert_dev(is_monotone(offset1));

    // Consistency check
    casadi_assert_dev(!offset2.empty());
    casadi_assert_dev(offset2.front()==0);
    casadi_assert_dev(offset2.back()==x.size2());
    casadi_assert_dev(is_monotone(offset2));

    // Number of outputs
    casadi_int n = offset1.size()-1;

    // Return value
    std::vector< Matrix<Scalar> > ret;

    // Caveat: this is a very silly implementation
    for (casadi_int i=0; i<n; ++i) {
      ret.push_back(x(Slice(offset1[i], offset1[i+1]), Slice(offset2[i], offset2[i+1])));
    }

    return ret;
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::dot(const Matrix<Scalar> &x,
                                     const Matrix<Scalar> &y) {
    casadi_assert(x.size()==y.size(), "dot: Dimension mismatch");
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
    for (casadi_int i=0; i<x.nnz(); ++i) {
      ret = ret && x->at(i)==1;
    }
    return ret;
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::any(const Matrix<Scalar>& x) {
    if (!x.is_dense()) return false;
    Scalar ret=0;
    for (casadi_int i=0; i<x.nnz(); ++i) {
      ret = ret || x->at(i)==1;
    }
    return ret;
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::norm_1(const Matrix<Scalar>& x) {
    return casadi_norm_1(x.nnz(), x.ptr());
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::norm_2(const Matrix<Scalar>& x) {
    if (x.is_vector()) {
      return norm_fro(x);
    } else {
      casadi_error("2-norms currently only supported for vectors. "
                   "Did you intend to calculate a Frobenius norms (norm_fro)?");
    }
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::norm_fro(const Matrix<Scalar>& x) {
    return casadi_norm_2(x.nnz(), x.ptr());
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
  void Matrix<Scalar>::
  qr_sparse(const Matrix<Scalar>& A,
    Matrix<Scalar>& V, Matrix<Scalar> &R, Matrix<Scalar>& beta,
    std::vector<casadi_int>& prinv, std::vector<casadi_int>& pc, bool amd) {
    // Calculate the pattern
    Sparsity spV, spR;
    A.sparsity().qr_sparse(spV, spR, prinv, pc, amd);
    // Calculate the nonzeros
    casadi_int nrow_ext = spV.size1(), ncol = spV.size2();
    V = nan(spV);
    R = nan(spR);
    beta = nan(ncol, 1);
    vector<Scalar> w(nrow_ext);
    casadi_qr(A.sparsity(), A.ptr(), get_ptr(w), spV, V.ptr(),
              spR, R.ptr(), beta.ptr(),
              get_ptr(prinv), get_ptr(pc));
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::
  qr_solve(const Matrix<Scalar>& b, const Matrix<Scalar>& v,
           const Matrix<Scalar>& r, const Matrix<Scalar>& beta,
           const std::vector<casadi_int>& prinv, const std::vector<casadi_int>& pc,
           bool tr) {
    // Get dimensions, check consistency
    casadi_int ncol = v.size2();
    casadi_int nrow = b.size1(), nrhs = b.size2();
    casadi_assert(r.size()==v.size(), "'r', 'v' dimension mismatch");
    casadi_assert(beta.is_vector() && beta.numel()==ncol, "'beta' has wrong dimension");
    casadi_assert(prinv.size()==r.size1(), "'pinv' has wrong dimension");
    // Work vector
    std::vector<Scalar> w(nrow+ncol);
    // Return value
    Matrix<Scalar> x = densify(b);
    casadi_qr_solve(x.ptr(), nrhs, tr, v.sparsity(), v.ptr(), r.sparsity(), r.ptr(),
                    beta.ptr(), get_ptr(prinv), get_ptr(pc), get_ptr(w));
    return x;
  }

  template<typename Scalar>
  void Matrix<Scalar>::qr(const Matrix<Scalar>& A,
                            Matrix<Scalar>& Q, Matrix<Scalar> &R) {
    // The following algorithm is taken from J. Demmel:
    // Applied Numerical Linear Algebra (algorithm 3.1.)
    casadi_assert(A.size1()>=A.size2(), "qr: fewer rows than columns");

    // compute Q and R column by column
    Q = R = Matrix<Scalar>();
    for (casadi_int i=0; i<A.size2(); ++i) {
      // Initialize qi to be the i-th column of *this
      Matrix<Scalar> ai = A(Slice(), i);
      Matrix<Scalar> qi = ai;
      // The i-th column of R
      Matrix<Scalar> ri = Matrix<Scalar>(A.size2(), 1);

      // subtract the projection of qi in the previous directions from ai
      for (casadi_int j=0; j<i; ++j) {

        // Get the j-th column of Q
        Matrix<Scalar> qj = Q(Slice(), j); // NOLINT(cppcoreguidelines-slicing)

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
  void Matrix<Scalar>::ldl(const Matrix<Scalar>& A, Matrix<Scalar> &D,
    Matrix<Scalar>& LT, std::vector<casadi_int>& p, bool amd) {
    // Symbolic factorization
    Sparsity Lt_sp = A.sparsity().ldl(p, amd);

    // Get dimension
    casadi_int n=A.size1();

    // Calculate entries in L and D
    vector<Scalar> D_nz(n), L_nz(Lt_sp.nnz()), w(n);
    casadi_ldl(A.sparsity(), get_ptr(A.nonzeros()), Lt_sp,
              get_ptr(L_nz), get_ptr(D_nz), get_ptr(p), get_ptr(w));

    // Assemble L and D
    LT = Matrix<Scalar>(Lt_sp, L_nz);
    D = D_nz;
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::
  ldl_solve(const Matrix<Scalar>& b, const Matrix<Scalar>& D, const Matrix<Scalar>& LT,
            const std::vector<casadi_int>& p) {
    // Get dimensions, check consistency
    casadi_int n = b.size1(), nrhs = b.size2();
    casadi_assert(p.size()==n, "'p' has wrong dimension");
    casadi_assert(LT.size1()==n && LT.size2()==n, "'LT' has wrong dimension");
    casadi_assert(D.is_vector() && D.numel()==n, "'D' has wrong dimension");
    // Solve for all right-hand-sides
    Matrix<Scalar> x = densify(b);
    std::vector<Scalar> w(n);
    casadi_ldl_solve(x.ptr(), nrhs, LT.sparsity(), LT.ptr(), D.ptr(), get_ptr(p), get_ptr(w));
    return x;
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::nullspace(const Matrix<Scalar>& A) {
    Matrix<Scalar> X = A;
    casadi_int n = X.size1();
    casadi_int m = X.size2();
    casadi_assert(m>=n, "nullspace(): expecting a flat matrix (more columns than rows), "
                          "but got " + str(X.dim()) + ".");

    Matrix<Scalar> seed = DM::eye(m)(Slice(0, m), Slice(n, m)); // NOLINT(cppcoreguidelines-slicing)

    std::vector< Matrix<Scalar> > us;
    std::vector< Matrix<Scalar> > betas;

    Matrix<Scalar> beta;

    for (casadi_int i=0;i<n;++i) {
      Matrix<Scalar> x = X(i, Slice(i, m)); // NOLINT(cppcoreguidelines-slicing)
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

    for (casadi_int i=n-1;i>=0;--i) {
      seed(Slice(i, m), Slice(0, m-n)) -=
        betas[i]*mtimes(us[i].T(), mtimes(us[i], seed(Slice(i, m), Slice(0, m-n))));
    }

    return seed;

  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::chol(const Matrix<Scalar>& A) {
    // Perform an LDL transformation
    Matrix<Scalar> D, LT;
    std::vector<casadi_int> p;
    ldl(A, D, LT, p, false);
    // Add unit diagonal
    LT += Matrix<Scalar>::eye(D.size1());
    // Get the cholesky factor: R*R' = L*D*L' = (sqrt(D)*L')'*(sqrt(D)*L')
    return mtimes(diag(sqrt(D)), LT);
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::solve(const Matrix<Scalar>& a, const Matrix<Scalar>& b) {
    // check dimensions
    casadi_assert(a.size1() == b.size1(), "solve Ax=b: dimension mismatch: b has "
                          + str(b.size1()) + " rows while A has " + str(a.size1()) + ".");
    casadi_assert(a.size1() == a.size2(), "solve: A not square but " + str(a.dim()));

    if (a.is_tril()) {
      // forward substitution if lower triangular
      Matrix<Scalar> x = b;
      const casadi_int*  Arow = a.row();
      const casadi_int*  Acolind = a.colind();
      const std::vector<Scalar> & Adata = a.nonzeros();
      for (casadi_int i=0; i<a.size2(); ++i) { // loop over columns forwards
        for (casadi_int k=0; k<b.size2(); ++k) { // for every right hand side
          if (!x.has_nz(i, k)) continue;
          x(i, k) /= a(i, i);
          for (casadi_int kk=Acolind[i+1]-1; kk>=Acolind[i] && Arow[kk]>i; --kk) {
            casadi_int j = Arow[kk];
            x(j, k) -= Adata[kk]*x(i, k);
          }
        }
      }
      return x;
    } else if (a.is_triu()) {
      // backward substitution if upper triangular
      Matrix<Scalar> x = b;
      const casadi_int*  Arow = a.row();
      const casadi_int*  Acolind = a.colind();
      const std::vector<Scalar> & Adata = a.nonzeros();
      for (casadi_int i=a.size2()-1; i>=0; --i) { // loop over columns backwards
        for (casadi_int k=0; k<b.size2(); ++k) { // for every right hand side
          if (!x.has_nz(i, k)) continue;
          x(i, k) /= a(i, i);
          for (casadi_int kk=Acolind[i]; kk<Acolind[i+1] && Arow[kk]<i; ++kk) {
            casadi_int j = Arow[kk];
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
      std::vector<casadi_int> rowperm, colperm, rowblock, colblock;
      std::vector<casadi_int> coarse_rowblock, coarse_colblock;
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
        xperm = mtimes(inv_minor(Aperm), bperm);

      } else {

        // Make a QR factorization
        Matrix<Scalar> Q, R;
        qr(Aperm, Q, R);

        // Solve the factorized system (note that solve will now be fast since it is triangular)
        xperm = solve(R, mtimes(Q.T(), bperm));
      }

      // get the inverted column permutation
      std::vector<casadi_int> inv_colperm(colperm.size());
      for (casadi_int k=0; k<colperm.size(); ++k)
        inv_colperm[colperm[k]] = k;

      // Permute back the solution and return
      Matrix<Scalar> x = xperm(inv_colperm, Slice()); // NOLINT(cppcoreguidelines-slicing)
      return x;
    }
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::
  solve(const Matrix<Scalar>& a, const Matrix<Scalar>& b,
           const std::string& lsolver, const Dict& dict) {
    casadi_error("'solve' with plugin not defined for " + type_name());
    return Matrix<Scalar>();
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::
  inv(const Matrix<Scalar>& a) {
    return solve(a, Matrix<Scalar>::eye(a.size1()));
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::
  inv(const Matrix<Scalar>& a,
           const std::string& lsolver, const Dict& dict) {
    casadi_error("'inv' with plugin not defined for " + type_name());
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
    casadi_error("'solve' not defined for " + type_name());
    return Matrix<Scalar>();
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::
  expm_const(const Matrix<Scalar>& A, const Matrix<Scalar>& t) {
    casadi_error("'solve' not defined for " + type_name());
    return Matrix<Scalar>();
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::
  expm(const Matrix<Scalar>& A) {
    casadi_error("'solve' not defined for " + type_name());
    return Matrix<Scalar>();
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::kron(const Matrix<Scalar>& a, const Matrix<Scalar>& b) {
    std::vector<Scalar> ret(a.nnz()*b.nnz());
    casadi_kron(get_ptr(a), a.sparsity(), get_ptr(b), b.sparsity(), get_ptr(ret));

    Sparsity sp_ret = Sparsity::kron(a.sparsity(), b.sparsity());
    return Matrix<Scalar>(sp_ret, ret, false);
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::diag(const Matrix<Scalar>& A) {
    // Nonzero mapping
    std::vector<casadi_int> mapping;
    // Get the sparsity
    Sparsity sp = A.sparsity().get_diag(mapping);

    Matrix<Scalar> ret = zeros(sp);

    for (casadi_int k=0; k<mapping.size(); k++) ret.nz(k) = A.nz(mapping[k]);
    return ret;
  }

  /** \brief   Construct a matrix with given block on the diagonal */
  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::diagcat(const std::vector< Matrix<Scalar> > &A) {
    std::vector<Scalar> data;

    std::vector<Sparsity> sp;
    for (casadi_int i=0;i<A.size();++i) {
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
    casadi_int elA=0, elB=0;
    for (casadi_int k=0; k<mapping.size(); ++k) {
      if (mapping[k]==1) {
        ret.nonzeros()[k] = A.nonzeros()[elA++];
      } else if (mapping[k]==2) {
        ret.nonzeros()[k] = B.nonzeros()[elB++];
      } else {
        casadi_error("Pattern intersection not empty");
      }
    }

    casadi_assert_dev(A.nnz()==elA);
    casadi_assert_dev(B.nnz()==elB);

    return ret;
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::polyval(const Matrix<Scalar>& p, const Matrix<Scalar>& x) {
    casadi_assert(p.is_dense(), "polynomial coefficients vector must be dense");
    casadi_assert(p.is_vector() && p.nnz()>0, "polynomial coefficients must be a vector");
    Matrix<Scalar> ret = x;
    for (auto&& e : ret.nonzeros()) {
      e = casadi_polyval(p.ptr(), p.numel()-1, e);
    }
    return ret;
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::norm_inf_mul(const Matrix<Scalar>& x,
                                                  const Matrix<Scalar>& y) {
    casadi_assert(y.size1()==x.size2(), "Dimension error. Got " + x.dim()
                          + " times " + y.dim() + ".");

    // Allocate work vectors
    std::vector<Scalar> dwork(x.size1());
    std::vector<casadi_int> iwork(x.size1()+1+y.size2());

    // Call C runtime
    return casadi_norm_inf_mul(x.ptr(), x.sparsity(), y.ptr(), y.sparsity(),
                               get_ptr(dwork), get_ptr(iwork));
  }

  template<typename Scalar>
  void Matrix<Scalar>::expand(const Matrix<Scalar>& ex,
                                Matrix<Scalar> &weights, Matrix<Scalar>& terms) {
    casadi_error("'expand' not defined for " + type_name());
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::pw_const(const Matrix<Scalar>& ex,
                                              const Matrix<Scalar>& tval,
                                              const Matrix<Scalar>& val) {
    casadi_error("'pw_const' not defined for " + type_name());
    return Matrix<Scalar>();
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::pw_lin(const Matrix<Scalar>& ex,
                                            const Matrix<Scalar>& tval,
                                            const Matrix<Scalar>& val) {
    casadi_error("'pw_lin' not defined for " + type_name());
    return Matrix<Scalar>();
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::if_else(const Matrix<Scalar> &cond,
                                             const Matrix<Scalar> &if_true,
                                             const Matrix<Scalar> &if_false,
                                             bool short_circuit) {
    return if_else_zero(cond, if_true) + if_else_zero(!cond, if_false);
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::conditional(const Matrix<Scalar>& ind,
                                                 const std::vector<Matrix<Scalar> >& x,
                                                 const Matrix<Scalar>& x_default,
                                                 bool short_circuit) {
    casadi_assert(!short_circuit,
      "Short-circuiting 'conditional' not supported for " + type_name());
    casadi_assert(ind.is_scalar(true),
      "conditional: first argument must be scalar. Got " + ind.dim()+ " instead.");

    Matrix<Scalar> ret = x_default;
    for (casadi_int k=0; k<x.size(); ++k) {
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
                   const Matrix<Scalar> &b, casadi_int order) {
    return gauss_quadrature(f, x, a, b, order, Matrix<Scalar>());
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::gauss_quadrature(const Matrix<Scalar>& f,
                                                      const Matrix<Scalar>& x,
                                                      const Matrix<Scalar>& a,
                                                      const Matrix<Scalar>& b, casadi_int order,
                                                      const Matrix<Scalar>& w) {
    casadi_error("'gauss_quadrature' not defined for " + type_name());
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
    casadi_error("'substitute' not defined for " + type_name());
    return Matrix<Scalar>();
  }

  template<typename Scalar>
  std::vector<Matrix<Scalar> >
  Matrix<Scalar>::substitute(const std::vector<Matrix<Scalar> >& ex,
                               const std::vector<Matrix<Scalar> >& v,
                               const std::vector<Matrix<Scalar> >& vdef) {
    casadi_error("'substitute' not defined for " + type_name());
    return std::vector<Matrix<Scalar> >();
  }

  template<typename Scalar>
  void Matrix<Scalar>::substitute_inplace(const std::vector<Matrix<Scalar> >& v,
                                           std::vector<Matrix<Scalar> >& vdef,
                                           std::vector<Matrix<Scalar> >& ex,
                                           bool reverse) {
    casadi_error("'substitute_inplace' not defined for " + type_name());
  }

  template<typename Scalar>
  bool Matrix<Scalar>::depends_on(const Matrix<Scalar> &x, const Matrix<Scalar> &arg) {
    casadi_error("'depends_on' not defined for " + type_name());
    return false;
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::
  jacobian(const Matrix<Scalar> &f, const Matrix<Scalar> &x, const Dict& opts) {
    casadi_error("'jacobian' not defined for " + type_name());
    return Matrix<Scalar>();
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::hessian(const Matrix<Scalar> &f,
                                             const Matrix<Scalar> &x,
                                             const Dict& opts) {
    casadi_error("'hessian' not defined for " + type_name());
    return Matrix<Scalar>();
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::hessian(const Matrix<Scalar> &f,
                                             const Matrix<Scalar> &x,
                                             Matrix<Scalar> &g,
                                              const Dict& opts) {
    casadi_error("'hessian' not defined for " + type_name());
    return Matrix<Scalar>();
  }

  template<typename Scalar>
  std::vector<std::vector<Matrix<Scalar> > >
  Matrix<Scalar>::
  forward(const std::vector<Matrix<Scalar> > &ex,
          const std::vector<Matrix<Scalar> > &arg,
          const std::vector<std::vector<Matrix<Scalar> > > &v,
          const Dict& opts) {
    casadi_error("'forward' not defined for " + type_name());
  }

  template<typename Scalar>
  std::vector<std::vector<Matrix<Scalar> > >
  Matrix<Scalar>::
  reverse(const std::vector<Matrix<Scalar> > &ex,
          const std::vector<Matrix<Scalar> > &arg,
          const std::vector<std::vector<Matrix<Scalar> > > &v,
          const Dict& opts) {
    casadi_error("'reverse' not defined for " + type_name());
  }

  template<typename Scalar>
  std::vector<bool>
  Matrix<Scalar>::which_depends(const Matrix<Scalar> &expr, const Matrix<Scalar> &var,
      casadi_int order, bool tr) {
    casadi_error("'which_depends' not defined for " + type_name());
    return std::vector<bool>();
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::taylor(const Matrix<Scalar>& f,
                                            const Matrix<Scalar>& x,
                                            const Matrix<Scalar>& a, casadi_int order) {
    casadi_error("'taylor' not defined for " + type_name());
    return Matrix<Scalar>();
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::mtaylor(const Matrix<Scalar>& f,
                                             const Matrix<Scalar>& x,
                                             const Matrix<Scalar>& a, casadi_int order) {
    casadi_error("'mtaylor' not defined for " + type_name());
    return Matrix<Scalar>();
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::mtaylor(const Matrix<Scalar>& f,
                                             const Matrix<Scalar>& x,
                                             const Matrix<Scalar>& a, casadi_int order,
                                             const std::vector<casadi_int>&order_contributions) {
    casadi_error("'mtaylor' not defined for " + type_name());
    return Matrix<Scalar>();
  }

  template<typename Scalar>
  casadi_int Matrix<Scalar>::n_nodes(const Matrix<Scalar>& x) {
    casadi_error("'n_nodes' not defined for " + type_name());
    return 0;
  }

  template<typename Scalar>
  std::string
  Matrix<Scalar>::print_operator(const Matrix<Scalar>& x,
                                   const std::vector<std::string>& args) {
    casadi_error("'print_operator' not defined for " + type_name());
    return std::string();
  }

  template<typename Scalar>
  std::vector<Matrix<Scalar> > Matrix<Scalar>::symvar(const Matrix<Scalar>& x) {
    casadi_error("'symvar' not defined for " + type_name());
    return std::vector<Matrix<Scalar> >();
  }

  template<typename Scalar>
  void Matrix<Scalar>::shared(std::vector<Matrix<Scalar> >& ex,
                                       std::vector<Matrix<Scalar> >& v,
                                       std::vector<Matrix<Scalar> >& vdef,
                                       const std::string& v_prefix,
                                       const std::string& v_suffix) {
    casadi_error("'shared' not defined for " + type_name());
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::poly_coeff(const Matrix<Scalar>& f,
                                                const Matrix<Scalar>&x) {
    casadi_error("'poly_coeff' not defined for " + type_name());
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::poly_roots(const Matrix<Scalar>& p) {
    casadi_error("'poly_roots' not defined for " + type_name());
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::eig_symbolic(const Matrix<Scalar>& m) {
    casadi_error("'eig_symbolic' not defined for " + type_name());
  }

  template<typename Scalar>
  DM Matrix<Scalar>::evalf(const Matrix<Scalar>& m) {
    Function f("f", std::vector<SX>{}, std::vector<SX>{m});
    return f(std::vector<DM>{})[0];
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::sparsify(const Matrix<Scalar>& x, double tol) {
    // Quick return if there are no entries to be removed
    bool remove_nothing = true;
    for (auto it=x.nonzeros().begin(); it!=x.nonzeros().end() && remove_nothing; ++it) {
      remove_nothing = !casadi_limits<Scalar>::is_almost_zero(*it, tol);
    }
    if (remove_nothing) return x;

    // Get the current sparsity pattern
    casadi_int size1 = x.size1();
    casadi_int size2 = x.size2();
    const casadi_int* colind = x.colind();
    const casadi_int* row = x.row();

    // Construct the new sparsity pattern
    std::vector<casadi_int> new_colind(1, 0), new_row;
    std::vector<Scalar> new_data;

    // Loop over the columns
    for (casadi_int cc=0; cc<size2; ++cc) {
      // Loop over existing nonzeros
      for (casadi_int el=colind[cc]; el<colind[cc+1]; ++el) {
        // If it is not known to be a zero
        if (!casadi_limits<Scalar>::is_almost_zero(x->at(el), tol)) {
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
    casadi_error("'get_input' not defined for " + type_name());
  }

  template<typename Scalar>
  std::vector<Matrix<Scalar> > Matrix<Scalar>::get_free(const Function& f) {
    casadi_error("'get_free' not defined for " + type_name());
  }

  template<typename Scalar>
  Matrix<Scalar>::operator double() const {
    casadi_assert_dev(is_scalar());
    return static_cast<double>(scalar());
  }

  template<typename Scalar>
  Matrix<Scalar>::operator casadi_int() const {
    casadi_assert_dev(is_scalar());
    return static_cast<casadi_int>(scalar());
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::_sym(const std::string& name, const Sparsity& sp) {
    casadi_error("'sym' not defined for " + type_name());
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::rand(const Sparsity& sp) { // NOLINT(runtime/threadsafe_fn)

    casadi_error("'rand' not defined for " + type_name());
  }

  template<typename Scalar>
  std::string Matrix<Scalar>::serialize() const {
    std::stringstream ss;
    serialize(ss);
    return ss.str();
  }

  template<typename Scalar>
  void Matrix<Scalar>::serialize(SerializingStream& s) const {
    s.pack("Matrix::sparsity", sparsity());
    s.pack("Matrix::nonzeros", nonzeros());
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::deserialize(DeserializingStream& s) {
    Sparsity sp;
    s.unpack("Matrix::sparsity", sp);
    std::vector<Scalar> nz;
    s.unpack("Matrix::nonzeros", nz);
    return Matrix<Scalar>(sp, nz, false);
  }

  template<typename Scalar>
  void Matrix<Scalar>::serialize(std::ostream &stream) const {
    SerializingStream s(stream);
    serialize(s);
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::deserialize(std::istream &stream) {
    DeserializingStream s(stream);
    return Matrix<Scalar>::deserialize(s);
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::deserialize(const std::string& s) {
    std::stringstream ss;
    ss << s;
    return deserialize(ss);
  }


} // namespace casadi

#endif // CASADI_MATRIX_IMPL_HPP
