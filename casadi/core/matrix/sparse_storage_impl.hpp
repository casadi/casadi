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


#ifndef CASADI_SPARSE_STORAGE_IMPL_HPP
#define CASADI_SPARSE_STORAGE_IMPL_HPP

/// \cond INTERNAL

namespace casadi {
  // Implementations

  template<typename DataType>
  const DataType& SparseStorage<DataType>::elem(int rr, int cc) const {
    int ind = sparsity().getNZ(rr, cc);
    if (ind==-1)
      return casadi_limits<DataType>::zero;
    else
      return at(ind);
  }

  template<typename DataType>
  DataType& SparseStorage<DataType>::elem(int rr, int cc) {
    int oldsize = sparsity().size();
    int ind = sparsityRef().getNZ(rr, cc);
    if (oldsize != sparsity().size())
      data().insert(begin()+ind, DataType(0));
    return at(ind);
  }

  template<typename DataType>
  SparseStorage<DataType>::SparseStorage() : sparsity_(Sparsity::sparse(0, 0)) {
  }

  template<typename DataType>
  SparseStorage<DataType>::SparseStorage(const SparseStorage<DataType>& m) :
    sparsity_(m.sparsity_), data_(m.data_) {}

  template<typename DataType>
  SparseStorage<DataType>::SparseStorage(const std::vector<DataType>& x) :
    sparsity_(Sparsity::dense(x.size(), 1)), data_(x) {}

  template<typename DataType>
  SparseStorage<DataType>::SparseStorage(const std::vector<DataType>& x, int nrow, int ncol)
      : sparsity_(Sparsity::dense(nrow, ncol)), data_(x) {
    casadi_assert_message(x.size() == nrow*ncol,
                          "Dimension mismatch." << std::endl
                          << "You supplied a vector of length " << x.size()
                          << ", but " << nrow << " x " << ncol << " = " << nrow*ncol);
  }

  template<typename DataType>
  SparseStorage<DataType>& SparseStorage<DataType>::operator=(const SparseStorage<DataType>& m) {
    sparsity_ = m.sparsity_;
    data_ = m.data_;
    return *this;
  }

  template<typename DataType>
  const std::vector<int>& SparseStorage<DataType>::row() const {
    return sparsity().row();
  }

  template<typename DataType>
  const std::vector<int>& SparseStorage<DataType>::colind() const {
    return sparsity_.colind();
  }

  template<typename DataType>
  int SparseStorage<DataType>::row(int el) const {
    return sparsity_.row(el);
  }

  template<typename DataType>
  int SparseStorage<DataType>::colind(int col) const {
    return sparsity_.colind(col);
  }

  template<typename DataType>
  void SparseStorage<DataType>::reserve(int nnz) {
    reserve(nnz, sparsity_.size2());
  }

  template<typename DataType>
  void SparseStorage<DataType>::reserve(int nnz, int ncol) {
    data().reserve(nnz);
    sparsity_.reserve(nnz, ncol);
  }

  template<typename DataType>
  void SparseStorage<DataType>::resize(int nrow, int ncol) {
    sparsity_.resize(nrow, ncol);
  }

  template<typename DataType>
  void SparseStorage<DataType>::clear() {
    sparsity_ = Sparsity::sparse(0, 0);
    data().clear();
  }

  template<typename DataType>
  SparseStorage<DataType>::SparseStorage(const std::vector< std::vector<DataType> >& d) {
    // Get dimensions
    int nrow=d.size();
    int ncol=d.empty() ? 1 : d.front().size();

    // Assert consistency
    for (int rr=0; rr<nrow; ++rr) {
      casadi_assert_message(ncol==d[rr].size(),
        "SparseStorage<DataType>::SparseStorage(const std::vector< std::vector<DataType> >& d): "
        "shape mismatch" << std::endl <<
        "Attempting to construct a matrix from a nested list." << std::endl <<
        "I got convinced that the desired size is ("<< nrow << " x " << ncol << " ), "
        "but now I encounter a vector of size (" <<
        d[rr].size() <<  " )" << std::endl);
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
  SparseStorage<DataType>::SparseStorage(const Sparsity& sparsity, const DataType& val) :
    sparsity_(sparsity), data_(sparsity.size(), val) {}

  template<typename DataType>
  SparseStorage<DataType>::SparseStorage(const Sparsity& sparsity, const std::vector<DataType>& d)
      : sparsity_(sparsity), data_(d) {
    casadi_assert_message(sparsity.size()==d.size(),
                          "Size mismatch." << std::endl << "You supplied a sparsity of "
                          << sparsity_.dimString()
                          << ", but the supplied vector is of length " << d.size());
  }

  template<typename DataType>
  Sparsity& SparseStorage<DataType>::sparsityRef() {
    sparsity_.makeUnique(); // NOTE: Remove?
    return sparsity_;
  }

  template<typename DataType>
  std::vector<DataType>& SparseStorage<DataType>::data() {
    return data_;
  }

  template<typename DataType>
  const std::vector<DataType>& SparseStorage<DataType>::data() const {
    return data_;
  }

} // namespace casadi

/// \endcond

#endif // CASADI_SPARSE_STORAGE_IMPL_HPP
