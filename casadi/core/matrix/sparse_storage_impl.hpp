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
  DataType& SparseStorage<DataType>::elem(int rr, int cc) {
    int oldsize = sparsity().nnz();
    int ind = sparsityRef().addNZ(rr, cc);
    if (oldsize != sparsity().nnz())
      data().insert(data().begin()+ind, DataType(0));
    return data().at(ind);
  }

  template<typename DataType>
  SparseStorage<DataType>::SparseStorage(const Sparsity& sparsity, const DataType& val) :
    sparsity_(sparsity), data_(sparsity.nnz(), val) {}

  template<typename DataType>
  SparseStorage<DataType>::SparseStorage() : sparsity_(Sparsity(0, 0)) {
  }

  template<typename DataType>
  SparseStorage<DataType>::SparseStorage(const SparseStorage<DataType>& m) :
    sparsity_(m.sparsity_), data_(m.data_) {}

  template<typename DataType>
  SparseStorage<DataType>& SparseStorage<DataType>::operator=(const SparseStorage<DataType>& m) {
    sparsity_ = m.sparsity_;
    data_ = m.data_;
    return *this;
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
    sparsity_ = Sparsity(0, 0);
    data().clear();
  }

  template<typename DataType>
  std::vector<DataType>& SparseStorage<DataType>::data() {
    return data_;
  }

  template<typename DataType>
  const std::vector<DataType>& SparseStorage<DataType>::data() const {
    return data_;
  }

  template<typename DataType>
  Sparsity& SparseStorage<DataType>::sparsityRef() {
    sparsity_.makeUnique(); // NOTE: Remove?
    return sparsity_;
  }

} // namespace casadi

/// \endcond

#endif // CASADI_SPARSE_STORAGE_IMPL_HPP
