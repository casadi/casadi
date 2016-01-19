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


#ifndef CASADI_SPARSE_STORAGE_HPP
#define CASADI_SPARSE_STORAGE_HPP

#include <vector>
#include <typeinfo>
#include "exception.hpp"

/// \cond INTERNAL

namespace casadi {

  template<typename DataType>
  class CASADI_EXPORT SparseStorage {
  public:

    /** \brief  constructors */
    /// empty 0-by-0 matrix constructor
    SparseStorage();

    /// Copy constructor
    SparseStorage(const SparseStorage<DataType>& m);

    ///@{
    /// Sparse matrix with a given sparsity
    explicit SparseStorage(const Sparsity& sparsity, const DataType& val=DataType(0));
    ///@}

    /// Assignment (normal)
    SparseStorage<DataType>& operator=(const SparseStorage<DataType>& m);

    /// Convert to scalar type
    const DataType scalar() const;

    /// Scalar type
    typedef DataType ScalarType;

    /// get a reference to an element
    DataType& elem(int rr, int cc);

    /// Returns true if the matrix has a non-zero at location rr, cc
    bool has_nz(int rr, int cc) const { return sparsity().has_nz(rr, cc); }

    // Get the sparsity pattern
    void clear();
    void resize(int nrow, int ncol);
    void reserve(int nnz);
    void reserve(int nnz, int ncol);

    /// Access the non-zero elements
    std::vector<DataType>& nonzeros();

    /// Const access the non-zero elements
    const std::vector<DataType>& nonzeros() const;

    /// Const access the sparsity - reference to data member
    const Sparsity& sparsity() const { return sparsity_; }

    /// Access the sparsity, make a copy if there are multiple references to it
    Sparsity& sparsityRef();

  private:
    /// Sparsity of the matrix in a compressed column storage (CCS) format
    Sparsity sparsity_;

    /// Nonzero elements
    std::vector<DataType> nonzeros_;

  };

} // namespace casadi

/// \endcond

#endif // CASADI_SPARSE_STORAGE_HPP
