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
#include "../casadi_exception.hpp"

/// \cond INTERNAL

namespace casadi {

  template<typename DataType>
  class CASADI_CORE_EXPORT SparseStorage {
  public:

    /** \brief  constructors */
    /// empty 0-by-0 matrix constructor
    SparseStorage();

    /// Copy constructor
    SparseStorage(const SparseStorage<DataType>& m);

#ifndef SWIG
    /// Assignment (normal)
    SparseStorage<DataType>& operator=(const SparseStorage<DataType>& m);
#endif // SWIG

    /// Dense matrix constructor with data given as vector of vectors
    explicit SparseStorage(const std::vector< std::vector<DataType> >& m);

    ///@{
    /// Sparse matrix with a given sparsity
    explicit SparseStorage(const Sparsity& sparsity, const DataType& val=DataType(0));
    ///@}

    /// Sparse matrix with a given sparsity and non-zero elements.
    SparseStorage(const Sparsity& sparsity, const std::vector<DataType>& d);

    /** \brief Check if the dimensions and colind, row vectors are compatible.
     * \param complete  set to true to also check elementwise
     * throws an error as possible result
     */
    void sanityCheck(bool complete=false) const;

    /// Construct from a vector
    /**
     * Thanks to implicit conversion, you can pretend that SparseStorage(const SXElement& x); exists.
     * Note: above remark applies only to C++, not python or octave interfaces
     */
    SparseStorage(const std::vector<DataType>& x);

    /// Construct dense matrix from a vector with the elements in column major ordering
    SparseStorage(const std::vector<DataType>& x, int nrow, int ncol);

    /// Convert to scalar type
    const DataType toScalar() const;

    /// Scalar type
    typedef DataType ScalarType;

    /// \cond INTERNAL
    /// Expose iterators
    typedef typename std::vector<DataType>::iterator iterator;
    typedef typename std::vector<DataType>::const_iterator const_iterator;
    typedef typename std::vector<DataType>::reverse_iterator reverse_iterator;
    typedef typename std::vector<DataType>::const_reverse_iterator const_reverse_iterator;

    /// References
    typedef DataType& reference;
    typedef const DataType& const_reference;

    /// Get iterators to beginning and end
    iterator begin() { return data().begin();}
    const_iterator begin() const { return data().begin();}
    reverse_iterator rbegin() { return data().rbegin();}
    const_reverse_iterator rbegin() const { return data().rbegin();}
    iterator end() { return data().end();}
    const_iterator end() const { return data().end();}
    reverse_iterator rend() { return data().rend();}
    const_reverse_iterator rend() const { return data().rend();}

    /// Get references to beginning and end
    reference front() { return data().front();}
    const_reference front() const { return data().front();}
    reference back() { return data().back();}
    const_reference back() const { return data().back();}
    /// \endcond

    /** \brief  Create a matrix from a matrix with a different type of matrix entries
     *          (assuming that the scalar conversion is valid) */
    template<typename A>
    SparseStorage(const SparseStorage<A>& x) :
        sparsity_(x.sparsity()), data_(std::vector<DataType>(x.size())) {
      copy(x.begin(), x.end(), begin());
    }

    /** \brief  Create an expression from an stl vector  */
    template<typename A>
    SparseStorage(const std::vector<A>& x) :
      sparsity_(Sparsity::dense(x.size(), 1)), data_(std::vector<DataType>(x.size())) {
      copy(x.begin(), x.end(), begin());
    }

    /** \brief  Create a non-vector expression from an stl vector */
    template<typename A>
    SparseStorage(const std::vector<A>& x,  int nrow, int ncol) :
      sparsity_(Sparsity::dense(nrow, ncol)), data_(std::vector<DataType>(x.size())) {
      if (x.size() != nrow*ncol)
        throw CasadiException("SparseStorage::SparseStorage(const std::vector<DataType>& x, "
                              "int n, int m): dimension mismatch");
      copy(x.begin(), x.end(), begin());
    }


#ifndef SWIG
    /// Get a non-zero element
    inline const DataType& at(int k) const {
      return const_cast<SparseStorage<DataType>*>(this)->at(k);
    }

    /// Access a non-zero element
    inline DataType& at(int k) {
      try {
        if (k<0) k+=sparsity_.size();
        return data().at(k);
      } catch(std::out_of_range& ex) {
        std::stringstream ss;
        ss << "Out of range error in SparseStorage<>::at: " << k
           << " not in range [0, " << sparsity_.size() << ")";
        throw CasadiException(ss.str());
      }
    }
#else // SWIG
    /// Access a non-zero element
    DataType at(int k) {
      try {
        if (k<0) k+=sparsity_.size();
        return data().at(k);
      } catch(std::out_of_range& ex) {
        std::stringstream ss;
        ss << "Out of range error in SparseStorage<>::at: " << k
           << " not in range [0, " << size() << ")";
        throw CasadiException(ss.str());
      }
    }
#endif // SWIG

#ifndef SWIG
    /// get an element
    const DataType& elem(int rr, int cc=0) const;

    /// get a reference to an element
    DataType& elem(int rr, int cc=0);
#else // SWIG
    /// Access a non-zero element
    DataType elem(int rr, int cc=0) { return elem(rr, cc);}
#endif // SWIG

    /// get an element, do not allocate
    const DataType getElement(int rr, int cc=0) const { return elem(rr, cc);}

    /// Returns true if the matrix has a non-zero at location rr, cc
    bool hasNZ(int rr, int cc) const { return sparsity().hasNZ(rr, cc); }

    // Get the sparsity pattern
    const std::vector<int>& row() const;
    const std::vector<int>& colind() const;
    int row(int el) const;
    int colind(int col) const;
    void clear();
    void resize(int nrow, int ncol);
    void reserve(int nnz);
    void reserve(int nnz, int ncol);

    /// Access the non-zero elements
    std::vector<DataType>& data();

    /// Const access the non-zero elements
    const std::vector<DataType>& data() const;

    /// Const access the sparsity - reference to data member
    const Sparsity& sparsity() const { return sparsity_; }

    /// Access the sparsity, make a copy if there are multiple references to it
    Sparsity& sparsityRef();

  private:
    /// Sparsity of the matrix in a compressed column storage (CCS) format
    Sparsity sparsity_;

    /// Nonzero elements
    std::vector<DataType> data_;

  };

} // namespace casadi

/// \endcond

#endif // CASADI_SPARSE_STORAGE_HPP
