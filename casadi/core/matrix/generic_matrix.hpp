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


#ifndef CASADI_GENERIC_MATRIX_HPP
#define CASADI_GENERIC_MATRIX_HPP

#include "slice.hpp"
#include "submatrix.hpp"
#include "nonzeros.hpp"
#include "sparsity.hpp"
#include "sparsity_tools.hpp"
#include "../casadi_math.hpp"

namespace casadi {

  /** Sparsity format for getting and setting inputs and outputs */
  enum SparsityType {SPARSE, SPARSESYM, DENSE, DENSESYM, DENSETRANS};

  /** \brief Matrix base class

      This is a common base class for MX and Matrix<>, introducing a uniform syntax and implementing
      common functionality using the curiously recurring template pattern (CRTP) idiom.\n

      The class is designed with the idea that "everything is a matrix",
      that is, also scalars and vectors.\n
      This philosophy makes it easy to use and to interface in particularly
      with Python and Matlab/Octave.\n

      The syntax tries to stay as close as possible to the ublas syntax
      when it comes to vector/matrix operations.\n

      Index starts with 0.\n
      Index vec happens as follows: (rr, cc) -> k = rr+cc*size1()\n
      Vectors are column vectors.\n

      The storage format is Compressed Column Storage (CCS),
      similar to that used for sparse matrices in Matlab, \n
      but unlike this format, we do allow for elements to be structurally non-zero
      but numerically zero.\n

      The sparsity pattern, which is reference counted and cached,
      can be accessed with Sparsity& sparsity()\n

      \author Joel Andersson
      \date 2012
  */
  template<typename MatType>
  class CASADI_CORE_EXPORT GenericMatrix {
  public:

    /** \brief Get the number of (structural) non-zero elements */
    int size() const;

    /** \brief Get the number of non-zeros in the lower triangular half */
    int sizeL() const;

    /** \brief Get the number of non-zeros in the upper triangular half */
    int sizeU() const;

    /** \brief Get get the number of non-zeros on the diagonal */
    int sizeD() const;

    /** \brief Get the number of elements */
    int numel() const;

    /** \brief Get the first dimension (i.e. number of rows) */
    int size1() const;

    /** \brief Get the second dimension (i.e. number of columns) */
    int size2() const;

    /** \brief Get the number if non-zeros for a given sparsity pattern */
    int size(SparsityType sp) const;

    /** \brief Get string representation of dimensions.
        The representation is (nrow x ncol = numel | size)
    */
    std::string dimString() const;

#ifndef SWIG
    /** \brief  Get the shape */
    std::pair<int, int> shape() const;
#endif

    /** \brief Check if the sparsity is empty, i.e. if one of the dimensions is zero
     * (or optionally both dimensions) */
    bool isEmpty(bool both=false) const { return sparsity().isEmpty(both);}

    /** \brief  Check if the matrix expression is dense */
    bool isDense() const { return sparsity().isDense();}

    /** \brief  Check if the matrix expression is scalar */
    bool isScalar(bool scalar_and_dense=false) const;

    /** \brief  Check if the matrix expression is square */
    bool isSquare() const { return sparsity().isSquare();}

    /** \brief  Check if the matrix is a vector (i.e. size2()==1) */
    bool isVector() const { return sparsity().isVector();}

    /** \brief Check if the matrix is upper triangular */
    bool isTriu() const { return sparsity().isTriu();}

    /** \brief Check if the matrix is lower triangular */
    bool isTril() const { return sparsity().isTril();}

    /** \brief Get the sparsity pattern */
    const Sparsity& sparsity() const;

    /** \brief Access the sparsity, make a copy if there are multiple references to it */
    Sparsity& sparsityRef();

#ifndef SWIG
    /** \brief  Get vector nonzero or slice of nonzeros */
    template<typename K>
    const MatType operator[](const K& k) const
    { return static_cast<const MatType*>(this)->getNZ(k); }

    /** \brief  Access vector nonzero or slice of nonzeros */
    template<typename K>
    NonZeros<MatType, K> operator[](const K& k)
    { return NonZeros<MatType, K>(static_cast<MatType&>(*this), k); }

    /** \brief  Get vector element or slice */
    template<typename RR>
    const MatType operator()(const RR& rr) const
    { return static_cast<const MatType*>(this)->sub(rr, 0);}

    /** \brief  Get Sparsity slice */
    const MatType operator()(const Sparsity& sp) const
    { return static_cast<const MatType*>(this)->sub(sp); }

    /** \brief  Get Matrix element or slice */
    template<typename RR, typename CC>
    const MatType operator()(const RR& rr, const CC& cc) const
    { return static_cast<const MatType*>(this)->sub(rr, cc); }

    /** \brief  Access vector element or slice */
    template<typename RR>
    SubMatrix<MatType, RR, int> operator()(const RR& rr)
    { return SubMatrix<MatType, RR, int>(static_cast<MatType&>(*this), rr, 0); }

    /** \brief  Access Sparsity slice */
    SubMatrix<MatType, Sparsity, int> operator()(const Sparsity& sp)
    { return SubMatrix<MatType, Sparsity, int>(static_cast<MatType&>(*this), sp, 0); }

    /** \brief  Access Matrix element or slice */
    template<typename RR, typename CC>
    SubMatrix<MatType, RR, CC> operator()(const RR& rr, const CC& cc)
    { return SubMatrix<MatType, RR, CC>(static_cast<MatType&>(*this), rr, cc); }
#endif // SWIG

    /** @name Construct symbolic primitives
        The "sym" function is intended to work in a similar way as "sym" used
        in the Symbolic Toolbox for Matlab but instead creating a
        CasADi symbolic primitive.
    */
    ///@{

    /** \brief Create an nrow-by-ncol symbolic primitive */
    static MatType sym(const std::string& name, int nrow=1, int ncol=1)
    { return sym(name, Sparsity::dense(nrow, ncol));}

    /** \brief  Construct a symbolic primitive with given dimensions */
    static MatType sym(const std::string& name, const std::pair<int, int> &rc)
    { return sym(name, rc.first, rc.second);}

    /** \brief Create symbolic primitive with a given sparsity pattern */
    static MatType sym(const std::string& name, const Sparsity& sp);

    /** \brief Create a vector of length p with with matrices
     * with symbolic primitives of given sparsity */
    static std::vector<MatType > sym(const std::string& name, const Sparsity& sp, int p);

    /** \brief Create a vector of length p with nrow-by-ncol symbolic primitives */
    static std::vector<MatType > sym(const std::string& name, int nrow, int ncol, int p)
    { return sym(name, Sparsity::dense(nrow, ncol), p);}

    /** \brief Create a vector of length r of vectors of length p with
     * symbolic primitives with given sparsity*/
    static std::vector<std::vector<MatType> >
      sym(const std::string& name, const Sparsity& sp, int p, int r);

    /** \brief Create a vector of length r of vectors of length p
     * with nrow-by-ncol symbolic primitives */
    static std::vector<std::vector<MatType> >
      sym(const std::string& name, int nrow, int ncol, int p, int r)
    { return sym(name, Sparsity::dense(nrow, ncol), p, r);}
    ///@}

    ///@{
    /** \brief  create a sparse matrix with all zeros */
    static MatType sparse(int nrow=1, int ncol=1) { return MatType(Sparsity::sparse(nrow, ncol));}
    static MatType sparse(const std::pair<int, int>& rc) { return sparse(rc.first, rc.second);}
    ///@}

    ///@{
    /** \brief Create a dense matrix or a matrix with specified sparsity with all entries zero */
    static MatType zeros(int nrow=1, int ncol=1) { return zeros(Sparsity::dense(nrow, ncol)); }
    static MatType zeros(const Sparsity& sp) { return MatType(sp, 0);}
    static MatType zeros(const std::pair<int, int>& rc) { return zeros(rc.first, rc.second);}
    ///@}

    ///@{
    /** \brief Create a dense matrix or a matrix with specified sparsity with all entries one */
    static MatType ones(int nrow=1, int ncol=1) { return ones(Sparsity::dense(nrow, ncol)); }
    static MatType ones(const Sparsity& sp) { return MatType(sp, 1);}
    static MatType ones(const std::pair<int, int>& rc) { return ones(rc.first, rc.second);}
    ///@}

    /** \brief Matrix-matrix multiplication.
     * Attempts to identify quick returns on matrix-level and
     * delegates to MatType::mul_full if no such quick returns are found.
     */
    MatType mul_smart(const MatType& y, const Sparsity& sp_z) const;

  };

#ifndef SWIG
#ifdef casadi_core_implementation
  // Implementations

  template<typename MatType>
  const Sparsity& GenericMatrix<MatType>::sparsity() const {
    return static_cast<const MatType*>(this)->sparsity();
  }

  template<typename MatType>
  Sparsity& GenericMatrix<MatType>::sparsityRef() {
    return static_cast<MatType*>(this)->sparsityRef();
  }

  template<typename MatType>
  int GenericMatrix<MatType>::size() const {
    return sparsity().size();
  }

  template<typename MatType>
  int GenericMatrix<MatType>::sizeL() const {
    return sparsity().sizeL();
  }

  template<typename MatType>
  int GenericMatrix<MatType>::sizeU() const {
    return sparsity().sizeU();
  }

  template<typename MatType>
  int GenericMatrix<MatType>::sizeD() const {
    return sparsity().sizeD();
  }

  template<typename MatType>
  int GenericMatrix<MatType>::numel() const {
    return sparsity().numel();
  }

  template<typename MatType>
  int GenericMatrix<MatType>::size1() const {
    return sparsity().size1();
  }

  template<typename MatType>
  int GenericMatrix<MatType>::size2() const {
    return sparsity().size2();
  }

  template<typename MatType>
  std::pair<int, int> GenericMatrix<MatType>::shape() const {
    return sparsity().shape();
  }

  template<typename MatType>
  std::string GenericMatrix<MatType>::dimString() const {
    return sparsity().dimString();
  }

  template<typename MatType>
  bool GenericMatrix<MatType>::isScalar(bool scalar_and_dense) const {
    return sparsity().isScalar(scalar_and_dense);
  }

  template<typename MatType>
  MatType GenericMatrix<MatType>::mul_smart(const MatType& y, const Sparsity &sp_z) const {
    const MatType& x = *static_cast<const MatType*>(this);

    if (!(x.isScalar() || y.isScalar())) {
      casadi_assert_message(size2()==y.size1(),
                            "Matrix product with incompatible dimensions. Lhs is "
                            << dimString() << " and rhs is " << y.dimString() << ".");
    }

    // Check if we can simplify the product
    if (x.isIdentity()) {
      return y;
    } else if (y.isIdentity()) {
      return x;
    } else if (x.isZero() || y.isZero()) {
      // See if one of the arguments can be used as result
      if (y.size()==0 && x.size2()==x.size1()) {
        return y;
      } else if (x.size()==0 && y.size2()==y.size1()) {
        return x;
      } else {
        if (y.size()==0 || x.size()==0 || x.isEmpty() || y.isEmpty()) {
          return MatType::sparse(x.size1(), y.size2());
        } else {
          return MatType::zeros(x.size1(), y.size2());
        }
      }
    } else if (x.isScalar() || y.isScalar()) {
      return x*y;
    } else {
      return x.mul_full(y, sp_z);
    }
  }

  template<typename MatType>
  int GenericMatrix<MatType>::size(SparsityType sp) const {
    if (sp==SPARSE) {
      return size();
    } else if (sp==SPARSESYM) {
      return sizeU();
    } else if (sp==DENSE) {
      return numel();
    } else if (sp==DENSESYM) {
      return (numel()+size2())/2;
    } else {
      throw CasadiException("Matrix<T>::size(Sparsity): unknown sparsity");
    }
  }

#endif
#endif // SWIG

#ifdef casadi_core_implementation
  template<typename MatType>
  std::vector<MatType> GenericMatrix<MatType>::sym(const std::string& name,
                                                   const Sparsity& sp, int p) {
    std::vector<MatType> ret(p);
    std::stringstream ss;
    for (int k=0; k<p; ++k) {
      ss.str("");
      ss << name << k;
      ret[k] = sym(ss.str(), sp);
    }
    return ret;
  }

  template<typename MatType>
  std::vector<std::vector<MatType> > GenericMatrix<MatType>::sym(const std::string& name,
                                                                 const Sparsity& sp, int p, int r) {
    std::vector<std::vector<MatType> > ret(r);
    for (int k=0; k<r; ++k) {
      std::stringstream ss;
      ss << name << "_" << k;
      ret[k] = sym(ss.str(), sp, p);
    }
    return ret;
  }

  template<typename MatType>
  MatType GenericMatrix<MatType>::sym(const std::string& name, const Sparsity& sp) {
    throw CasadiException("\"sym\" not defined for instantiation");
  }
#endif

} // namespace casadi

#endif // CASADI_GENERIC_MATRIX_HPP
