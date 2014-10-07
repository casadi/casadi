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


#ifndef CASADI_MATRIX_HPP
#define CASADI_MATRIX_HPP

#include <vector>
#include <typeinfo>
#include "../casadi_exception.hpp"
#include "../printable_object.hpp"
#include "../casadi_limits.hpp"
#include "../std_vector_tools.hpp"
#include "generic_matrix.hpp"
#include "generic_expression.hpp"

namespace casadi {

/// \cond INTERNAL
  template<typename DataType>
  struct NonZero {
    int k; // Non-zero index into matrix
    int i; // Row into matrix
    int j; // Col into matrix
    DataType el;  // Element
  };


  template<typename DataType>
  class CASADI_CORE_EXPORT NonZeroIterator :
        public std::iterator< std::forward_iterator_tag, NonZero<DataType> > {
  public:
    NonZeroIterator(const Matrix<DataType> & m);

#ifndef SWIG
    NonZeroIterator<DataType>& operator++();
#endif // SWIG

    NonZero<DataType>& operator*();

    NonZeroIterator<DataType> begin();
    NonZeroIterator<DataType> end();
    bool operator==(const NonZeroIterator<DataType>& rhs);

  private:
    Matrix<DataType> m_;
    NonZero<DataType> nz;
  };
/// \endcond

/// \cond CLUTTER
  ///@{
  /** \brief Get typename */
  template <typename DataType> inline std::string matrixName()
  { return std::string("Matrix<") + typeid(DataType).name() + std::string(">");}
  template<> inline std::string matrixName<double>() { return "DMatrix"; }
  template<> inline std::string matrixName<int>() { return "IMatrix"; }
  ///@}
/// \endcond

  /** \brief Sparse matrix class. SX and DMatrix are specializations.

      General sparse matrix class that is designed with the idea that "everything is a matrix",
      that is, also scalars and vectors.\n
      This philosophy makes it easy to use and to interface in particularly
      with Python and Matlab/Octave.\n

      Index starts with 0.\n
      Index vec happens as follows: (rr, cc) -> k = rr+cc*size1()\n
      Vectors are column vectors.\n

      The storage format is Compressed Column Storage (CCS), similar to that used for
      sparse matrices in Matlab, \n
      but unlike this format, we do allow for elements to be structurally non-zero
      but numerically zero.\n

      Matrix<DataType> is polymorphic with a std::vector<DataType> that contain
      all non-identical-zero elements.\n
      The sparsity can be accessed with Sparsity& sparsity()\n

      \author Joel Andersson
      \date 2010-2014
  */
  template<typename DataType>
  class CASADI_CORE_EXPORT Matrix :
        public GenericExpression<Matrix<DataType> >,
        public GenericMatrix<Matrix<DataType> >,
        public PrintableObject<Matrix<DataType> > {
  public:

    /** \brief  constructors */
    /// empty 0-by-0 matrix constructor
    Matrix();

    /// Copy constructor
    Matrix(const Matrix<DataType>& m);

#ifndef SWIG
    /// Assignment (normal)
    Matrix<DataType>& operator=(const Matrix<DataType>& m);
#endif // SWIG

    /// Dense matrix constructor with data given as vector of vectors
    explicit Matrix(const std::vector< std::vector<DataType> >& m);

    ///@{
    /// Sparse matrix with a given sparsity
    explicit Matrix(const Sparsity& sparsity, const DataType& val=DataType(0));
    ///@}

    /// Sparse matrix with a given sparsity and non-zero elements.
    Matrix(const Sparsity& sparsity, const std::vector<DataType>& d);

    /** \brief Check if the dimensions and colind, row vectors are compatible.
     * \param complete  set to true to also check elementwise
     * throws an error as possible result
     */
    void sanityCheck(bool complete=false) const;

    /// This constructor enables implicit type conversion from a numeric type
    Matrix(double val);

    /// Construct from a vector
    /**
     * Thanks to implicit conversion, you can pretend that Matrix(const SXElement& x); exists.
     * Note: above remark applies only to C++, not python or octave interfaces
     */
    Matrix(const std::vector<DataType>& x);

    /// Construct dense matrix from a vector with the elements in column major ordering
    Matrix(const std::vector<DataType>& x, int nrow, int ncol);

    /// Convert to scalar type
    const DataType toScalar() const;

    /// Scalar type
    typedef DataType ScalarType;

#ifndef SWIG
    /// Base class
    typedef GenericMatrix<Matrix<DataType> > B;

    /// Expose base class functions
    using B::size;
    using B::sizeL;
    using B::sizeU;
    using B::numel;
    using B::size1;
    using B::size2;
    using B::shape;
    using B::isEmpty;
    using B::isScalar;
    using B::isDense;
    using B::isVector;
    using B::isTril;
    using B::isTriu;
    using B::dimString;
    using B::sym;
    using B::sparse;
    using B::zeros;
    using B::ones;
    using B::operator[];
    using B::operator();

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
#endif // SWIG

    /** \brief  Create a matrix from a matrix with a different type of matrix entries
     * (assuming that the scalar conversion is valid) */
    template<typename A>
    Matrix(const Matrix<A>& x) : sparsity_(x.sparsity()), data_(std::vector<DataType>(x.size())) {
      copy(x.begin(), x.end(), begin());
    }

    /** \brief  Create an expression from an stl vector  */
    template<typename A>
    Matrix(const std::vector<A>& x) : sparsity_(Sparsity::dense(x.size(), 1)),
        data_(std::vector<DataType>(x.size())) {
      copy(x.begin(), x.end(), begin());
    }

    /** \brief  Create a non-vector expression from an stl vector */
    template<typename A>
    Matrix(const std::vector<A>& x,  int nrow, int ncol) : sparsity_(Sparsity::dense(nrow, ncol)),
        data_(std::vector<DataType>(x.size())) {
      if (x.size() != nrow*ncol)
          throw CasadiException("Matrix::Matrix(const std::vector<DataType>& x, "
                                " int n, int m): dimension mismatch");
      copy(x.begin(), x.end(), begin());
    }


#ifndef SWIG
    /// Get a non-zero element
    inline const DataType& at(int k) const {
      return const_cast<Matrix<DataType>*>(this)->at(k);
    }

    /// Access a non-zero element
    inline DataType& at(int k) {
      try {
        if (k<0) k+=size();
        return data().at(k);
      } catch(std::out_of_range& /* unnamed */) {
        std::stringstream ss;
        ss << "Out of range error in Matrix<>::at: " << k << " not in range [0, " << size() << ")";
        throw CasadiException(ss.str());
      }
    }
#else // SWIG
    /// Access a non-zero element
    DataType at(int k) {
      try {
        if (k<0) k+=size();
        return data().at(k);
      } catch(std::out_of_range& /* unnamed */) {
        std::stringstream ss;
        ss << "Out of range error in Matrix<>::at: " << k << " not in range [0, " << size() << ")";
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

    /// Returns the truth value of a Matrix
    bool __nonzero__() const;

    /// \cond INTERNAL
    ///@{
    /// Get a submatrix
    const Matrix<DataType> sub(int rr, int cc) const;
    const Matrix<DataType> sub(const std::vector<int>& rr, int cc) const
    { return sub(rr, std::vector<int>(1, cc));}
    const Matrix<DataType> sub(int rr, const std::vector<int>& cc) const
    { return sub(std::vector<int>(1, rr), cc);}
    const Matrix<DataType> sub(const std::vector<int>& rr, const std::vector<int>& cc) const;
    const Matrix<DataType> sub(const Slice& rr, const std::vector<int>& cc) const
    { return sub(rr.getAll(size1()), cc);}
    const Matrix<DataType> sub(const std::vector<int>& rr, const Slice& cc) const
    { return sub(rr, cc.getAll(size2()));}
    const Matrix<DataType> sub(const Slice& rr, const Slice& cc) const
    { return sub(rr.getAll(size1()), cc.getAll(size2()));}
    const Matrix<DataType> sub(const Slice& rr, int cc) const
    { return sub(rr.getAll(size1()), std::vector<int>(1, cc));}
    const Matrix<DataType> sub(int rr, const Slice& cc) const
    { return sub(std::vector<int>(1, rr), cc.getAll(size2()));}
    const Matrix<DataType> sub(const Matrix<int>& rr, const std::vector<int>& cc) const;
    const Matrix<DataType> sub(const Matrix<int>& rr, int cc) const { return sub(rr, Slice(cc));}
    const Matrix<DataType> sub(const std::vector<int>& rr, const Matrix<int>& cc) const;
    const Matrix<DataType> sub(int rr, const Matrix<int>& cc) const { return sub(Slice(rr), cc);}
    const Matrix<DataType> sub(const Matrix<int>& rr, const Slice& cc) const
    {return sub(rr, cc.getAll(size2()));}
    const Matrix<DataType> sub(const Slice& rr, const Matrix<int>& cc) const
    {return sub(rr.getAll(size1()), cc);}
    const Matrix<DataType> sub(const Matrix<int>& rr, const Matrix<int>& cc) const;
    const Matrix<DataType> sub(const Sparsity& sp, int dummy = 0) const;
    ///@}

    ///@{
    /// Set a submatrix
    void setSub(const Matrix<DataType>& m, int rr, int cc);
    void setSub(const Matrix<DataType>& m, const std::vector<int>& rr, int cc)
    { setSub(m, rr, std::vector<int>(1, cc));}
    void setSub(const Matrix<DataType>& m, int rr, const std::vector<int>& cc)
    { setSub(m, std::vector<int>(1, rr), cc);}
    void setSub(const Matrix<DataType>& m, const std::vector<int>& rr, const std::vector<int>& cc);
    void setSub(const Matrix<DataType>& m, const Slice& rr, const std::vector<int>& cc)
    { setSub(m, rr.getAll(size1()), cc);}
    void setSub(const Matrix<DataType>& m, const std::vector<int>& rr, const Slice& cc)
    { setSub(m, rr, cc.getAll(size2()));}
    void setSub(const Matrix<DataType>& m, const Slice& rr, const Slice& cc)
    { setSub(m, rr.getAll(size1()), cc.getAll(size2()));}
    void setSub(const Matrix<DataType>& m, const Matrix<int>& rr, const std::vector<int>& cc);
    void setSub(const Matrix<DataType>& m, const Matrix<int>& rr, int cc)
    { setSub(m, rr, std::vector<int>(1, cc)); }
    void setSub(const Matrix<DataType>& m, const std::vector<int>& rr, const Matrix<int>& cc);
    void setSub(const Matrix<DataType>& m, int rr, const Matrix<int>& cc)
    { setSub(m, std::vector<int>(1, rr), cc); }
    void setSub(const Matrix<DataType>& m, const Matrix<int>& rr, const Slice& cc)
    {setSub(m, rr, cc.getAll(size2()));}
    void setSub(const Matrix<DataType>& m, const Slice& rr, const Matrix<int>& cc)
    {setSub(m, rr.getAll(size1()), cc);}
    void setSub(const Matrix<DataType>& m, const Matrix<int>& rr, const Matrix<int>& cc);
    void setSub(const Matrix<DataType>& m, const Slice& rr, int cc)
    {setSub(m, rr.getAll(size1()), std::vector<int>(1, cc));}
    void setSub(const Matrix<DataType>& m, const int rr, const Slice& cc)
    {setSub(m, std::vector<int>(1, rr), cc.getAll(size2()));}
    void setSub(const Matrix<DataType>& m, const Sparsity& sp, int dummy);

    ///@}



    ///@{
    /// Add a submatrix to an existing matrix (TODO: remove memory allocation)
    template<typename RR, typename CC>
    void addSub(const Matrix<DataType>& m, RR rr, CC cc) { setSub(m+sub(rr, cc), rr, cc);}
    ///@}

    ///@{
    /// Retrieve a submatrix (TODO: remove memory allocation)
    template<typename RR, typename CC>
    void getSub(Matrix<DataType>& m, RR rr, CC cc) { m = sub(rr, cc);}
    ///@}
    /// \endcond

    ///@{
    /// Get a set of nonzeros
    const Matrix<DataType> getNZ(int k) const { return at(k);}
    const Matrix<DataType> getNZ(const std::vector<int>& k) const;
    const Matrix<DataType> getNZ(const Slice& k) const { return getNZ(k.getAll(size()));}
    const Matrix<DataType> getNZ(const Matrix<int>& k) const;
    ///@}

    ///@{
    /// Set a set of nonzeros
    void setNZ(int k, const Matrix<DataType>& m);
    void setNZ(const std::vector<int>& k, const Matrix<DataType>& m);
    void setNZ(const Slice& k, const Matrix<DataType>& m) { setNZ(k.getAll(size()), m);}
    void setNZ(const Matrix<int>& k, const Matrix<DataType>& m);
    ///@}

    /// \cond INTERNAL
    /// Append a matrix vertically (NOTE: only efficient if vector)
    void append(const Matrix<DataType>& y);

    /// Append a matrix horizontally
    void appendColumns(const Matrix<DataType>& y);
    /// \endcond

    /// \cond INTERNAL
    ///@{
    /// Indexing for interfaced languages
    /// get a non-zero
    const Matrix<DataType> nz_indexed_one_based(int k) const { return this->operator[](k-1);}
    const Matrix<DataType> nz_indexed_zero_based(int k) const { return this->operator[](k);}
    const Matrix<DataType> nz_indexed_one_based(const Matrix<int>& k) const
    { return this->operator[](k-1);}
    const Matrix<DataType> nz_indexed_zero_based(const Matrix<int>& k) const
    { return this->operator[](k);}
    const Matrix<DataType> indexed_one_based(const Matrix<int>& k) const
    { return this->operator[](k-1);}
    const Matrix<DataType> indexed_zero_based(const Matrix<int>& k) const
    { return this->operator[](k);}
    const Matrix<DataType> nz_indexed(const Slice &k) const
    { return this->operator[](k);}
    const Matrix<DataType> nz_indexed(const IndexList &k) const {
      return (*this)[k.getAll(size())];
    }

    /// get a matrix element
    const Matrix<DataType> indexed_one_based(int rr, int cc) const { return (*this)(rr-1, cc-1);}
    const Matrix<DataType> indexed_zero_based(int rr, int cc) const { return (*this)(rr, cc);}
    const Matrix<DataType> indexed(const Slice &rr, const Slice &cc) const
    { return (*this)(rr, cc); }
    const Matrix<DataType> indexed(const IndexList &rr, const IndexList &cc) const {
      return (*this)(rr.getAll(size1()), cc.getAll(size2()));
    }
    const Matrix<DataType> indexed(const Slice &rr, const Matrix<int>& cc) const
    { return (*this)(rr, cc); }
    const Matrix<DataType> indexed(const Matrix<int>& rr, const IndexList &cc) const {
      return (*this)(rr, cc.getAll(size2()));
    }
    const Matrix<DataType> indexed(const Matrix<int>& rr, const Slice &cc) const
    { return (*this)(rr, cc); }
    const Matrix<DataType> indexed(const IndexList& rr, const Matrix<int> &cc) const {
      return (*this)(rr.getAll(size1()), cc);
    }
    const Matrix<DataType> indexed(const Matrix<int>& rr, const Matrix<int>& cc) const {
      return (*this)(rr, cc);
    }
    const Matrix<DataType> indexed(const Sparsity &sp) const { return (*this)(sp); }

    /// Get a vector element
    const Matrix<DataType> indexed_one_based(int rr) const {
      casadi_assert_message(isDense() && isVector(),
                            "Matrix must be a dense vector, but got "
                            << dimString() << ".");
      return (*this)(rr-1);
    }
    const Matrix<DataType> indexed_zero_based(int rr) const {
      casadi_assert_message(isDense() && isVector(),
                            "Matrix must be a dense vector, but got " << dimString() << ".");
      return (*this)(rr);
    }
    const Matrix<DataType> indexed(const Slice &rr) const {
      casadi_assert_message(isDense() && isVector(),
                            "Matrix must be a dense vector, but got " << dimString() << ".");
      return (*this)(rr);
    }
    const Matrix<DataType> indexed(const IndexList &rr) const {
      casadi_assert_message(isDense() && isVector(),
                            "Matrix must be a dense vector, but got " << dimString() << ".");
      return (*this)(rr.getAll(size1()));
    }

    /// set a non-zero
    void nz_indexed_one_based_assignment(int k, const DataType & m) { at(k-1) = m;}
    void nz_indexed_zero_based_assignment(int k, const DataType & m) { at(k) = m;}
    void nz_indexed_assignment(const Slice &k, const Matrix<DataType>& m) { (*this)[k] = m;}
    void indexed_one_based_assignment(const Matrix<int> &k, const Matrix<DataType>& m)
    { (*this)[k-1] = m;}
    void indexed_zero_based_assignment(const Matrix<int> &k, const Matrix<DataType>& m)
    { (*this)[k] = m;}
    void nz_indexed_one_based_assignment(const Matrix<int> &k, const Matrix<DataType>& m)
    { (*this)[k-1] = m;}
    void nz_indexed_zero_based_assignment(const Matrix<int> &k, const Matrix<DataType>& m)
    { (*this)[k] = m;}
    void nz_indexed_assignment(const IndexList &k, const Matrix<DataType>& m) {
      (*this)[k.getAll(size())] = m;
    }

    /// set a matrix element
    void indexed_one_based_assignment(int rr, int cc, const DataType & m) { elem(rr-1, cc-1) = m;}
    void indexed_zero_based_assignment(int rr, int cc, const DataType & m) { elem(rr, cc) = m;}
    void indexed_assignment(const Slice &rr, const Slice &cc, const Matrix<DataType>& m)
    { (*this)(rr, cc) = m; }
    void indexed_assignment(const IndexList &rr, const IndexList &cc, const Matrix<DataType>& m) {
      (*this)(rr.getAll(size1()), cc.getAll(size2())) = m;
    }
    void indexed_assignment(const Slice &rr, const Matrix<int>& cc, const Matrix<DataType>& m) {
      (*this)(rr, cc) = m;
    }
    void indexed_assignment(const Matrix<int>& rr, const Slice &cc, const Matrix<DataType>& m) {
      (*this)(rr, cc) = m;
    }
    void indexed_assignment(const Matrix<int> &rr, const IndexList& cc, const Matrix<DataType>& m) {
      (*this)(rr, cc.getAll(size2())) = m;
    }
    void indexed_assignment(const IndexList& rr, const Matrix<int> &cc, const Matrix<DataType>& m) {
      (*this)(rr.getAll(size1()), cc) = m;
    }
    void indexed_assignment(const Matrix<int>& rr, const Matrix<int>& cc,
                            const Matrix<DataType>& m) {
      (*this)(rr, cc) = m;
    }
    void indexed_assignment(const Sparsity &sp, const Matrix<DataType>& m) {
      // (*this)(sp) = m;   // VC2010 compiler errors
          SubMatrix<Matrix<DataType>, Sparsity, int> temp(*this, sp, 0);
          temp = m;
    }
    ///@}

    /// set a vector element
    void indexed_one_based_assignment(int rr, const DataType & m) {
      casadi_assert_message(isDense() && isVector(),
                            "Matrix must be a dense vector, but got " << dimString() << ".");
      elem(rr-1) = m;
    }
    void indexed_zero_based_assignment(int rr, const DataType & m) {
      casadi_assert_message(isDense() && isVector(),
                            "Matrix must be a dense vector, but got " << dimString() << ".");
      elem(rr) = m;
    }
    void indexed_assignment(const Slice &rr, const Matrix<DataType>& m) {
      casadi_assert_message(isDense() && isVector(),
                            "Matrix must be a dense vector, but got " << dimString() << ".");
      (*this)(rr, Slice(0)) = m;
    }
    void indexed_assignment(const IndexList &rr, const Matrix<DataType>& m) {
      (*this)(rr.getAll(size1()), 0) = m;
    }
    /// \endcond

    /// Set all elements to zero
    void setZero();

    /// Set all elements to a value
    void setAll(const DataType& val);

    /** \brief Set sparse */
    Matrix<DataType> setSparse(const Sparsity& sp, bool intersect=false) const;

    /// Make the matrix dense
    void densify(const DataType& val = 0);

    /** \brief  Make a matrix sparse by removing numerical zeros smaller
     * in absolute value than a specified tolerance */
    void sparsify(double tol=0);

    Matrix<DataType> operator+() const;
    Matrix<DataType> operator-() const;

    /// \cond INTERNAL
    ///@{
    /** \brief  Create nodes by their ID */
    static Matrix<DataType> binary(int op, const Matrix<DataType> &x, const Matrix<DataType> &y);
    static Matrix<DataType> unary(int op, const Matrix<DataType> &x);
    static Matrix<DataType> scalar_matrix(int op,
                                          const Matrix<DataType> &x, const Matrix<DataType> &y);
    static Matrix<DataType> matrix_scalar(int op,
                                          const Matrix<DataType> &x, const Matrix<DataType> &y);
    static Matrix<DataType> matrix_matrix(int op,
                                          const Matrix<DataType> &x, const Matrix<DataType> &y);
    ///@}
    /// \endcond

    /// \cond INTERNAL
    ///@{
    /// Elementwise operations -- Octave/Python naming
    Matrix<DataType> __add__(const Matrix<DataType> &y) const;
    Matrix<DataType> __sub__(const Matrix<DataType> &y) const;
    Matrix<DataType> __mul__(const Matrix<DataType> &y) const;
    Matrix<DataType> __div__(const Matrix<DataType> &y) const;
    Matrix<DataType> __lt__(const Matrix<DataType> &y) const;
    Matrix<DataType> __le__(const Matrix<DataType> &y) const;
    Matrix<DataType> __gt__(const Matrix<DataType> &y) const { return y.__lt__(*this);}
    Matrix<DataType> __ge__(const Matrix<DataType> &y) const { return y.__le__(*this);}
    Matrix<DataType> __eq__(const Matrix<DataType> &y) const;
    Matrix<DataType> __ne__(const Matrix<DataType> &y) const;
    Matrix<DataType> __truediv__(const Matrix<DataType> &y) const {return __div__(y);}
    Matrix<DataType> __pow__(const Matrix<DataType> &y) const;
    Matrix<DataType> __constpow__(const Matrix<DataType> &y) const;
    Matrix<DataType> __mpower__(const Matrix<DataType> &y) const;
    Matrix<DataType> __mrdivide__(const Matrix<DataType> &y) const;
    ///@}

    /// Matrix-matrix product
    Matrix<DataType> mul_full(const Matrix<DataType> &y, const Sparsity & sp_z=Sparsity()) const;

    /// Matrix-matrix product
    Matrix<DataType> mul(const Matrix<DataType> &y, const Sparsity & sp_z=Sparsity()) const;

    /// Matrix-matrix product, no memory allocation: z += mul(x, y)
    static void mul_no_alloc_nn(const Matrix<DataType> &x, const Matrix<DataType>& y,
                                Matrix<DataType>& z);

    /// Matrix-matrix product, no memory allocation: z += mul(x, y), with work vector
    static void mul_no_alloc_nn(const Matrix<DataType> &x, const Matrix<DataType>& y,
                                Matrix<DataType>& z, std::vector<DataType>& work);

    /// Matrix-matrix product, no memory allocation: z += mul(trans(x), y)
    static void mul_no_alloc_tn(const Matrix<DataType> &trans_x, const Matrix<DataType> &y,
                                Matrix<DataType>& z);

    /// Matrix-matrix product, no memory allocation: z += mul(x, trans(y))
    static void mul_no_alloc_nt(const Matrix<DataType>& x, const Matrix<DataType> &trans_y,
                                Matrix<DataType>& z);

    /// Matrix-vector product, no memory allocation: z += mul(trans(x), y)
    static void mul_no_alloc_tn(const Matrix<DataType>& trans_x, const std::vector<DataType> &y,
                                std::vector<DataType>& z);

    /// vector-matrix product, no memory allocation: z += mul(x, y)
    static void mul_no_alloc_nn(const Matrix<DataType>& x, const std::vector<DataType> &y,
                                std::vector<DataType>& z);

    /** \brief Propagate sparsity using 0-1 logic through a matrix product,
     * no memory allocation: <tt>z = mul(trans(x), y)</tt>
     */
    template<bool Fwd>
    static void mul_sparsity(Matrix<DataType> &x_trans, Matrix<DataType> &y, Matrix<DataType>& z);

    /// Calculates inner_prod(x, mul(A, x)) without memory allocation
    static DataType quad_form(const std::vector<DataType>& x, const Matrix<DataType>& A);
    /// \endcond

    /// Transpose the matrix
    Matrix<DataType> trans() const;

#ifndef SWIG
    /// Transpose the matrix (shorthand)
    Matrix<DataType> T() const { return trans();}
#endif

    ///@{

    ///@{
    /// Operations defined in the standard namespace for unambiguous access and Numpy compatibility
    Matrix<DataType> sin() const;
    Matrix<DataType> cos() const;
    Matrix<DataType> tan() const;
    Matrix<DataType> arcsin() const;
    Matrix<DataType> arccos() const;
    Matrix<DataType> arctan() const;
    Matrix<DataType> exp() const;
    Matrix<DataType> log() const;
    Matrix<DataType> sqrt() const;
    Matrix<DataType> floor() const;
    Matrix<DataType> ceil() const;
    Matrix<DataType> fmod(const Matrix<DataType>& y) const;
    Matrix<DataType> fabs() const;
    Matrix<DataType> sign() const;
    Matrix<DataType> __copysign__(const Matrix<DataType>& y) const;
    Matrix<DataType> erfinv() const;
    Matrix<DataType> fmin(const Matrix<DataType>& y) const;
    Matrix<DataType> fmax(const Matrix<DataType>& y) const;
    Matrix<DataType> erf() const;
    Matrix<DataType> sinh() const;
    Matrix<DataType> cosh() const;
    Matrix<DataType> tanh() const;
    Matrix<DataType> arcsinh() const;
    Matrix<DataType> arccosh() const;
    Matrix<DataType> arctanh() const;
    Matrix<DataType> arctan2(const Matrix<DataType>& y) const;
    Matrix<DataType> log10() const;
    Matrix<DataType> printme(const Matrix<DataType>& y) const;
    Matrix<DataType> logic_not() const;
    Matrix<DataType> logic_and(const Matrix<DataType>& y) const;
    Matrix<DataType> logic_or(const Matrix<DataType>& y) const;
    Matrix<DataType> if_else_zero(const Matrix<DataType>& y) const;
    ///@}

    /** \brief Set or reset the maximum number of calls to the
     * printing function when printing an expression */
    static void setMaxNumCallsInPrint(long num=10000);

    /** \brief Get the maximum number of calls to the printing
     * function when printing an expression */
    static long getMaxNumCallsInPrint();

    /** \brief Set or reset the depth to which equalities are being checked for simplifications */
    static void setEqualityCheckingDepth(int eq_depth=1);

    /** \brief Get the depth to which equalities are being checked for simplifications */
    static int getEqualityCheckingDepth();

    /// Get name of the class
    static std::string className();

    /// Print a description of the object
    void print(std::ostream &stream=std::cout, bool trailing_newline=true) const;

    /// Print a representation of the object
    void repr(std::ostream &stream=std::cout, bool trailing_newline=true) const;

    /// Print scalar
    void printScalar(std::ostream &stream=std::cout, bool trailing_newline=true) const;

    /// Print vector-style
    void printVector(std::ostream &stream=std::cout, bool trailing_newline=true) const;

    /// Print dense matrix-stype
    void printDense(std::ostream &stream=std::cout, bool trailing_newline=true) const;

    /// Print sparse matrix style
    void printSparse(std::ostream &stream=std::cout, bool trailing_newline=true) const;

    // Get the sparsity pattern
    const std::vector<int>& row() const;
    const std::vector<int>& colind() const;
    int row(int el) const;
    int colind(int col) const;
    void clear();
    void resize(int nrow, int ncol);
    void reserve(int nnz);
    void reserve(int nnz, int ncol);

    /** \brief Erase a submatrix
        Erase rows and/or columns of a matrix */
    void erase(const std::vector<int>& rr, const std::vector<int>& cc);

    /** \brief Remove cols or rows
        Rremove/delete rows and/or columns of a matrix */
    void remove(const std::vector<int>& rr, const std::vector<int>& cc);

    /** \brief Enlarge matrix
        Make the matrix larger by inserting empty rows and columns,
        keeping the existing non-zeros */
    void enlarge(int nrow, int ncol, const std::vector<int>& rr, const std::vector<int>& cc);

    /// Access the non-zero elements
    std::vector<DataType>& data();

    /// Const access the non-zero elements
    const std::vector<DataType>& data() const;

    /// \cond INTERNAL
    /// Get a pointer to the data
    DataType* ptr() { return isEmpty() ? static_cast<DataType*>(0) : &front();}

    /// Get a const pointer to the data
    const DataType* ptr() const { return isEmpty() ? static_cast<const DataType*>(0) : &front();}
    /// \endcond

    /// Const access the sparsity - reference to data member
    const Sparsity& sparsity() const { return sparsity_; }

    /// Access the sparsity, make a copy if there are multiple references to it
    Sparsity& sparsityRef();

    /// \cond INTERNAL
    /** \brief  Set the non-zero elements, scalar */
    void set(DataType val, SparsityType sp=SPARSE);

    /** \brief  Get the non-zero elements, scalar */
    void get(DataType& val, SparsityType sp=SPARSE) const;

    /** \brief  Set the non-zero elements, vector */
    void set(const std::vector<DataType>& val, SparsityType sp=SPARSE);

    /** \brief  Get the non-zero elements, vector */
    void get(std::vector<DataType>& val, SparsityType sp=SPARSE) const;

    /** \brief  Set the non-zero elements, Matrix */
    void set(const Matrix<DataType>& val, SparsityType sp=SPARSE);

    /** \brief  Get the non-zero elements, Matrix */
    void get(Matrix<DataType>& val, SparsityType sp=SPARSE) const;

#ifdef SWIG
    %rename(get) getStridedArray;
    %rename(set) setArray;
#endif

    /** \brief  Get the non-zero elements, array */
    void getArray(DataType* val, int len, SparsityType sp=SPARSE) const;

    /** \brief  Set the non-zero elements, array */
    void setArray(const DataType* val, int len, SparsityType sp=SPARSE);

    /** \brief  Get the non-zero elements, array, sparse and correct length */
    void getArray(DataType* val) const;

    /** \brief  Set the non-zero elements, array, sparse and correct length */
    void setArray(const DataType* val);

    /** \brief  Get the non-zero elements, strided array */
    void getStridedArray(DataType* val, int len, int stride1, int stride2,
                         SparsityType sp=SPARSE) const;

#ifndef SWIG
    /** \brief  Legacy - use getArray instead */
    void get(DataType* val, SparsityType sp=SPARSE) const;

    /** \brief  Legacy - use setArray instead */
    void set(const DataType* val, SparsityType sp=SPARSE);

    /** Bitwise set, reinterpreting the data as a bvec_t array */
    void setZeroBV();

    /** Bitwise set, reinterpreting the data as a bvec_t array */
    void setBV(const Matrix<DataType>& val);

    /** Bitwise set, reinterpreting the data as a bvec_t array */
    void getBV(Matrix<DataType>& val) const { val.setBV(*this);}

    /** Bitwise or, reinterpreting the data as a bvec_t array */
    void borBV(const Matrix<DataType>& val);

    /** \brief Bitwise get the non-zero elements, array */
    void getArrayBV(bvec_t* val, int len) const;

    /** \brief Bitwise set the non-zero elements, array */
    void setArrayBV(const bvec_t* val, int len);

    /** \brief Bitwise or the non-zero elements, array */
    void borArrayBV(const bvec_t* val, int len);

#endif
/// \endcond

    /** \brief  Save the result to the LAPACK banded format -- see LAPACK documentation
        kl:    The number of subdiagonals in res
        ku:    The number of superdiagonals in res
        ldres: The leading dimension in res
        res:   The number of superdiagonals */
    void getBand(int kl, int ku, int ldres, DataType *res) const;

    /* \brief Construct a sparse matrix from triplet form
     * Default matrix size is max(col) x max(row)
     */
    ///@{
    static Matrix<DataType> triplet(const std::vector<int>& row, const std::vector<int>& col,
                                    const std::vector<DataType>& d);
    static Matrix<DataType> triplet(const std::vector<int>& row, const std::vector<int>& col,
                                    const std::vector<DataType>& d, int nrow, int ncol);
    static Matrix<DataType> triplet(const std::vector<int>& row, const std::vector<int>& col,
                                    const std::vector<DataType>& d, const std::pair<int, int>& rc);
    ///@}

    ///@{
    /** \brief  create a matrix with all inf */
    static Matrix<DataType> inf(const Sparsity& sp);
    static Matrix<DataType> inf(int nrow=1, int ncol=1);
    static Matrix<DataType> inf(const std::pair<int, int>& rc);
    ///@}

    ///@{
    /** \brief  create a matrix with all nan */
    static Matrix<DataType> nan(const Sparsity& sp);
    static Matrix<DataType> nan(int nrow=1, int ncol=1);
    static Matrix<DataType> nan(const std::pair<int, int>& rc);
    ///@}

    ///@{
    /** \brief  create a matrix by repeating an existing matrix */
    static Matrix<DataType> repmat(const DataType& x, const Sparsity& sp);
    static Matrix<DataType> repmat(const Matrix<DataType>& x, const Sparsity& sp);
    static Matrix<DataType> repmat(const Matrix<DataType>& x, int nrow, int ncol=1);
    static Matrix<DataType> repmat(const Matrix<DataType>& x, const std::pair<int, int>& rc);
    ///@}

    /** \brief  create an n-by-n identity matrix */
    static Matrix<DataType> eye(int ncol);

    /// Checks if expression does not contain NaN or Inf
    bool isRegular() const;

    /** \brief Check if smooth */
    bool isSmooth() const;

    /** \brief Check if symbolic (Dense)
        Sparse matrices invariable return false
    */
    bool isSymbolic() const;

    /** \brief Check if symbolic
        Sparse matrices can return true if all non-zero elements are symbolic
    */
    bool isSymbolicSparse() const;

    /** \brief Check if the matrix is constant (note that false negative answers are possible)*/
    bool isConstant() const;

    /** \brief Check if the matrix is integer-valued
     * (note that false negative answers are possible)*/
    bool isInteger() const;

    /** \brief  check if the matrix is 0 (note that false negative answers are possible)*/
    bool isZero() const;

    /** \brief  check if the matrix is 1 (note that false negative answers are possible)*/
    bool isOne() const;

    /** \brief  check if the matrix is -1 (note that false negative answers are possible)*/
    bool isMinusOne() const;

    /** \brief  check if the matrix is an identity matrix (note that false negative answers
     * are possible)*/
    bool isIdentity() const;

    /** \brief Check if two expressions are equal
     *  May give false negatives
     *
     *  Note: does not work when CasadiOptions.setSimplificationOnTheFly(False) was called
     */
    bool isEqual(const Matrix<DataType> &ex2) const;

    /** \brief  Check if the matrix has any zero entries which are not structural zeros */
    bool hasNonStructuralZeros() const;

    /** \brief Get double value (only if constant) */
    double getValue() const;

    /** \brief Get name (only if symbolic scalar) */
    std::string getName() const;

    // @{
    /// Set the 'precision, width & scientific' used in printing and serializing to streams
    static void setPrecision(int precision) { stream_precision_ = precision; }
    static void setWidth(int width) { stream_width_ = width; }
    static void setScientific(bool scientific) { stream_scientific_ = scientific; }
    // @}

  private:
    /// Sparsity of the matrix in a compressed column storage (CCS) format
    Sparsity sparsity_;

    /// Nonzero elements
    std::vector<DataType> data_;

    /// Precision used in streams
    static int stream_precision_;
    static int stream_width_;
    static bool stream_scientific_;

  };

} // namespace casadi

// Typedefs/template initializations
namespace casadi {

  typedef Matrix<int> IMatrix;
  typedef Matrix<double> DMatrix;
  typedef std::vector<Matrix<double> > DMatrixVector;
  typedef std::vector< std::vector<Matrix<double> > > DMatrixVectorVector;

  /// \cond INTERNAL
  typedef DMatrix* DMatrixPtr;
  typedef std::vector<DMatrixPtr> DMatrixPtrV;
  typedef std::vector<DMatrixPtrV> DMatrixPtrVV;
  /// \endcond
} // namespace casadi

#ifdef casadi_core_implementation
#include "matrix_impl.hpp"
#endif

#endif // CASADI_MATRIX_HPP
