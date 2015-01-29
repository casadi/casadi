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
#include "../runtime/runtime.hpp"
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
  class CASADI_EXPORT NonZeroIterator :
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
  class CASADI_EXPORT Matrix :
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

    /** \brief Create a sparse matrix with all structural zeros */
    Matrix(int nrow, int ncol);

#ifndef SWIG
    /** \brief Create a sparse matrix with all structural zeros */
    explicit Matrix(const std::pair<int, int>& rc);
#endif // SWIG

    /** \brief Sparse matrix with a given sparsity and zero entries
        Alias for Matrix::zeros(sparsity)
     */
    explicit Matrix(const Sparsity& sp);

    /** \brief Construct matrix with a given sparsity and nonzeros */
    Matrix(const Sparsity& sp, const Matrix<DataType>& d);

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
     * Note: above remark applies only to C++, not Python or MATLAB interfaces
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
    using B::nnz;
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
    using B::colind;
    using B::row;
    using B::dimString;
    using B::sym;
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
    Matrix(const Matrix<A>& x) : sparsity_(x.sparsity()), data_(std::vector<DataType>(x.nnz())) {
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
        if (k<0) k+=nnz();
        return data().at(k);
      } catch(std::out_of_range& /* unnamed */) {
        std::stringstream ss;
        ss << "Out of range error in Matrix<>::at: " << k << " not in range [0, " << nnz() << ")";
        throw CasadiException(ss.str());
      }
    }
#else // SWIG
    /// Access a non-zero element
    DataType at(int k) {
      try {
        if (k<0) k+=nnz();
        return data().at(k);
      } catch(std::out_of_range& /* unnamed */) {
        std::stringstream ss;
        ss << "Out of range error in Matrix<>::at: " << k << " not in range [0, " << nnz() << ")";
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

    /// Is the Matrix a Slice (only for IMatrix)
    bool isSlice(bool ind1=false) const;

    ///  Convert to Slice (only for IMatrix)
    Slice toSlice(bool ind1=false) const;

    ///@{
    /// Get a submatrix, single argument
    const Matrix<DataType> getSub(bool ind1, const Slice& rr) const;
    const Matrix<DataType> getSub(bool ind1, const Matrix<int>& rr) const;
    const Matrix<DataType> getSub(bool ind1, const Sparsity& sp) const;
    ///@}

    /// Get a submatrix, two arguments
    ///@{
    const Matrix<DataType> getSub(bool ind1, const Slice& rr, const Slice& cc) const;
    const Matrix<DataType> getSub(bool ind1, const Slice& rr, const Matrix<int>& cc) const;
    const Matrix<DataType> getSub(bool ind1, const Matrix<int>& rr, const Slice& cc) const;
    const Matrix<DataType> getSub(bool ind1, const Matrix<int>& rr, const Matrix<int>& cc) const;
    ///@}

    ///@{
    /// Set a submatrix, single argument
    void setSub(const Matrix<DataType>& m, bool ind1, const Slice& rr);
    void setSub(const Matrix<DataType>& m, bool ind1, const Matrix<int>& rr);
    void setSub(const Matrix<DataType>& m, bool ind1, const Sparsity& sp);
    ///@}

    ///@{
    /// Set a submatrix, two arguments
    void setSub(const Matrix<DataType>& m, bool ind1, const Slice& rr, const Slice& cc);
    void setSub(const Matrix<DataType>& m, bool ind1, const Slice& rr, const Matrix<int>& cc);
    void setSub(const Matrix<DataType>& m, bool ind1, const Matrix<int>& rr, const Slice& cc);
    void setSub(const Matrix<DataType>& m, bool ind1, const Matrix<int>& rr, const Matrix<int>& cc);
    ///@}

    ///@{
    /// Add a submatrix to an existing matrix (TODO: remove memory allocation)
    template<typename RR, typename CC>
    void addSub(const Matrix<DataType>& m, RR rr, CC cc, bool ind1) {
      setSub(m+sub(rr, cc, ind1), rr, cc, ind1);
    }
    ///@}

    ///@{
    /// Get a set of nonzeros
    const Matrix<DataType> getNZ(bool ind1, const Slice& k) const;
    const Matrix<DataType> getNZ(bool ind1, const Matrix<int>& k) const;
    ///@}

    ///@{
    /// Set a set of nonzeros
    void setNZ(const Matrix<DataType>& m, bool ind1, const Slice& k);
    void setNZ(const Matrix<DataType>& m, bool ind1, const Matrix<int>& k);
    ///@}

    /// Append a matrix vertically (NOTE: only efficient if vector)
    void append(const Matrix<DataType>& y);

    /// Append a matrix horizontally
    void appendColumns(const Matrix<DataType>& y);

    /// Set all elements to zero
    void setZero();

    /// Set all elements to a value
    void setAll(const DataType& val);

    /** \brief Set sparse */
    Matrix<DataType> setSparse(const Sparsity& sp, bool intersect=false) const;

    /// Make the matrix dense
    void makeDense(const DataType& val = 0);

    /** \brief  Make a matrix sparse by removing numerical zeros smaller
     * in absolute value than a specified tolerance */
    void makeSparse(double tol=0);

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

    /// \cond CLUTTER
    ///@{
    /// Functions called by the corresponding friend functions -- MATLAB naming
    Matrix<DataType> zz_plus(const Matrix<DataType> &y) const;
    Matrix<DataType> zz_minus(const Matrix<DataType> &y) const;
    Matrix<DataType> zz_times(const Matrix<DataType> &y) const;
    Matrix<DataType> zz_rdivide(const Matrix<DataType> &y) const;
    Matrix<DataType> zz_lt(const Matrix<DataType> &y) const;
    Matrix<DataType> zz_le(const Matrix<DataType> &y) const;
    Matrix<DataType> zz_eq(const Matrix<DataType> &y) const;
    Matrix<DataType> zz_ne(const Matrix<DataType> &y) const;
    Matrix<DataType> __truediv__(const Matrix<DataType> &y) const {return zz_rdivide(y);}
    Matrix<DataType> zz_power(const Matrix<DataType> &y) const;
    Matrix<DataType> __constpow__(const Matrix<DataType> &y) const;
    Matrix<DataType> zz_mpower(const Matrix<DataType> &y) const;
    Matrix<DataType> __mrdivide__(const Matrix<DataType> &y) const;
    bool zz_isEqual(const Matrix<DataType> &ex2, int depth=0) const;
    void zz_expand(Matrix<DataType> &weights, Matrix<DataType>& terms) const;
    Matrix<DataType> zz_pw_const(const Matrix<DataType> &tval, const Matrix<DataType> &val) const;
    Matrix<DataType> zz_pw_lin(const Matrix<DataType> &tval, const Matrix<DataType> &val) const;
    Matrix<DataType> zz_if_else(const Matrix<DataType> &if_true,
                                const Matrix<DataType> &if_false) const;
    Matrix<DataType> zz_heaviside() const;
    Matrix<DataType> zz_rectangle() const;
    Matrix<DataType> zz_triangle() const;
    Matrix<DataType> zz_ramp() const;
    Matrix<DataType> zz_gauss_quadrature(const Matrix<DataType> &x, const Matrix<DataType> &a,
                                         const Matrix<DataType> &b, int order=5) const {
      return zz_gauss_quadrature(x, a, b, order, Matrix<DataType>());
    }
    Matrix<DataType> zz_gauss_quadrature(const Matrix<DataType> &x, const Matrix<DataType> &a,
                                         const Matrix<DataType> &b, int order,
                                         const Matrix<DataType>& w) const;
    Matrix<DataType> zz_simplify() const;
    Matrix<DataType> zz_substitute(const Matrix<DataType>& v, const Matrix<DataType>& vdef) const;
    static std::vector<Matrix<DataType> > zz_substitute(const std::vector<Matrix<DataType> >& ex,
                                                        const std::vector<Matrix<DataType> >& v,
                                                        const std::vector<Matrix<DataType> >& vdef);
    static void zz_substituteInPlace(const std::vector<Matrix<DataType> >& v,
                                     std::vector<Matrix<DataType> >& vdef,
                                     std::vector<Matrix<DataType> >& ex,
                                     bool reverse=false);
    Matrix<DataType> zz_spy() const;
    bool zz_dependsOn(const Matrix<DataType> &arg) const;
    std::vector<Matrix<DataType> > zz_getSymbols() const;
    static std::vector<Matrix<DataType> > zz_getSymbols(const std::vector<Matrix<DataType> >& e);
    Matrix<DataType> zz_jacobian(const Matrix<DataType> &arg) const;
    Matrix<DataType> zz_gradient(const Matrix<DataType> &arg) const;
    Matrix<DataType> zz_tangent(const Matrix<DataType> &arg) const;
    Matrix<DataType> zz_hessian(const Matrix<DataType> &arg) const;
    void zz_hessian(const Matrix<DataType> &arg, Matrix<DataType> &H, Matrix<DataType> &g) const;
    Matrix<DataType> zz_jacobianTimesVector(const Matrix<DataType> &arg, const Matrix<DataType> &v,
                                            bool transpose_jacobian=false) const;
    Matrix<DataType> zz_taylor(const Matrix<DataType>& x,
                               const Matrix<DataType>& a=0, int order=1) const;
    Matrix<DataType> zz_mtaylor(const Matrix<DataType>& x,
                                const Matrix<DataType>& a, int order=1) const;
    Matrix<DataType> zz_mtaylor(const Matrix<DataType>& x, const Matrix<DataType>& a, int order,
                                const std::vector<int>& order_contributions) const;
    int zz_countNodes() const;
    std::string zz_getOperatorRepresentation(const std::vector<std::string>& args) const;
    static void zz_extractShared(std::vector<Matrix<DataType> >& ex,
                                 std::vector<Matrix<DataType> >& v,
                                 std::vector<Matrix<DataType> >& vdef,
                                 const std::string& v_prefix="v_",
                                 const std::string& v_suffix="");
    void zz_printCompact(std::ostream &stream=CASADI_COUT) const;
    Matrix<DataType> zz_poly_coeff(const Matrix<DataType>&x) const;
    Matrix<DataType> zz_poly_roots() const;
    Matrix<DataType> zz_eig_symbolic() const;
    Matrix<DataType> zz_sparsify(double tol=0) const;
    ///@}

    /** \brief Propagate sparsity using 0-1 logic through a matrix product,
     * no memory allocation: <tt>z = mul(x, y)</tt> with work vector
     */
    template<bool Fwd>
    static void mul_sparsity(Matrix<DataType> &x, Matrix<DataType> &y, Matrix<DataType>& z,
                             std::vector<DataType>& work);

    /// Calculates inner_prod(x, mul(A, x)) without memory allocation
    static DataType quad_form(const std::vector<DataType>& x, const Matrix<DataType>& A);
    /// \endcond

    /// Transpose the matrix
    Matrix<DataType> T() const;
    ///@{

    /// \cond CLUTTER

    ///@{
    /// Operations called by the corresponding friend functions, MATLAB naming convention
    Matrix<DataType> zz_sin() const;
    Matrix<DataType> zz_cos() const;
    Matrix<DataType> zz_tan() const;
    Matrix<DataType> zz_asin() const;
    Matrix<DataType> zz_acos() const;
    Matrix<DataType> zz_atan() const;
    Matrix<DataType> zz_exp() const;
    Matrix<DataType> zz_log() const;
    Matrix<DataType> zz_sqrt() const;
    Matrix<DataType> zz_floor() const;
    Matrix<DataType> zz_ceil() const;
    Matrix<DataType> zz_mod(const Matrix<DataType>& y) const;
    Matrix<DataType> zz_abs() const;
    Matrix<DataType> zz_sign() const;
    Matrix<DataType> __copysign__(const Matrix<DataType>& y) const;
    Matrix<DataType> zz_erfinv() const;
    Matrix<DataType> zz_min(const Matrix<DataType>& y) const;
    Matrix<DataType> zz_max(const Matrix<DataType>& y) const;
    Matrix<DataType> zz_erf() const;
    Matrix<DataType> zz_sinh() const;
    Matrix<DataType> zz_cosh() const;
    Matrix<DataType> zz_tanh() const;
    Matrix<DataType> zz_asinh() const;
    Matrix<DataType> zz_acosh() const;
    Matrix<DataType> zz_atanh() const;
    Matrix<DataType> zz_atan2(const Matrix<DataType>& y) const;
    Matrix<DataType> zz_log10() const;
    Matrix<DataType> printme(const Matrix<DataType>& y) const;
    Matrix<DataType> zz_not() const;
    Matrix<DataType> zz_and(const Matrix<DataType>& y) const;
    Matrix<DataType> zz_or(const Matrix<DataType>& y) const;
    Matrix<DataType> zz_if_else_zero(const Matrix<DataType>& y) const;
    Matrix<DataType> zz_mtimes(const Matrix<DataType> &y) const;
    Matrix<DataType> zz_mtimes(const Matrix<DataType> &y, const Matrix<DataType> &z) const;
    Matrix<DataType> zz_det() const;
    Matrix<DataType> zz_sumAll() const;
    Matrix<DataType> zz_sumCols() const;
    Matrix<DataType> zz_sumRows() const;
    Matrix<DataType> zz_adj() const;
    Matrix<DataType> zz_inv() const;
    Matrix<DataType> zz_cofactor(int i, int j) const;
    Matrix<DataType> zz_getMinor(int i, int j) const;
    Matrix<DataType> zz_reshape(int nrow, int ncol) const;
    Matrix<DataType> zz_reshape(const Sparsity& sp) const;
    Matrix<DataType> zz_trace() const;
    Matrix<DataType> zz_vecNZ() const;
    static Matrix<DataType> zz_blockcat(const std::vector< std::vector<Matrix<DataType> > > &v);
    static Matrix<DataType> zz_horzcat(const std::vector<Matrix<DataType> > &v);
    std::vector<Matrix<DataType> > zz_horzsplit(const std::vector<int>& offset) const;
    static Matrix<DataType> zz_vertcat(const std::vector<Matrix<DataType> > &v);
    std::vector< Matrix<DataType> > zz_vertsplit(const std::vector<int>& offset) const;
    std::vector< Matrix<DataType> > zz_diagsplit(const std::vector<int>& offset1,
                                                 const std::vector<int>& offset2) const;
    Matrix<DataType> zz_inner_prod(const Matrix<DataType> &y) const;
    Matrix<DataType> zz_outer_prod(const Matrix<DataType> &y) const;
    Matrix<DataType> zz_all() const;
    Matrix<DataType> zz_any() const;
    Matrix<DataType> zz_norm_1() const;
    Matrix<DataType> zz_norm_2() const;
    Matrix<DataType> zz_norm_F() const;
    Matrix<DataType> zz_norm_inf() const;
    void zz_qr(Matrix<DataType>& Q, Matrix<DataType> &R) const;
    Matrix<DataType> zz_nullspace() const;
    Matrix<DataType> zz_solve(const Matrix<DataType>& b) const;
    Matrix<DataType> zz_pinv() const;
    Matrix<DataType> zz_kron(const Matrix<DataType>& b) const;
    Matrix<DataType> zz_repmat(int n, int m=1) const;
    Matrix<DataType> zz_repmat(const Sparsity& sp) const;
    Matrix<DataType> zz_diag() const;
    static Matrix<DataType> zz_diagcat(const std::vector< Matrix<DataType> > &A);
    Matrix<DataType> zz_unite(const Matrix<DataType>& B) const;
    Matrix<DataType> zz_polyval(const Matrix<DataType>& x) const;
    void zz_addMultiple(const std::vector<DataType>& v,
                        std::vector<DataType>& res, bool trans_A=false) const;
    Matrix<DataType> zz_project(const Sparsity& sparsity) const;
    Matrix<DataType> zz_norm_inf_mul(const Matrix<DataType> &y) const;
    ///@}

    /// \endcond

/**
\ingroup expression_tools
@{
*/
#if !defined(SWIG) || defined(DOXYGEN)
    /** \brief Matrix adjoint */
    inline friend Matrix<DataType> adj(const Matrix<DataType>& A) { return A.zz_adj();}

    /** \brief Get the (i,j) minor matrix */
    inline friend Matrix<DataType> getMinor(const Matrix<DataType> &x, int i, int j) {
      return x.zz_getMinor(i, j);
    }

    /** \brief Get the (i,j) cofactor matrix */
    inline friend Matrix<DataType> cofactor(const Matrix<DataType> &x, int i, int j) {
      return x.zz_cofactor(i, j);
    }

    /** \brief  QR factorization using the modified Gram-Schmidt algorithm
     * More stable than the classical Gram-Schmidt, but may break down if the rows of A
     * are nearly linearly dependent
     * See J. Demmel: Applied Numerical Linear Algebra (algorithm 3.1.).
     * Note that in SWIG, Q and R are returned by value. */
    inline friend void qr(const Matrix<DataType>& A, Matrix<DataType>& Q, Matrix<DataType>& R) {
      return A.zz_qr(Q, R);
    }

    /** \brief Create a new matrix with a given sparsity pattern but with the
     * nonzeros taken from an existing matrix */
    inline friend Matrix<DataType> project(const Matrix<DataType>& A, const Sparsity& sp) {
      return A.zz_project(sp);
    }

    /// Returns true only if every element in the matrix is true
    inline friend Matrix<DataType> all(const Matrix<DataType> &x) { return x.zz_all();}

    /// Returns true if any element in the matrix is true
    inline friend Matrix<DataType> any(const Matrix<DataType> &x) { return x.zz_any();}

    /** Inf-norm of a Matrix-Matrix product */
    inline friend Matrix<DataType>
      norm_inf_mul(const Matrix<DataType> &x, const Matrix<DataType> &y) {
      return x.zz_norm_inf_mul(y);
    }

    /// same as: res += mul(A, v)
    inline friend void addMultiple(const Matrix<DataType>& A,
                                   const std::vector<DataType>& v,
                                   std::vector<DataType>& res,
                                   bool trans_A=false) {
      return A.zz_addMultiple(v, res, trans_A);
    }

    /** \brief  Make a matrix sparse by removing numerical zeros*/
    inline friend Matrix<DataType> sparsify(const Matrix<DataType>& A, double tol=0) {
      return A.zz_sparsify(tol);
    }
#endif // !SWIG || DOXYGEN
/** @} */

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
    void print(std::ostream &stream=CASADI_COUT, bool trailing_newline=true) const;

    /// Print a representation of the object
    void repr(std::ostream &stream=CASADI_COUT, bool trailing_newline=true) const;

    /// Print scalar
    void printScalar(std::ostream &stream=CASADI_COUT, bool trailing_newline=true) const;

    /// Print vector-style
    void printVector(std::ostream &stream=CASADI_COUT, bool trailing_newline=true) const;

    /// Print dense matrix-stype
    void printDense(std::ostream &stream=CASADI_COUT, bool trailing_newline=true) const;

    /// Print sparse matrix style
    void printSparse(std::ostream &stream=CASADI_COUT, bool trailing_newline=true) const;

    void clear();
    void resize(int nrow, int ncol);
    void reserve(int nnz);
    void reserve(int nnz, int ncol);

    /** \brief Erase a submatrix (leaving structural zeros in its place)
        Erase rows and/or columns of a matrix */
    void erase(const std::vector<int>& rr, const std::vector<int>& cc, bool ind1=false);

    /** \brief Erase a submatrix (leaving structural zeros in its place)
        Erase elements of a matrix */
    void erase(const std::vector<int>& rr, bool ind1=false);

    /** \brief Remove columns and rows
        Remove/delete rows and/or columns of a matrix */
    void remove(const std::vector<int>& rr, const std::vector<int>& cc);

    /** \brief Enlarge matrix
        Make the matrix larger by inserting empty rows and columns,
        keeping the existing non-zeros */
    void enlarge(int nrow, int ncol,
                 const std::vector<int>& rr, const std::vector<int>& cc, bool ind1=false);

    /// Access the non-zero elements
    std::vector<DataType>& data();

    /// Const access the non-zero elements
    const std::vector<DataType>& data() const;

#ifndef SWIG
    /// \cond INTERNAL
    /// Get a pointer to the data
    DataType* ptr() { return isEmpty() ? static_cast<DataType*>(0) : &front();}
    friend inline DataType* getPtr(Matrix<DataType>& v) { return v.ptr();}

    /// Get a const pointer to the data
    const DataType* ptr() const { return isEmpty() ? static_cast<const DataType*>(0) : &front();}
    friend inline const DataType* getPtr(const Matrix<DataType>& v) { return v.ptr();}
    /// \endcond
#endif // SWIG

    /// Get the non-zero elements
    Matrix<DataType> nonzeros() const { return data();}

    /// Const access the sparsity - reference to data member
    const Sparsity& sparsity() const { return sparsity_; }

    /// Access the sparsity, make a copy if there are multiple references to it
    Sparsity& sparsityRef();

    /// \cond INTERNAL
    /** \brief  Set the non-zero elements, scalar */
    void set(DataType val, SparsityType sp=SP_SPARSE);

    /** \brief  Get the non-zero elements, scalar */
    void get(DataType& val, SparsityType sp=SP_SPARSE) const;

    /** \brief  Set the non-zero elements, vector */
    void set(const std::vector<DataType>& val, SparsityType sp=SP_SPARSE);

    /** \brief  Get the non-zero elements, vector */
    void get(std::vector<DataType>& val, SparsityType sp=SP_SPARSE) const;

    /** \brief  Set the non-zero elements, Matrix */
    void set(const Matrix<DataType>& val, SparsityType sp=SP_SPARSE);

    /** \brief  Get the non-zero elements, Matrix */
    void get(Matrix<DataType>& val, SparsityType sp=SP_SPARSE) const;

#ifdef SWIG
    %rename(get) getStridedArray;
    %rename(set) setArray;
#endif

    /** \brief  Get the non-zero elements, array */
    void getArray(DataType* val, int len, SparsityType sp=SP_SPARSE) const;

    /** \brief  Set the non-zero elements, array */
    void setArray(const DataType* val, int len, SparsityType sp=SP_SPARSE);

    /** \brief  Get the non-zero elements, array, sparse and correct length */
    void getArray(DataType* val) const;

    /** \brief  Set the non-zero elements, array, sparse and correct length */
    void setArray(const DataType* val);

    /** \brief  Get the non-zero elements, strided array */
    void getStridedArray(DataType* val, int len, int stride1, int stride2,
                         SparsityType sp=SP_SPARSE) const;

#ifndef SWIG
    /** \brief  Legacy - use getArray instead */
    void get(DataType* val, SparsityType sp=SP_SPARSE) const;

    /** \brief  Legacy - use setArray instead */
    void set(const DataType* val, SparsityType sp=SP_SPARSE);

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

    /** \brief  create an n-by-n identity matrix */
    static Matrix<DataType> eye(int ncol);

    /** \brief Returns a number that is unique for a given symbolic scalar
     * 
     * Only defined if symbolic scalar. 
     */
    long getElementHash() const;

    /// Checks if expression does not contain NaN or Inf
    bool isRegular() const;

    /** \brief Check if smooth */
    bool isSmooth() const;

    /** \brief Check if SX is a leaf of the SX graph
    
        Only defined if symbolic scalar. 
    */
    bool isLeaf() const;

    /** \brief Check whether a binary SX is commutative
    
        Only defined if symbolic scalar. 
    */
    bool isCommutative() const;

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

    /** \brief  Check if the matrix has any zero entries which are not structural zeros */
    bool hasNonStructuralZeros() const;

    /** \brief Get double value (only if constant) */
    double getValue() const;

    /** \brief Get name (only if symbolic scalar) */
    std::string getName() const;

    /** \brief Get expressions of the children of the expression 
        Only defined if symbolic scalar. 
        Wraps SXElement SXElement::getDep(int ch=0) const.
     */
    Matrix<DataType> getDep(int ch=0) const;

    /** \brief Get the number of dependencies of a binary SXElement
        Only defined if symbolic scalar. 
    */
    int getNdeps() const;

    // @{
    /// Set the 'precision, width & scientific' used in printing and serializing to streams
    static void setPrecision(int precision) { stream_precision_ = precision; }
    static void setWidth(int width) { stream_width_ = width; }
    static void setScientific(bool scientific) { stream_scientific_ = scientific; }
    // @}

#ifndef SWIG
    /// Sparse matrix with a given sparsity with all values same
    Matrix(const Sparsity& sp, const DataType& val, bool dummy);

    /// Sparse matrix with a given sparsity and non-zero elements.
    Matrix(const Sparsity& sp, const std::vector<DataType>& d, bool dummy);

  private:
    /// Sparsity of the matrix in a compressed column storage (CCS) format
    Sparsity sparsity_;

    /// Nonzero elements
    std::vector<DataType> data_;

    /// Precision used in streams
    static int stream_precision_;
    static int stream_width_;
    static bool stream_scientific_;
#endif // SWIG
  };

  // Template specialization declarations
  template<> bool Matrix<int>::isSlice(bool ind1) const;
  template<> Slice Matrix<int>::toSlice(bool ind1) const;

  // Typedefs initializations
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

#ifdef casadi_implementation
#include "matrix_impl.hpp"
#endif

#endif // CASADI_MATRIX_HPP
