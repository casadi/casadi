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
#include "exception.hpp"
#include "printable_object.hpp"
#include "casadi_limits.hpp"
#include "std_vector_tools.hpp"
#include "runtime/runtime.hpp"
#include "generic_matrix.hpp"
#include "generic_expression.hpp"

namespace casadi {

  /** \brief Empty Base
      This class is extended in SWIG.
  */
  struct CASADI_EXPORT MatrixCommon {};

/// \cond CLUTTER
  ///@{
  /** \brief Get typename */
  template <typename DataType> inline std::string matrixName()
  { return std::string("Matrix<") + typeid(DataType).name() + std::string(">");}
  template<> inline std::string matrixName<double>() { return "DM"; }
  template<> inline std::string matrixName<int>() { return "IM"; }
  ///@}
/// \endcond

  /** \brief Sparse matrix class. SX and DM are specializations.

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
    public MatrixCommon,
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

    /** \brief Create a sparse matrix with all structural zeros */
    Matrix(int nrow, int ncol);

#ifndef SWIG
    /** \brief Create a sparse matrix with all structural zeros */
    explicit Matrix(const std::pair<int, int>& rc);

    /** \brief  Access functions of the node */
    std::vector<DataType>* operator->() { return &data_;}

    /** \brief  Const access functions of the node */
    const std::vector<DataType>* operator->() const { return &data_;}
#endif // SWIG

    /** \brief Create a sparse matrix from a sparsity pattern.
        Same as Matrix::ones(sparsity)
     */
    explicit Matrix(const Sparsity& sp);

    /** \brief Construct matrix with a given sparsity and nonzeros */
    Matrix(const Sparsity& sp, const Matrix<DataType>& d);

    /** \brief Check if the dimensions and colind, row vectors are compatible.
     * \param complete  set to true to also check elementwise
     * throws an error as possible result
     */
    void sanity_check(bool complete=false) const;

    /// This constructor enables implicit type conversion from a numeric type
    Matrix(double val);

    /// Dense matrix constructor with data given as vector of vectors
    explicit Matrix(const std::vector< std::vector<double> >& m);

    /** \brief Create a matrix from another matrix with a different entry type
     *  Assumes that the scalar conversion is valid.
     */
    template<typename A>
    Matrix(const Matrix<A>& x) : sparsity_(x.sparsity()), data_(std::vector<DataType>(x.nnz())) {
      copy(x->begin(), x->end(), data_.begin());
    }

    /** \brief  Create an expression from a vector  */
    template<typename A>
    Matrix(const std::vector<A>& x) : sparsity_(Sparsity::dense(x.size(), 1)),
      data_(std::vector<DataType>(x.size())) {
        copy(x.begin(), x.end(), data_.begin());
    }

#ifndef SWIG
    /// Construct from a vector
    Matrix(const std::vector<DataType>& x);

    /// Convert to scalar type
    const DataType toScalar() const;

    /// Scalar type
    typedef DataType ScalarType;

    /// Base class
    typedef GenericMatrix<Matrix<DataType> > B;

    /// Expose base class functions
    using B::nnz;
    using B::nnz_lower;
    using B::nnz_upper;
    using B::numel;
    using B::size1;
    using B::size2;
    using B::size;
    using B::is_empty;
    using B::is_scalar;
    using B::is_dense;
    using B::is_vector;
    using B::is_row;
    using B::is_column;
    using B::is_tril;
    using B::is_triu;
    using B::colind;
    using B::row;
    using B::dim;
    using B::sym;
    using B::zeros;
    using B::ones;
    using B::operator[];
    using B::operator();
    using B::horzsplit;
    using B::vertsplit;
    using B::diagsplit;
    using B::mtimes;

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

    /// get an element
    const DataType& elem(int rr, int cc=0) const;

    /// get a reference to an element
    DataType& elem(int rr, int cc=0);

    /// get an element, do not allocate
    const DataType getElement(int rr, int cc=0) const { return elem(rr, cc);}
#endif // SWIG

    /// Returns true if the matrix has a non-zero at location rr, cc
    bool hasNZ(int rr, int cc) const { return sparsity().hasNZ(rr, cc); }

    /// Returns the truth value of a Matrix
    bool __nonzero__() const;

    /// Is the Matrix a Slice (only for IM)
    bool isSlice(bool ind1=false) const;

    ///  Convert to Slice (only for IM)
    Slice toSlice(bool ind1=false) const;

    ///@{
    /// Get a submatrix, single argument
    void get(Matrix<DataType>& SWIG_OUTPUT(m), bool ind1, const Slice& rr) const;
    void get(Matrix<DataType>& SWIG_OUTPUT(m), bool ind1, const Matrix<int>& rr) const;
    void get(Matrix<DataType>& SWIG_OUTPUT(m), bool ind1, const Sparsity& sp) const;
    ///@}

    /// Get a submatrix, two arguments
    ///@{
    void get(Matrix<DataType>& SWIG_OUTPUT(m), bool ind1,
                const Slice& rr, const Slice& cc) const;
    void get(Matrix<DataType>& SWIG_OUTPUT(m), bool ind1,
                const Slice& rr, const Matrix<int>& cc) const;
    void get(Matrix<DataType>& SWIG_OUTPUT(m), bool ind1,
                const Matrix<int>& rr, const Slice& cc) const;
    void get(Matrix<DataType>& SWIG_OUTPUT(m), bool ind1,
                const Matrix<int>& rr, const Matrix<int>& cc) const;
    ///@}

    ///@{
    /// Set a submatrix, single argument
    void set(const Matrix<DataType>& m, bool ind1, const Slice& rr);
    void set(const Matrix<DataType>& m, bool ind1, const Matrix<int>& rr);
    void set(const Matrix<DataType>& m, bool ind1, const Sparsity& sp);
    ///@}

    ///@{
    /// Set a submatrix, two arguments
    void set(const Matrix<DataType>& m, bool ind1, const Slice& rr, const Slice& cc);
    void set(const Matrix<DataType>& m, bool ind1, const Slice& rr, const Matrix<int>& cc);
    void set(const Matrix<DataType>& m, bool ind1, const Matrix<int>& rr, const Slice& cc);
    void set(const Matrix<DataType>& m, bool ind1, const Matrix<int>& rr, const Matrix<int>& cc);
    ///@}

    ///@{
    /// Get a set of nonzeros
    void getNZ(Matrix<DataType>& SWIG_OUTPUT(m), bool ind1, const Slice& k) const;
    void getNZ(Matrix<DataType>& SWIG_OUTPUT(m), bool ind1, const Matrix<int>& k) const;
    ///@}

    ///@{
    /// Set a set of nonzeros
    void setNZ(const Matrix<DataType>& m, bool ind1, const Slice& k);
    void setNZ(const Matrix<DataType>& m, bool ind1, const Matrix<int>& k);
    ///@}

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

#ifndef SWIG
    /// \cond CLUTTER
    ///@{
    /// Functions called by friend functions defined for GenericExpression
    static bool is_equal(const Matrix<DataType> &x, const Matrix<DataType> &y, int depth=0);
    ///@}

    ///@{
    /// Functions called by friend functions defined for GenericMatrix
    static Matrix<DataType> simplify(const Matrix<DataType> &x);
    static Matrix<DataType> jacobian(const Matrix<DataType> &f, const Matrix<DataType> &x,
                                     bool symmetric=false);
    static Matrix<DataType> gradient(const Matrix<DataType> &f, const Matrix<DataType> &x);
    static Matrix<DataType> tangent(const Matrix<DataType> &f, const Matrix<DataType> &x);
    static Matrix<DataType> hessian(const Matrix<DataType> &f, const Matrix<DataType> &x);
    static Matrix<DataType> hessian(const Matrix<DataType> &f, const Matrix<DataType> &x,
                                    Matrix<DataType>& g);
    static Matrix<DataType>
      substitute(const Matrix<DataType>& ex,
                 const Matrix<DataType>& v,
                 const Matrix<DataType>& vdef);
    static std::vector<Matrix<DataType> >
      substitute(const std::vector<Matrix<DataType> >& ex,
                 const std::vector<Matrix<DataType> >& v,
                 const std::vector<Matrix<DataType> >& vdef);
    static void substituteInPlace(const std::vector<Matrix<DataType> >& v,
                                  std::vector<Matrix<DataType> >& vdef,
                                  std::vector<Matrix<DataType> >& ex,
                                  bool revers);
    static Matrix<DataType> pinv(const Matrix<DataType> &x);
    static Matrix<DataType> pinv(const Matrix<DataType> &A,
                                 const std::string& lsolver, const Dict& opts);
    static Matrix<DataType> solve(const Matrix<DataType> &A, const Matrix<DataType>& b);
    static Matrix<DataType> solve(const Matrix<DataType> &A, const Matrix<DataType>& b,
                                  const std::string& lsolver, const Dict& opts);
    static int countNodes(const Matrix<DataType> &x);
    static std::string print_operator(const Matrix<DataType> &x,
                                      const std::vector<std::string>& args);
    static void extractShared(std::vector<Matrix<DataType> >& ex,
                              std::vector<Matrix<DataType> >& v,
                              std::vector<Matrix<DataType> >& vdef,
                              const std::string& v_prefix,
                              const std::string& v_suffix);
    static Matrix<DataType> _bilin(const Matrix<DataType>& A,
                                   const Matrix<DataType>& x,
                                   const Matrix<DataType>& y);
    static Matrix<DataType> _rank1(const Matrix<DataType>& A,
                                   const Matrix<DataType>& alpha,
                                   const Matrix<DataType>& x,
                                   const Matrix<DataType>& y);
    static Matrix<DataType> if_else(const Matrix<DataType> &x,
                                    const Matrix<DataType> &if_true,
                                    const Matrix<DataType> &if_false,
                                    bool short_circuit);
    static Matrix<DataType> conditional(const Matrix<DataType> &ind,
                                        const std::vector<Matrix<DataType> > &x,
                                        const Matrix<DataType> &x_default,
                                        bool short_circuit);
    static bool dependsOn(const Matrix<DataType> &x, const Matrix<DataType> &arg);
    static Matrix<DataType> mpower(const Matrix<DataType> &x, const Matrix<DataType> &y);
    static Matrix<DataType> mrdivide(const Matrix<DataType> &x, const Matrix<DataType> &y);
    static Matrix<DataType> mldivide(const Matrix<DataType> &x, const Matrix<DataType> &y);
    static std::vector<Matrix<DataType> > symvar(const Matrix<DataType> &x);
    static Matrix<DataType> det(const Matrix<DataType> &x);
    static Matrix<DataType> inv(const Matrix<DataType> &x);
    static Matrix<DataType> trace(const Matrix<DataType> &x);
    static Matrix<DataType> norm_1(const Matrix<DataType> &x);
    static Matrix<DataType> norm_2(const Matrix<DataType> &x);
    static Matrix<DataType> norm_F(const Matrix<DataType> &x);
    static Matrix<DataType> norm_inf(const Matrix<DataType> &x);
    static Matrix<DataType> sumCols(const Matrix<DataType> &x);
    static Matrix<DataType> sumRows(const Matrix<DataType> &x);
    static Matrix<DataType> dot(const Matrix<DataType> &x, const Matrix<DataType> &y);
    static Matrix<DataType> nullspace(const Matrix<DataType> &x);
    static Matrix<DataType> diag(const Matrix<DataType> &x);
    static Matrix<DataType> unite(const Matrix<DataType> &A, const Matrix<DataType>& B);
    static Matrix<DataType> project(const Matrix<DataType> &x,
                                    const Sparsity& sp, bool intersect=false);
    static Matrix<DataType> polyval(const Matrix<DataType> &p, const Matrix<DataType>& x);
    static Matrix<DataType> densify(const Matrix<DataType> &x, const Matrix<DataType>& val);
    static Matrix<DataType> densify(const Matrix<DataType> &x);
    ///@}

    ///@{
    /// Functions called by friend functions defined for SparsityInterface
    static Matrix<DataType> blockcat(const std::vector< std::vector<Matrix<DataType> > > &v);
    static Matrix<DataType> horzcat(const std::vector<Matrix<DataType> > &v);
    static std::vector<Matrix<DataType> >
      horzsplit(const Matrix<DataType>& x,
                const std::vector<int>& offset);
    static Matrix<DataType> vertcat(const std::vector<Matrix<DataType> > &v);
    static std::vector< Matrix<DataType> >
      vertsplit(const Matrix<DataType>& x,
                const std::vector<int>& offset);
    static std::vector< Matrix<DataType> >
      diagsplit(const Matrix<DataType>& x,
                const std::vector<int>& offset1,
                const std::vector<int>& offset2);
    static Matrix<DataType> reshape(const Matrix<DataType> &x, int nrow, int ncol);
    static Matrix<DataType> reshape(const Matrix<DataType> &x, const Sparsity& sp);
    static Matrix<DataType> vecNZ(const Matrix<DataType> &x);
    static Matrix<DataType> kron(const Matrix<DataType> &x, const Matrix<DataType>& y);
    static Matrix<DataType> mtimes(const Matrix<DataType> &x, const Matrix<DataType> &y);
    static Matrix<DataType> mac(const Matrix<DataType> &x,
                                const Matrix<DataType> &y,
                                const Matrix<DataType> &z);
    ///@}

    ///@{
    /// Functions called by friend functions defined here
    static Matrix<DataType> sparsify(const Matrix<DataType> &x, double tol=0);
    static void expand(const Matrix<DataType>& x,
                       Matrix<DataType>& weights,
                       Matrix<DataType>& terms);
    static Matrix<DataType> pw_const(const Matrix<DataType> &t,
                                     const Matrix<DataType> &tval, const Matrix<DataType> &val);
    static Matrix<DataType> pw_lin(const Matrix<DataType> &t,
                                   const Matrix<DataType> &tval, const Matrix<DataType> &val);
    static Matrix<DataType> heaviside(const Matrix<DataType> &x);
    static Matrix<DataType> rectangle(const Matrix<DataType> &x);
    static Matrix<DataType> triangle(const Matrix<DataType> &x);
    static Matrix<DataType> ramp(const Matrix<DataType> &x);
    static Matrix<DataType> gauss_quadrature(const Matrix<DataType> &f,
                                             const Matrix<DataType> &x, const Matrix<DataType> &a,
                                             const Matrix<DataType> &b, int order=5);
    static Matrix<DataType> gauss_quadrature(const Matrix<DataType> &f,
                                             const Matrix<DataType> &x, const Matrix<DataType> &a,
                                             const Matrix<DataType> &b, int order,
                                             const Matrix<DataType>& w);
    static Matrix<DataType> jtimes(const Matrix<DataType> &ex, const Matrix<DataType> &arg,
                                   const Matrix<DataType> &v, bool tr=false);
    static std::vector<bool> nl_var(const Matrix<DataType> &expr, const Matrix<DataType> &var);
    static Matrix<DataType> taylor(const Matrix<DataType>& ex, const Matrix<DataType>& x,
                                   const Matrix<DataType>& a, int order);
    static Matrix<DataType> mtaylor(const Matrix<DataType>& ex, const Matrix<DataType>& x,
                                    const Matrix<DataType>& a, int order);
    static Matrix<DataType> mtaylor(const Matrix<DataType>& ex,
                                    const Matrix<DataType>& x, const Matrix<DataType>& a, int order,
                                    const std::vector<int>& order_contributions);
    static Matrix<DataType> poly_coeff(const Matrix<DataType>& ex, const Matrix<DataType>&x);
    static Matrix<DataType> poly_roots(const Matrix<DataType>& p);
    static Matrix<DataType> eig_symbolic(const Matrix<DataType>& m);
    static void qr(const Matrix<DataType>& A, Matrix<DataType>& Q, Matrix<DataType>& R);
    static Matrix<DataType> all(const Matrix<DataType>& x);
    static Matrix<DataType> any(const Matrix<DataType>& x);
    static Matrix<DataType> adj(const Matrix<DataType>& x);
    static Matrix<DataType> getMinor(const Matrix<DataType>& x, int i, int j);
    static Matrix<DataType> cofactor(const Matrix<DataType>& A, int i, int j);
    static Matrix<DataType> chol(const Matrix<DataType>& A);
    static Matrix<DataType> norm_inf_mul(const Matrix<DataType>& x, const Matrix<DataType> &y);
    static Matrix<DataType> diagcat(const std::vector< Matrix<DataType> > &A);
    ///@}
    /// \endcond
#endif // SWIG

    Matrix<DataType> printme(const Matrix<DataType>& y) const;

    /// Transpose the matrix
    Matrix<DataType> T() const;

#if !defined(SWIG) || defined(DOXYGEN)
/**
\ingroup expression_tools
@{
*/
    /** \brief Matrix adjoint
    */
    friend inline Matrix<DataType> adj(const Matrix<DataType>& A) {
      return Matrix<DataType>::adj(A);
    }

    /** \brief Get the (i,j) minor matrix
     */
    friend inline Matrix<DataType> getMinor(const Matrix<DataType> &x, int i, int j) {
      return Matrix<DataType>::getMinor(x, i, j);
    }

    /** \brief Get the (i,j) cofactor matrix
    */
    friend inline Matrix<DataType> cofactor(const Matrix<DataType> &x, int i, int j) {
      return Matrix<DataType>::cofactor(x, i, j);
    }

    /** \brief  QR factorization using the modified Gram-Schmidt algorithm
     * More stable than the classical Gram-Schmidt, but may break down if the rows of A
     * are nearly linearly dependent
     * See J. Demmel: Applied Numerical Linear Algebra (algorithm 3.1.).
     * Note that in SWIG, Q and R are returned by value.
     */
    friend inline void qr(const Matrix<DataType>& A, Matrix<DataType>& Q, Matrix<DataType>& R) {
      return Matrix<DataType>::qr(A, Q, R);
    }

    /** \brief Obtain a Cholesky factorisation of a matrix
     * Returns an upper triangular R such that R'R = A.
     * Matrix A must be positive definite.
     *
     * At the moment, the algorithm is dense (Cholesky-Banachiewicz).
     * There is an open ticket #1212 to make it sparse.
     */
    friend inline Matrix<DataType> chol(const Matrix<DataType>& A) {
      return Matrix<DataType>::chol(A);
    }

    /** \brief Returns true only if any element in the matrix is true
     */
    friend inline Matrix<DataType> any(const Matrix<DataType> &x) {
      return Matrix<DataType>::any(x);
    }

    /** \brief Returns true only if every element in the matrix is true
     */
    friend inline Matrix<DataType> all(const Matrix<DataType> &x) {
      return Matrix<DataType>::all(x);
    }

    /** \brief Inf-norm of a Matrix-Matrix product
    */
    friend inline Matrix<DataType>
      norm_inf_mul(const Matrix<DataType> &x, const Matrix<DataType> &y) {
      return Matrix<DataType>::norm_inf_mul(x, y);
    }

    /** \brief  Make a matrix sparse by removing numerical zeros
    */
    friend inline Matrix<DataType>
      sparsify(const Matrix<DataType>& A, double tol=0) {
      return Matrix<DataType>::sparsify(A, tol);
    }

    /** \brief  Expand the expression as a weighted sum (with constant weights)
     */
    friend inline void expand(const Matrix<DataType>& ex, Matrix<DataType> &weights,
                              Matrix<DataType>& terms) {
      Matrix<DataType>::expand(ex, weights, terms);
    }

    /** \brief Create a piecewise constant function
        Create a piecewise constant function with n=val.size() intervals

        Inputs:
        \param t a scalar variable (e.g. time)
        \param tval vector with the discrete values of t at the interval transitions (length n-1)
        \param val vector with the value of the function for each interval (length n)
    */
    friend inline Matrix<DataType> pw_const(const Matrix<DataType> &t,
                                            const Matrix<DataType> &tval,
                                            const Matrix<DataType> &val) {
      return Matrix<DataType>::pw_const(t, tval, val);
    }

    /** Create a piecewise linear function
        Create a piecewise linear function:

        Inputs:
        \brief t a scalar variable (e.g. time)
        \brief tval vector with the the discrete values of t (monotonically increasing)
        \brief val vector with the corresponding function values (same length as tval)
    */
    friend inline Matrix<DataType>
      pw_lin(const Matrix<DataType> &t, const Matrix<DataType> &tval,
             const Matrix<DataType> &val) {
      return Matrix<DataType>::pw_lin(t, tval, val);
    }

    /**  \brief Heaviside function
     *
     * \f[
     * \begin {cases}
     * H(x) = 0 & x<0 \\
     * H(x) = 1/2 & x=0 \\
     * H(x) = 1 & x>0 \\
     * \end {cases}
     * \f]
     */
    friend inline Matrix<DataType> heaviside(const Matrix<DataType> &x) {
      return Matrix<DataType>::heaviside(x);
    }

    /**
     * \brief rectangle function
     *
     * \f[
     * \begin {cases}
     * \Pi(x) = 1     & |x| < 1/2 \\
     * \Pi(x) = 1/2   & |x| = 1/2  \\
     * \Pi(x) = 0     & |x| > 1/2  \\
     * \end {cases}
     * \f]
     *
     * Also called: gate function, block function, band function, pulse function, window function
     */
    friend inline Matrix<DataType> rectangle(const Matrix<DataType> &x) {
      return Matrix<DataType>::rectangle(x);
    }

    /**
     * \brief triangle function
     *
     * \f[
     * \begin {cases}
     * \Lambda(x) = 0 &    |x| >= 1  \\
     * \Lambda(x) = 1-|x| &  |x| < 1
     * \end {cases}
     * \f]
     *
     */
    friend inline Matrix<DataType> triangle(const Matrix<DataType> &x) {
      return Matrix<DataType>::triangle(x);
    }

    /**
     * \brief ramp function
     *
     *
     * \f[
     * \begin {cases}
     *  R(x) = 0   & x <= 1 \\
     *  R(x) = x   & x > 1 \\
     * \end {cases}
     * \f]
     *
     * Also called: slope function
     */
    friend inline Matrix<DataType> ramp(const Matrix<DataType> &x) {
      return Matrix<DataType>::ramp(x);
    }

    ///@{
    /** \brief  Integrate f from a to b using Gaussian quadrature with n points */
    friend inline Matrix<DataType>
      gauss_quadrature(const Matrix<DataType> &f, const Matrix<DataType> &x,
                       const Matrix<DataType> &a, const Matrix<DataType> &b,
                       int order=5) {
      return Matrix<DataType>::gauss_quadrature(f, x, a, b, order);
    }
    friend inline Matrix<DataType>
      gauss_quadrature(const Matrix<DataType> &f, const Matrix<DataType> &x,
                       const Matrix<DataType> &a, const Matrix<DataType> &b,
                       int order, const Matrix<DataType>& w) {
      return Matrix<DataType>::gauss_quadrature(f, x, a, b, order, w);
    }
    ///@}

    ///@{
    /**
     * \brief univariate Taylor series expansion
     *
     * Calculate the Taylor expansion of expression 'ex' up to order 'order' with
     * respect to variable 'x' around the point 'a'
     *
     * \f$(x)=f(a)+f'(a)(x-a)+f''(a)\frac {(x-a)^2}{2!}+f'''(a)\frac{(x-a)^3}{3!}+\ldots\f$
     *
     * Example usage:
     * \code
     * taylor(sin(x), x)
     * \endcode
     * \verbatim >>   x \endverbatim
     */
    friend inline Matrix<DataType> taylor(const Matrix<DataType>& ex, const Matrix<DataType>& x,
                                          const Matrix<DataType>& a, int order=1) {
      return Matrix<DataType>::taylor(ex, x, a, order);
    }
    friend inline Matrix<DataType> taylor(const Matrix<DataType>& ex, const Matrix<DataType>& x) {
      return Matrix<DataType>::taylor(ex, x, 0, 1);
    }
    ///@}

    /**
     * \brief multivariate Taylor series expansion
     *
     * Do Taylor expansions until the aggregated order of a term is equal to 'order'.
     * The aggregated order of \f$x^n y^m\f$ equals \f$n+m\f$.
     *
     */
    friend inline Matrix<DataType> mtaylor(const Matrix<DataType>& ex, const Matrix<DataType>& x,
                                           const Matrix<DataType>& a, int order=1) {
      return Matrix<DataType>::mtaylor(ex, x, a, order);
    }

    /**
     * \brief multivariate Taylor series expansion
     *
     * Do Taylor expansions until the aggregated order of a term is equal to 'order'.
     * The aggregated order of \f$x^n y^m\f$ equals \f$n+m\f$.
     *
     * The argument order_contributions can denote how match each variable contributes
     * to the aggregated order. If x=[x, y] and order_contributions=[1, 2], then the
     * aggregated order of \f$x^n y^m\f$ equals \f$1n+2m\f$.
     *
     * Example usage
     *
     * \code
     * taylor(sin(x+y),[x, y],[a, b], 1)
     * \endcode
     * \f$ \sin(b+a)+\cos(b+a)(x-a)+\cos(b+a)(y-b) \f$
     * \code
     * taylor(sin(x+y),[x, y],[0, 0], 4)
     * \endcode
     * \f$  y+x-(x^3+3y x^2+3 y^2 x+y^3)/6  \f$
     * \code
     * taylor(sin(x+y),[x, y],[0, 0], 4,[1, 2])
     * \endcode
     * \f$  (-3 x^2 y-x^3)/6+y+x \f$
     *
     */
    friend inline Matrix<DataType> mtaylor(const Matrix<DataType>& ex, const Matrix<DataType>& x,
                                           const Matrix<DataType>& a, int order,
                                           const std::vector<int>& order_contributions) {
      return Matrix<DataType>::mtaylor(ex, x, a, order, order_contributions);
    }

    /** \brief extracts polynomial coefficients from an expression
     *
     * \param ex Scalar expression that represents a polynomial
     * \param x  Scalar symbol that the polynomial is build up with
     */
    friend inline Matrix<DataType> poly_coeff(const Matrix<DataType>& f,
                                              const Matrix<DataType>& x) {
      return Matrix<DataType>::poly_coeff(f, x);
    }

    /** \brief Attempts to find the roots of a polynomial
     *
     *  This will only work for polynomials up to order 3
     *  It is assumed that the roots are real.
     *
     */
    friend inline Matrix<DataType> poly_roots(const Matrix<DataType>& p) {
      return Matrix<DataType>::poly_roots(p);
    }

    /** \brief Attempts to find the eigenvalues of a symbolic matrix
     *  This will only work for up to 3x3 matrices
     */
    friend inline Matrix<DataType> eig_symbolic(const Matrix<DataType>& m) {
      return Matrix<DataType>::eig_symbolic(m);
    }
/** @} */
#endif

    /** \brief Set or reset the depth to which equalities are being checked for simplifications */
    static void setEqualityCheckingDepth(int eq_depth=1);

    /** \brief Get the depth to which equalities are being checked for simplifications */
    static int getEqualityCheckingDepth();

    /** \brief Get function input */
    static std::vector<Matrix<DataType> > get_input(const Function& f);

    ///@{
    /** \brief Jacobian expression */
    static Matrix<DataType> jac(const Function& f, int iind=0, int oind=0,
                  bool compact=false, bool symmetric=false);
    static Matrix<DataType> jac(const Function& f, const std::string & iname, int oind=0,
                  bool compact=false, bool symmetric=false);
    static Matrix<DataType> jac(const Function& f, int iind, const std::string& oname,
                  bool compact=false, bool symmetric=false);
    static Matrix<DataType> jac(const Function& f, const std::string& iname,
                                const std::string& oname,
                                bool compact=false, bool symmetric=false);
    ///@}

    ///@{
    /** \brief Gradient expression */
    static Matrix<DataType> grad(const Function& f, int iind=0, int oind=0);
    static Matrix<DataType> grad(const Function& f, const std::string& iname, int oind=0);
    static Matrix<DataType> grad(const Function& f, int iind, const std::string& oname);
    static Matrix<DataType> grad(const Function& f, const std::string& iname,
                                 const std::string& oname);
    ///@}

    ///@{
    /** \brief Tangent expression */
    static Matrix<DataType> tang(const Function& f, int iind=0, int oind=0);
    static Matrix<DataType> tang(const Function& f, const std::string& iname, int oind=0);
    static Matrix<DataType> tang(const Function& f, int iind, const std::string& oname);
    static Matrix<DataType> tang(const Function& f, const std::string& iname,
                                 const std::string& oname);
    ///@}

    ///@{
    /** \bried Hessian expression */
    static Matrix<DataType> hess(const Function& f, int iind=0, int oind=0);
    static Matrix<DataType> hess(const Function& f, const std::string& iname, int oind=0);
    static Matrix<DataType> hess(const Function& f, int iind, const std::string& oname);
    static Matrix<DataType> hess(const Function& f, const std::string& iname,
                                 const std::string& oname);
    ///@}

    /// Get name of the class
    static std::string className();

    /// Print a description of the object
    void print(std::ostream &stream=casadi::userOut(), bool trailing_newline=true) const;

    /// Get strings corresponding to the nonzeros and the interdependencies
    void printSplit(std::vector<std::string>& SWIG_OUTPUT(nz),
                    std::vector<std::string>& SWIG_OUTPUT(inter)) const;

    /// Print a representation of the object
    void repr(std::ostream &stream=casadi::userOut(), bool trailing_newline=true) const;

    /// Print scalar
    void printScalar(std::ostream &stream=casadi::userOut(), bool trailing_newline=true) const;

    /// Print vector-style
    void printVector(std::ostream &stream=casadi::userOut(), bool trailing_newline=true) const;

    /// Print dense matrix-stype
    void printDense(std::ostream &stream=casadi::userOut(), bool trailing_newline=true) const;

    /// Print sparse matrix style
    void printSparse(std::ostream &stream=casadi::userOut(), bool trailing_newline=true) const;

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

#ifndef SWIG
    /// Access the non-zero elements
    std::vector<DataType>& data();

    /// Const access the non-zero elements
    const std::vector<DataType>& data() const;

    /// \cond INTERNAL
    /// Get a pointer to the data
    DataType* ptr() {
      return is_empty() ? static_cast<DataType*>(0) : &data_.front();
    }
    friend inline DataType* getPtr(Matrix<DataType>& v) {
      return v.ptr();
    }

    /// Get a const pointer to the data
    const DataType* ptr() const {
      return is_empty() ? static_cast<const DataType*>(0) : &data_.front();
    }
    friend inline const DataType* getPtr(const Matrix<DataType>& v) {
      return v.ptr();
    }
    /// \endcond

    /// Const access the sparsity - reference to data member
    const Sparsity& sparsity() const { return sparsity_; }

#endif // SWIG

    /** \brief Get an owning reference to the sparsity pattern */
    Sparsity getSparsity() const { return sparsity();}

#ifndef SWIG
    /// \cond INTERNAL
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

    /** \brief  Save the result to the LAPACK banded format -- see LAPACK documentation
        kl:    The number of subdiagonals in res
        ku:    The number of superdiagonals in res
        ldres: The leading dimension in res
        res:   The number of superdiagonals */
    void getBand(int kl, int ku, int ldres, DataType *res) const;
/// \endcond
#endif

    /* \brief Construct a sparse matrix from triplet form
     * Default matrix size is max(col) x max(row)
     */
    ///@{
    static Matrix<DataType> triplet(const std::vector<int>& row, const std::vector<int>& col,
                                    const Matrix<DataType>& d);
    static Matrix<DataType> triplet(const std::vector<int>& row, const std::vector<int>& col,
                                    const Matrix<DataType>& d, int nrow, int ncol);
    static Matrix<DataType> triplet(const std::vector<int>& row, const std::vector<int>& col,
                                    const Matrix<DataType>& d, const std::pair<int, int>& rc);
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
    size_t element_hash() const;

    /// Checks if expression does not contain NaN or Inf
    bool is_regular() const;

    /** \brief Check if smooth */
    bool is_smooth() const;

    /** \brief Check if SX is a leaf of the SX graph

        Only defined if symbolic scalar.
    */
    bool is_leaf() const;

    /** \brief Check whether a binary SX is commutative

        Only defined if symbolic scalar.
    */
    bool is_commutative() const;

    /** \brief Check if symbolic (Dense)
        Sparse matrices invariable return false
    */
    bool is_symbolic() const;

    /** \brief Check if matrix can be used to define function inputs.
        Sparse matrices can return true if all non-zero elements are symbolic
    */
    bool is_valid_input() const;

    /// \cond INTERNAL
    /** \brief Detect duplicate symbolic expressions
        If there are symbolic primitives appearing more than once, the function will return
        true and the names of the duplicate expressions will be printed to userOut<true, PL_WARN>().
        Note: Will mark the node using SXElem::setTemp.
        Make sure to call resetInput() after usage.
    */
    bool has_duplicates();

    /** \brief Reset the marker for an input expression */
    void resetInput();
  /// \endcond

    /** \brief Check if the matrix is constant (note that false negative answers are possible)*/
    bool is_constant() const;

    /** \brief Check if the matrix is integer-valued
     * (note that false negative answers are possible)*/
    bool is_integer() const;

    /** \brief  check if the matrix is 0 (note that false negative answers are possible)*/
    bool is_zero() const;

    /** \brief  check if the matrix is 1 (note that false negative answers are possible)*/
    bool is_one() const;

    /** \brief  check if the matrix is -1 (note that false negative answers are possible)*/
    bool is_minus_one() const;

    /** \brief  check if the matrix is an identity matrix (note that false negative answers
     * are possible)*/
    bool is_identity() const;

    /** \brief  Check if the matrix has any zero entries which are not structural zeros */
    bool has_zeros() const;

    /** \brief Get double value (only if constant) */
    double getValue() const;

    /** \brief Get double value (particular nonzero) */
    double getValue(int k) const;

    /** \brief Set double value (only if constant) */
    void setValue(double m);

    /** \brief Set double value (particular nonzero) */
    void setValue(double m, int k);

    /** \brief Get double value (only if integer constant) */
    int getIntValue() const;

    /** \brief Get all nonzeros */
    std::vector<double> nonzeros() const;

    /** \brief Get all nonzeros */
    std::vector<int> nonzeros_int() const;

#ifndef SWIG
    /** \brief Type conversion to double */
    explicit operator double() const;

    /** \brief Type conversion to double vector */
    explicit operator std::vector<double>() const;

    /** \brief Type conversion to int */
    explicit operator int() const;
#endif // SWIG

    /** \brief Get name (only if symbolic scalar) */
    std::string getName() const;

    /** \brief Get expressions of the children of the expression
        Only defined if symbolic scalar.
        Wraps SXElem SXElem::dep(int ch=0) const.
     */
    Matrix<DataType> dep(int ch=0) const;

    /** \brief Get the number of dependencies of a binary SXElem
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

  ///@{
  /// Readability typedefs
  typedef Matrix<int> IM;
  typedef Matrix<double> DM;
  typedef std::vector<DM> DMVector;
  typedef std::vector<DMVector> DMVectorVector;
  typedef std::map<std::string, DM> DMDict;
  ///@}
} // namespace casadi

#endif // CASADI_MATRIX_HPP
