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
  template <typename Scalar> inline std::string matrixName()
  { return std::string("Matrix<") + typeid(Scalar).name() + std::string(">");}
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

      Matrix<Scalar> is polymorphic with a std::vector<Scalar> that contain
      all non-identical-zero elements.\n
      The sparsity can be accessed with Sparsity& sparsity()\n

      \author Joel Andersson
      \date 2010-2014
  */
  template<typename Scalar>
  class CASADI_EXPORT Matrix :
    public MatrixCommon,
    public GenericExpression<Matrix<Scalar> >,
    public GenericMatrix<Matrix<Scalar> >,
    public PrintableObject<Matrix<Scalar> > {
  public:

    /** \brief  constructors */
    /// empty 0-by-0 matrix constructor
    Matrix();

    /// Copy constructor
    Matrix(const Matrix<Scalar>& m);

#ifndef SWIG
    /// Assignment (normal)
    Matrix<Scalar>& operator=(const Matrix<Scalar>& m);
#endif // SWIG

    /** \brief Create a sparse matrix with all structural zeros */
    Matrix(int nrow, int ncol);

#ifndef SWIG
    /** \brief Create a sparse matrix with all structural zeros */
    explicit Matrix(const std::pair<int, int>& rc);

    /** \brief  Access functions of the node */
    std::vector<Scalar>* operator->() { return &nonzeros_;}

    /** \brief  Const access functions of the node */
    const std::vector<Scalar>* operator->() const { return &nonzeros_;}
#endif // SWIG

    /** \brief Create a sparse matrix from a sparsity pattern.
        Same as Matrix::ones(sparsity)
     */
    explicit Matrix(const Sparsity& sp);

    /** \brief Construct matrix with a given sparsity and nonzeros */
    Matrix(const Sparsity& sp, const Matrix<Scalar>& d);

    /** \brief Check if the dimensions and colind, row vectors are compatible.
     * \param complete  set to true to also check elementwise
     * throws an error as possible result
     */
    void sanity_check(bool complete=false) const;

    /// This constructor enables implicit type conversion from a numeric type
    Matrix(double val);

#if !(defined(SWIG) && defined(SWIGMATLAB))
    /// Dense matrix constructor with data given as vector of vectors
    explicit Matrix(const std::vector< std::vector<double> >& m);

    /** \brief  Create an expression from a vector  */
    template<typename A>
    Matrix(const std::vector<A>& x) : sparsity_(Sparsity::dense(x.size(), 1)),
      nonzeros_(std::vector<Scalar>(x.size())) {
        auto x_it = x.begin();
        for (auto&& d : nonzeros_) d = static_cast<Scalar>(*x_it++);
    }
#endif

    /** \brief Create a matrix from another matrix with a different entry type
     *  Assumes that the scalar conversion is valid.
     */
    template<typename A>
    Matrix(const Matrix<A>& x) : sparsity_(x.sparsity()), nonzeros_(std::vector<Scalar>(x.nnz())) {
      auto x_it = x->begin();
      for (auto&& d : nonzeros_) d = static_cast<Scalar>(*x_it++);
    }

#ifndef SWIG
    /// Construct from a vector
    Matrix(const std::vector<Scalar>& x);

    /// Convert to scalar type
    const Scalar scalar() const;

    /// Scalar type
    typedef Scalar ScalarType;

    /// Base class
    typedef GenericMatrix<Matrix<Scalar> > B;

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
    using B::nz;
    using B::operator();
    using B::horzsplit;
    using B::vertsplit;
    using B::diagsplit;
    using B::mtimes;
#endif // SWIG

    /// Returns true if the matrix has a non-zero at location rr, cc
    bool has_nz(int rr, int cc) const { return sparsity().has_nz(rr, cc); }

    /// Returns the truth value of a Matrix
    bool __nonzero__() const;

    ///@{
    /// Get a submatrix, single argument
    void get(Matrix<Scalar>& SWIG_OUTPUT(m), bool ind1, const Slice& rr) const;
    void get(Matrix<Scalar>& SWIG_OUTPUT(m), bool ind1, const Matrix<int>& rr) const;
    void get(Matrix<Scalar>& SWIG_OUTPUT(m), bool ind1, const Sparsity& sp) const;
    ///@}

    /// Get a submatrix, two arguments
    ///@{
    void get(Matrix<Scalar>& SWIG_OUTPUT(m), bool ind1,
                const Slice& rr, const Slice& cc) const;
    void get(Matrix<Scalar>& SWIG_OUTPUT(m), bool ind1,
                const Slice& rr, const Matrix<int>& cc) const;
    void get(Matrix<Scalar>& SWIG_OUTPUT(m), bool ind1,
                const Matrix<int>& rr, const Slice& cc) const;
    void get(Matrix<Scalar>& SWIG_OUTPUT(m), bool ind1,
                const Matrix<int>& rr, const Matrix<int>& cc) const;
    ///@}

    ///@{
    /// Set a submatrix, single argument
    void set(const Matrix<Scalar>& m, bool ind1, const Slice& rr);
    void set(const Matrix<Scalar>& m, bool ind1, const Matrix<int>& rr);
    void set(const Matrix<Scalar>& m, bool ind1, const Sparsity& sp);
    ///@}

    ///@{
    /// Set a submatrix, two arguments
    void set(const Matrix<Scalar>& m, bool ind1, const Slice& rr, const Slice& cc);
    void set(const Matrix<Scalar>& m, bool ind1, const Slice& rr, const Matrix<int>& cc);
    void set(const Matrix<Scalar>& m, bool ind1, const Matrix<int>& rr, const Slice& cc);
    void set(const Matrix<Scalar>& m, bool ind1, const Matrix<int>& rr, const Matrix<int>& cc);
    ///@}

    ///@{
    /// Get a set of nonzeros
    void get_nz(Matrix<Scalar>& SWIG_OUTPUT(m), bool ind1, const Slice& k) const;
    void get_nz(Matrix<Scalar>& SWIG_OUTPUT(m), bool ind1, const Matrix<int>& k) const;
    ///@}

    ///@{
    /// Set a set of nonzeros
    void set_nz(const Matrix<Scalar>& m, bool ind1, const Slice& k);
    void set_nz(const Matrix<Scalar>& m, bool ind1, const Matrix<int>& k);
    ///@}

    Matrix<Scalar> operator+() const;
    Matrix<Scalar> operator-() const;

    /// \cond INTERNAL
    ///@{
    /** \brief  Create nodes by their ID */
    static Matrix<Scalar> binary(int op, const Matrix<Scalar> &x, const Matrix<Scalar> &y);
    static Matrix<Scalar> unary(int op, const Matrix<Scalar> &x);
    static Matrix<Scalar> scalar_matrix(int op,
                                          const Matrix<Scalar> &x, const Matrix<Scalar> &y);
    static Matrix<Scalar> matrix_scalar(int op,
                                          const Matrix<Scalar> &x, const Matrix<Scalar> &y);
    static Matrix<Scalar> matrix_matrix(int op,
                                          const Matrix<Scalar> &x, const Matrix<Scalar> &y);
    ///@}
    /// \endcond

#ifndef SWIG
    /// \cond CLUTTER
    ///@{
    /// Functions called by friend functions defined for GenericExpression
    static bool is_equal(const Matrix<Scalar> &x, const Matrix<Scalar> &y, int depth=0);
    ///@}

    ///@{
    /// Functions called by friend functions defined for GenericMatrix
    static Matrix<Scalar> simplify(const Matrix<Scalar> &x);
    static Matrix<Scalar> jacobian(const Matrix<Scalar> &f, const Matrix<Scalar> &x,
                                     bool symmetric=false);
    static Matrix<Scalar> gradient(const Matrix<Scalar> &f, const Matrix<Scalar> &x);
    static Matrix<Scalar> tangent(const Matrix<Scalar> &f, const Matrix<Scalar> &x);
    static Matrix<Scalar> hessian(const Matrix<Scalar> &f, const Matrix<Scalar> &x);
    static Matrix<Scalar> hessian(const Matrix<Scalar> &f, const Matrix<Scalar> &x,
                                    Matrix<Scalar>& g);
    static Matrix<Scalar>
      substitute(const Matrix<Scalar>& ex,
                 const Matrix<Scalar>& v,
                 const Matrix<Scalar>& vdef);
    static std::vector<Matrix<Scalar> >
      substitute(const std::vector<Matrix<Scalar> >& ex,
                 const std::vector<Matrix<Scalar> >& v,
                 const std::vector<Matrix<Scalar> >& vdef);
    static void substitute_inplace(const std::vector<Matrix<Scalar> >& v,
                                  std::vector<Matrix<Scalar> >& vdef,
                                  std::vector<Matrix<Scalar> >& ex,
                                  bool revers);
    static Matrix<Scalar> pinv(const Matrix<Scalar> &x);
    static Matrix<Scalar> pinv(const Matrix<Scalar> &A,
                                 const std::string& lsolver, const Dict& opts);
    static Matrix<Scalar> solve(const Matrix<Scalar> &A, const Matrix<Scalar>& b);
    static Matrix<Scalar> solve(const Matrix<Scalar> &A, const Matrix<Scalar>& b,
                                  const std::string& lsolver, const Dict& opts);
    static int n_nodes(const Matrix<Scalar> &x);
    static std::string print_operator(const Matrix<Scalar> &x,
                                      const std::vector<std::string>& args);
    static void shared(std::vector<Matrix<Scalar> >& ex,
                              std::vector<Matrix<Scalar> >& v,
                              std::vector<Matrix<Scalar> >& vdef,
                              const std::string& v_prefix,
                              const std::string& v_suffix);
    static Matrix<Scalar> _bilin(const Matrix<Scalar>& A,
                                   const Matrix<Scalar>& x,
                                   const Matrix<Scalar>& y);
    static Matrix<Scalar> _rank1(const Matrix<Scalar>& A,
                                   const Matrix<Scalar>& alpha,
                                   const Matrix<Scalar>& x,
                                   const Matrix<Scalar>& y);
    static Matrix<Scalar> if_else(const Matrix<Scalar> &x,
                                    const Matrix<Scalar> &if_true,
                                    const Matrix<Scalar> &if_false,
                                    bool short_circuit);
    static Matrix<Scalar> conditional(const Matrix<Scalar> &ind,
                                        const std::vector<Matrix<Scalar> > &x,
                                        const Matrix<Scalar> &x_default,
                                        bool short_circuit);
    static bool depends_on(const Matrix<Scalar> &x, const Matrix<Scalar> &arg);
    static Matrix<Scalar> mpower(const Matrix<Scalar> &x, const Matrix<Scalar> &y);
    static Matrix<Scalar> mrdivide(const Matrix<Scalar> &x, const Matrix<Scalar> &y);
    static Matrix<Scalar> mldivide(const Matrix<Scalar> &x, const Matrix<Scalar> &y);
    static std::vector<Matrix<Scalar> > symvar(const Matrix<Scalar> &x);
    static Matrix<Scalar> det(const Matrix<Scalar> &x);
    static Matrix<Scalar> inv(const Matrix<Scalar> &x);
    static Matrix<Scalar> trace(const Matrix<Scalar> &x);
    static Matrix<Scalar> norm_1(const Matrix<Scalar> &x);
    static Matrix<Scalar> norm_2(const Matrix<Scalar> &x);
    static Matrix<Scalar> norm_fro(const Matrix<Scalar> &x);
    static Matrix<Scalar> norm_inf(const Matrix<Scalar> &x);
    static Matrix<Scalar> sum2(const Matrix<Scalar> &x);
    static Matrix<Scalar> sum1(const Matrix<Scalar> &x);
    static Matrix<Scalar> dot(const Matrix<Scalar> &x, const Matrix<Scalar> &y);
    static Matrix<Scalar> nullspace(const Matrix<Scalar> &x);
    static Matrix<Scalar> diag(const Matrix<Scalar> &x);
    static Matrix<Scalar> unite(const Matrix<Scalar> &A, const Matrix<Scalar>& B);
    static Matrix<Scalar> project(const Matrix<Scalar> &x,
                                    const Sparsity& sp, bool intersect=false);
    static Matrix<Scalar> polyval(const Matrix<Scalar> &p, const Matrix<Scalar>& x);
    static Matrix<Scalar> densify(const Matrix<Scalar> &x, const Matrix<Scalar>& val);
    static Matrix<Scalar> densify(const Matrix<Scalar> &x);
    ///@}

    ///@{
    /// Functions called by friend functions defined for SparsityInterface
    static Matrix<Scalar> blockcat(const std::vector< std::vector<Matrix<Scalar> > > &v);
    static Matrix<Scalar> horzcat(const std::vector<Matrix<Scalar> > &v);
    static std::vector<Matrix<Scalar> >
      horzsplit(const Matrix<Scalar>& x,
                const std::vector<int>& offset);
    static Matrix<Scalar> vertcat(const std::vector<Matrix<Scalar> > &v);
    static std::vector< Matrix<Scalar> >
      vertsplit(const Matrix<Scalar>& x,
                const std::vector<int>& offset);
    static std::vector< Matrix<Scalar> >
      diagsplit(const Matrix<Scalar>& x,
                const std::vector<int>& offset1,
                const std::vector<int>& offset2);
    static Matrix<Scalar> reshape(const Matrix<Scalar> &x, int nrow, int ncol);
    static Matrix<Scalar> reshape(const Matrix<Scalar> &x, const Sparsity& sp);
    static Matrix<Scalar> kron(const Matrix<Scalar> &x, const Matrix<Scalar>& y);
    static Matrix<Scalar> mtimes(const Matrix<Scalar> &x, const Matrix<Scalar> &y);
    static Matrix<Scalar> mac(const Matrix<Scalar> &x,
                                const Matrix<Scalar> &y,
                                const Matrix<Scalar> &z);
    ///@}

    ///@{
    /// Functions called by friend functions defined here
    static Matrix<Scalar> sparsify(const Matrix<Scalar> &x, double tol=0);
    static void expand(const Matrix<Scalar>& x,
                       Matrix<Scalar>& weights,
                       Matrix<Scalar>& terms);
    static Matrix<Scalar> pw_const(const Matrix<Scalar> &t,
                                     const Matrix<Scalar> &tval, const Matrix<Scalar> &val);
    static Matrix<Scalar> pw_lin(const Matrix<Scalar> &t,
                                   const Matrix<Scalar> &tval, const Matrix<Scalar> &val);
    static Matrix<Scalar> heaviside(const Matrix<Scalar> &x);
    static Matrix<Scalar> rectangle(const Matrix<Scalar> &x);
    static Matrix<Scalar> triangle(const Matrix<Scalar> &x);
    static Matrix<Scalar> ramp(const Matrix<Scalar> &x);
    static Matrix<Scalar> gauss_quadrature(const Matrix<Scalar> &f,
                                             const Matrix<Scalar> &x, const Matrix<Scalar> &a,
                                             const Matrix<Scalar> &b, int order=5);
    static Matrix<Scalar> gauss_quadrature(const Matrix<Scalar> &f,
                                             const Matrix<Scalar> &x, const Matrix<Scalar> &a,
                                             const Matrix<Scalar> &b, int order,
                                             const Matrix<Scalar>& w);
    static Matrix<Scalar> jtimes(const Matrix<Scalar> &ex, const Matrix<Scalar> &arg,
                                   const Matrix<Scalar> &v, bool tr=false);

#ifdef WITH_DEPRECATED_FEATURES
    static std::vector<bool> nl_var(const Matrix<Scalar> &expr, const Matrix<Scalar> &var);
#endif

    static std::vector<std::vector<Matrix<Scalar> > >
    forward(const std::vector<Matrix<Scalar> > &ex,
            const std::vector<Matrix<Scalar> > &arg,
            const std::vector<std::vector<Matrix<Scalar> > > &v,
            const Dict& opts = Dict());
    static std::vector<std::vector<Matrix<Scalar> > >
    reverse(const std::vector<Matrix<Scalar> > &ex,
            const std::vector<Matrix<Scalar> > &arg,
            const std::vector<std::vector<Matrix<Scalar> > > &v,
            const Dict& opts = Dict());

    static std::vector<bool> which_depends(const Matrix<Scalar> &expr, const Matrix<Scalar> &var,
        int order=1, bool tr=false);
    static Matrix<Scalar> taylor(const Matrix<Scalar>& ex, const Matrix<Scalar>& x,
                                   const Matrix<Scalar>& a, int order);
    static Matrix<Scalar> mtaylor(const Matrix<Scalar>& ex, const Matrix<Scalar>& x,
                                    const Matrix<Scalar>& a, int order);
    static Matrix<Scalar> mtaylor(const Matrix<Scalar>& ex,
                                    const Matrix<Scalar>& x, const Matrix<Scalar>& a, int order,
                                    const std::vector<int>& order_contributions);
    static Matrix<Scalar> poly_coeff(const Matrix<Scalar>& ex, const Matrix<Scalar>&x);
    static Matrix<Scalar> poly_roots(const Matrix<Scalar>& p);
    static Matrix<Scalar> eig_symbolic(const Matrix<Scalar>& m);
    static void qr(const Matrix<Scalar>& A, Matrix<Scalar>& Q, Matrix<Scalar>& R);
    static Matrix<Scalar> all(const Matrix<Scalar>& x);
    static Matrix<Scalar> any(const Matrix<Scalar>& x);
    static Matrix<Scalar> adj(const Matrix<Scalar>& x);
    static Matrix<Scalar> getMinor(const Matrix<Scalar>& x, int i, int j);
    static Matrix<Scalar> cofactor(const Matrix<Scalar>& A, int i, int j);
    static Matrix<Scalar> chol(const Matrix<Scalar>& A);
    static Matrix<Scalar> norm_inf_mul(const Matrix<Scalar>& x, const Matrix<Scalar> &y);
    static Matrix<Scalar> diagcat(const std::vector< Matrix<Scalar> > &A);
    ///@}
    /// \endcond
#endif // SWIG

    Matrix<Scalar> printme(const Matrix<Scalar>& y) const;

    /// Transpose the matrix
    Matrix<Scalar> T() const;

#if !defined(SWIG) || defined(DOXYGEN)
/**
\ingroup expression_tools
@{
*/
    /** \brief Matrix adjoint
    */
    friend inline Matrix<Scalar> adj(const Matrix<Scalar>& A) {
      return Matrix<Scalar>::adj(A);
    }

    /** \brief Get the (i,j) minor matrix
     */
    friend inline Matrix<Scalar> getMinor(const Matrix<Scalar> &x, int i, int j) {
      return Matrix<Scalar>::getMinor(x, i, j);
    }

    /** \brief Get the (i,j) cofactor matrix
    */
    friend inline Matrix<Scalar> cofactor(const Matrix<Scalar> &x, int i, int j) {
      return Matrix<Scalar>::cofactor(x, i, j);
    }

    /** \brief  QR factorization using the modified Gram-Schmidt algorithm
     * More stable than the classical Gram-Schmidt, but may break down if the rows of A
     * are nearly linearly dependent
     * See J. Demmel: Applied Numerical Linear Algebra (algorithm 3.1.).
     * Note that in SWIG, Q and R are returned by value.
     */
    friend inline void qr(const Matrix<Scalar>& A, Matrix<Scalar>& Q, Matrix<Scalar>& R) {
      return Matrix<Scalar>::qr(A, Q, R);
    }

    /** \brief Obtain a Cholesky factorisation of a matrix
     * Returns an upper triangular R such that R'R = A.
     * Matrix A must be positive definite.
     *
     * At the moment, the algorithm is dense (Cholesky-Banachiewicz).
     * There is an open ticket #1212 to make it sparse.
     */
    friend inline Matrix<Scalar> chol(const Matrix<Scalar>& A) {
      return Matrix<Scalar>::chol(A);
    }

    /** \brief Returns true only if any element in the matrix is true
     */
    friend inline Matrix<Scalar> any(const Matrix<Scalar> &x) {
      return Matrix<Scalar>::any(x);
    }

    /** \brief Returns true only if every element in the matrix is true
     */
    friend inline Matrix<Scalar> all(const Matrix<Scalar> &x) {
      return Matrix<Scalar>::all(x);
    }

    /** \brief Inf-norm of a Matrix-Matrix product
    */
    friend inline Matrix<Scalar>
      norm_inf_mul(const Matrix<Scalar> &x, const Matrix<Scalar> &y) {
      return Matrix<Scalar>::norm_inf_mul(x, y);
    }

    /** \brief  Make a matrix sparse by removing numerical zeros
    */
    friend inline Matrix<Scalar>
      sparsify(const Matrix<Scalar>& A, double tol=0) {
      return Matrix<Scalar>::sparsify(A, tol);
    }

    /** \brief  Expand the expression as a weighted sum (with constant weights)
     */
    friend inline void expand(const Matrix<Scalar>& ex, Matrix<Scalar> &weights,
                              Matrix<Scalar>& terms) {
      Matrix<Scalar>::expand(ex, weights, terms);
    }

    /** \brief Create a piecewise constant function
        Create a piecewise constant function with n=val.size() intervals

        Inputs:
        \param t a scalar variable (e.g. time)
        \param tval vector with the discrete values of t at the interval transitions (length n-1)
        \param val vector with the value of the function for each interval (length n)
    */
    friend inline Matrix<Scalar> pw_const(const Matrix<Scalar> &t,
                                            const Matrix<Scalar> &tval,
                                            const Matrix<Scalar> &val) {
      return Matrix<Scalar>::pw_const(t, tval, val);
    }

    /** Create a piecewise linear function
        Create a piecewise linear function:

        Inputs:
        \brief t a scalar variable (e.g. time)
        \brief tval vector with the the discrete values of t (monotonically increasing)
        \brief val vector with the corresponding function values (same length as tval)
    */
    friend inline Matrix<Scalar>
      pw_lin(const Matrix<Scalar> &t, const Matrix<Scalar> &tval,
             const Matrix<Scalar> &val) {
      return Matrix<Scalar>::pw_lin(t, tval, val);
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
    friend inline Matrix<Scalar> heaviside(const Matrix<Scalar> &x) {
      return Matrix<Scalar>::heaviside(x);
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
    friend inline Matrix<Scalar> rectangle(const Matrix<Scalar> &x) {
      return Matrix<Scalar>::rectangle(x);
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
    friend inline Matrix<Scalar> triangle(const Matrix<Scalar> &x) {
      return Matrix<Scalar>::triangle(x);
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
    friend inline Matrix<Scalar> ramp(const Matrix<Scalar> &x) {
      return Matrix<Scalar>::ramp(x);
    }

    ///@{
    /** \brief  Integrate f from a to b using Gaussian quadrature with n points */
    friend inline Matrix<Scalar>
      gauss_quadrature(const Matrix<Scalar> &f, const Matrix<Scalar> &x,
                       const Matrix<Scalar> &a, const Matrix<Scalar> &b,
                       int order=5) {
      return Matrix<Scalar>::gauss_quadrature(f, x, a, b, order);
    }
    friend inline Matrix<Scalar>
      gauss_quadrature(const Matrix<Scalar> &f, const Matrix<Scalar> &x,
                       const Matrix<Scalar> &a, const Matrix<Scalar> &b,
                       int order, const Matrix<Scalar>& w) {
      return Matrix<Scalar>::gauss_quadrature(f, x, a, b, order, w);
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
    friend inline Matrix<Scalar> taylor(const Matrix<Scalar>& ex, const Matrix<Scalar>& x,
                                          const Matrix<Scalar>& a, int order=1) {
      return Matrix<Scalar>::taylor(ex, x, a, order);
    }
    friend inline Matrix<Scalar> taylor(const Matrix<Scalar>& ex, const Matrix<Scalar>& x) {
      return Matrix<Scalar>::taylor(ex, x, 0, 1);
    }
    ///@}

    /**
     * \brief multivariate Taylor series expansion
     *
     * Do Taylor expansions until the aggregated order of a term is equal to 'order'.
     * The aggregated order of \f$x^n y^m\f$ equals \f$n+m\f$.
     *
     */
    friend inline Matrix<Scalar> mtaylor(const Matrix<Scalar>& ex, const Matrix<Scalar>& x,
                                           const Matrix<Scalar>& a, int order=1) {
      return Matrix<Scalar>::mtaylor(ex, x, a, order);
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
    friend inline Matrix<Scalar> mtaylor(const Matrix<Scalar>& ex, const Matrix<Scalar>& x,
                                           const Matrix<Scalar>& a, int order,
                                           const std::vector<int>& order_contributions) {
      return Matrix<Scalar>::mtaylor(ex, x, a, order, order_contributions);
    }

    /** \brief extracts polynomial coefficients from an expression
     *
     * \param ex Scalar expression that represents a polynomial
     * \param x  Scalar symbol that the polynomial is build up with
     */
    friend inline Matrix<Scalar> poly_coeff(const Matrix<Scalar>& f,
                                              const Matrix<Scalar>& x) {
      return Matrix<Scalar>::poly_coeff(f, x);
    }

    /** \brief Attempts to find the roots of a polynomial
     *
     *  This will only work for polynomials up to order 3
     *  It is assumed that the roots are real.
     *
     */
    friend inline Matrix<Scalar> poly_roots(const Matrix<Scalar>& p) {
      return Matrix<Scalar>::poly_roots(p);
    }

    /** \brief Attempts to find the eigenvalues of a symbolic matrix
     *  This will only work for up to 3x3 matrices
     */
    friend inline Matrix<Scalar> eig_symbolic(const Matrix<Scalar>& m) {
      return Matrix<Scalar>::eig_symbolic(m);
    }
/** @} */
#endif

    /** \brief Set or reset the depth to which equalities are being checked for simplifications */
    static void setEqualityCheckingDepth(int eq_depth=1);

    /** \brief Get the depth to which equalities are being checked for simplifications */
    static int getEqualityCheckingDepth();

    /** \brief Get function input */
    static std::vector<Matrix<Scalar> > get_input(const Function& f);

    /** \brief Get free */
    static std::vector<Matrix<Scalar> > get_free(const Function& f);

    ///@{
    /** \brief Jacobian expression */
    static Matrix<Scalar> jac(const Function& f, int iind=0, int oind=0,
                  bool compact=false, bool symmetric=false);
    static Matrix<Scalar> jac(const Function& f, const std::string & iname, int oind=0,
                  bool compact=false, bool symmetric=false);
    static Matrix<Scalar> jac(const Function& f, int iind, const std::string& oname,
                  bool compact=false, bool symmetric=false);
    static Matrix<Scalar> jac(const Function& f, const std::string& iname,
                                const std::string& oname,
                                bool compact=false, bool symmetric=false);
    ///@}

    ///@{
    /** \brief Gradient expression */
    static Matrix<Scalar> grad(const Function& f, int iind=0, int oind=0);
    static Matrix<Scalar> grad(const Function& f, const std::string& iname, int oind=0);
    static Matrix<Scalar> grad(const Function& f, int iind, const std::string& oname);
    static Matrix<Scalar> grad(const Function& f, const std::string& iname,
                                 const std::string& oname);
    ///@}

    ///@{
    /** \brief Tangent expression */
    static Matrix<Scalar> tang(const Function& f, int iind=0, int oind=0);
    static Matrix<Scalar> tang(const Function& f, const std::string& iname, int oind=0);
    static Matrix<Scalar> tang(const Function& f, int iind, const std::string& oname);
    static Matrix<Scalar> tang(const Function& f, const std::string& iname,
                                 const std::string& oname);
    ///@}

    ///@{
    /** \bried Hessian expression */
    static Matrix<Scalar> hess(const Function& f, int iind=0, int oind=0);
    static Matrix<Scalar> hess(const Function& f, const std::string& iname, int oind=0);
    static Matrix<Scalar> hess(const Function& f, int iind, const std::string& oname);
    static Matrix<Scalar> hess(const Function& f, const std::string& iname,
                                 const std::string& oname);
    ///@}

    /// Get name of the class
    static std::string type_name();

    /// Print a description of the object
    void print(std::ostream &stream=casadi::userOut(), bool trailing_newline=true) const;

    /// Get strings corresponding to the nonzeros and the interdependencies
    void print_split(std::vector<std::string>& SWIG_OUTPUT(nz),
                    std::vector<std::string>& SWIG_OUTPUT(inter)) const;

    /// Print a representation of the object
    void repr(std::ostream &stream=casadi::userOut(), bool trailing_newline=true) const;

    /// Print scalar
    void print_scalar(std::ostream &stream=casadi::userOut(), bool trailing_newline=true) const;

    /// Print vector-style
    void print_vector(std::ostream &stream=casadi::userOut(), bool trailing_newline=true) const;

    /// Print dense matrix-stype
    void print_dense(std::ostream &stream=casadi::userOut(), bool trailing_newline=true) const;

    /// Print sparse matrix style
    void print_sparse(std::ostream &stream=casadi::userOut(), bool trailing_newline=true) const;

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
    ///@{
    /// Access the non-zero elements
    std::vector<Scalar>& nonzeros() { return nonzeros_;}
    const std::vector<Scalar>& nonzeros() const { return nonzeros_;}
    ///@}

    ///@{
    /// Get a pointer to the data
    Scalar* ptr() { return nonzeros_.empty() ? 0 : &nonzeros_.front(); }
    const Scalar* ptr() const { return nonzeros_.empty() ? 0 : &nonzeros_.front(); }
    friend inline Scalar* get_ptr(Matrix<Scalar>& v) { return v.ptr(); }
    friend inline const Scalar* get_ptr(const Matrix<Scalar>& v) { return v.ptr(); }
    ///@}

    /// Const access the sparsity - reference to data member
    const Sparsity& sparsity() const { return sparsity_; }

#endif // SWIG

    /** \brief Get an owning reference to the sparsity pattern */
    Sparsity get_sparsity() const { return sparsity();}

    /* \brief Construct a sparse matrix from triplet form
     * Default matrix size is max(col) x max(row)
     */
    ///@{
    static Matrix<Scalar> triplet(const std::vector<int>& row, const std::vector<int>& col,
                                    const Matrix<Scalar>& d);
    static Matrix<Scalar> triplet(const std::vector<int>& row, const std::vector<int>& col,
                                    const Matrix<Scalar>& d, int nrow, int ncol);
    static Matrix<Scalar> triplet(const std::vector<int>& row, const std::vector<int>& col,
                                    const Matrix<Scalar>& d, const std::pair<int, int>& rc);
    ///@}

    ///@{
    /** \brief  create a matrix with all inf */
    static Matrix<Scalar> inf(const Sparsity& sp);
    static Matrix<Scalar> inf(int nrow=1, int ncol=1);
    static Matrix<Scalar> inf(const std::pair<int, int>& rc);
    ///@}

    ///@{
    /** \brief  create a matrix with all nan */
    static Matrix<Scalar> nan(const Sparsity& sp);
    static Matrix<Scalar> nan(int nrow=1, int ncol=1);
    static Matrix<Scalar> nan(const std::pair<int, int>& rc);
    ///@}

    /** \brief  create an n-by-n identity matrix */
    static Matrix<Scalar> eye(int ncol);

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

    /** \brief Get all nonzeros */
    std::vector<Scalar> get_nonzeros() const { return nonzeros_;}

#ifndef SWIG
    /** \brief Get all nonzeros */
    template<typename A>
    std::vector<A> get_nonzeros() const;
#endif // SWIG

    /** \brief Type conversion to double */
    explicit operator double() const;

    /** \brief Type conversion to int */
    explicit operator int() const;

#ifndef SWIG
    /** \brief Type conversion to a vector */
    template<typename A>
    explicit operator std::vector<A>() const;
#endif // SWIG

    /** \brief Get name (only if symbolic scalar) */
    std::string name() const;

    /** \brief Get expressions of the children of the expression
        Only defined if symbolic scalar.
        Wraps SXElem SXElem::dep(int ch=0) const.
     */
    Matrix<Scalar> dep(int ch=0) const;

    /** \brief Get the number of dependencies of a binary SXElem
        Only defined if symbolic scalar.
    */
    int n_dep() const;

    // @{
    /// Set the 'precision, width & scientific' used in printing and serializing to streams
    static void setPrecision(int precision) { stream_precision_ = precision; }
    static void setWidth(int width) { stream_width_ = width; }
    static void setScientific(bool scientific) { stream_scientific_ = scientific; }
    // @}

#ifndef SWIG
    /// Sparse matrix with a given sparsity with all values same
    Matrix(const Sparsity& sp, const Scalar& val, bool dummy);

    /// Sparse matrix with a given sparsity and non-zero elements.
    Matrix(const Sparsity& sp, const std::vector<Scalar>& d, bool dummy);

  private:
    /// Sparsity of the matrix in a compressed column storage (CCS) format
    Sparsity sparsity_;

    /// Nonzero elements
    std::vector<Scalar> nonzeros_;

    /// Precision used in streams
    static int stream_precision_;
    static int stream_width_;
    static bool stream_scientific_;
#endif // SWIG
  };

  ///@{
  /// Readability typedefs
  typedef Matrix<int> IM;
  typedef Matrix<double> DM;
  typedef std::vector<DM> DMVector;
  typedef std::vector<DMVector> DMVectorVector;
  typedef std::map<std::string, DM> DMDict;
  ///@}

  /// Is the IM a Slice
  bool CASADI_EXPORT is_slice(const IM& x, bool ind1=false);

  ///  Convert IM to Slice
  Slice CASADI_EXPORT to_slice(const IM& x, bool ind1=false);


} // namespace casadi

#endif // CASADI_MATRIX_HPP
