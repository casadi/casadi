/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010-2023 Joel Andersson, Joris Gillis, Moritz Diehl,
 *                            KU Leuven. All rights reserved.
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


#ifndef CASADI_MATRIX_DECL_HPP
#define CASADI_MATRIX_DECL_HPP

#include "matrix_fwd.hpp"
#include "exception.hpp"
#include "casadi_limits.hpp"
#include "runtime/casadi_runtime.hpp"
#include "generic_matrix.hpp"
#include "generic_expression.hpp"
#include "printable.hpp"

#include <random>
#include <typeinfo>
#include <vector>
#include <initializer_list>
#ifdef CASADI_WITH_THREAD
#ifdef CASADI_WITH_THREAD_MINGW
#include <mingw.mutex.h>
#else // CASADI_WITH_THREAD_MINGW
#include <mutex>
#endif // CASADI_WITH_THREAD_MINGW
#endif //CASADI_WITH_THREAD

namespace casadi {

  /** \brief Empty Base

      This class is extended in SWIG.

      \identifier{18c} */
  struct CASADI_EXPORT MatrixCommon {};

/// \cond CLUTTER
  ///@{
  /** \brief Get typename

      \identifier{18d} */
  template <typename Scalar> inline std::string matrixName()
  { return std::string("Matrix<") + typeid(Scalar).name() + std::string(">");}
  template<> inline std::string matrixName<double>() { return "DM"; }
  template<> inline std::string matrixName<casadi_int>() { return "IM"; }
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

      \identifier{18e} */
  template<typename Scalar>
  class CASADI_EXPORT Matrix :
    public MatrixCommon,
    public SWIG_IF_ELSE(GenericExpressionCommon, GenericExpression<Matrix<Scalar> >),
    public GenericMatrix<Matrix<Scalar> >,
    public SWIG_IF_ELSE(PrintableCommon, Printable<Matrix<Scalar> >) {
  public:

    /** \brief  constructors

        \identifier{18f} */
    /// empty 0-by-0 matrix constructor
    Matrix();

    /// Copy constructor
    Matrix(const Matrix<Scalar>& m);

#ifndef SWIG
    /// Assignment (normal)
    Matrix<Scalar>& operator=(const Matrix<Scalar>& m);
#endif // SWIG

    /** \brief Create a sparse matrix with all structural zeros

        \identifier{18g} */
    Matrix(casadi_int nrow, casadi_int ncol);

#ifndef SWIG
    /** \brief Create a sparse matrix with all structural zeros

        \identifier{18h} */
    explicit Matrix(const std::pair<casadi_int, casadi_int>& rc);

    /** \brief  Access functions of the node

        \identifier{18i} */
    std::vector<Scalar>* operator->();

    /** \brief  Const access functions of the node

        \identifier{18j} */
    const std::vector<Scalar>* operator->() const;
#endif // SWIG

    /** \brief Create a sparse matrix from a sparsity pattern.

        Same as Matrix::ones(sparsity)

        \identifier{18k} */
    explicit Matrix(const Sparsity& sp);

    /** \brief Construct matrix with a given sparsity and nonzeros

        \identifier{18l} */
    Matrix(const Sparsity& sp, const Matrix<Scalar>& d);

    /// This constructor enables implicit type conversion from a numeric type
    Matrix(double val);

#if !(defined(SWIG) && defined(SWIGMATLAB))
    /// Dense matrix constructor with data given as vector of vectors
    explicit Matrix(const std::vector< std::vector<double> >& m);

    /** \brief  Create an expression from a vector

        \identifier{18m} */
    template<typename A>
    Matrix(const std::vector<A>& x) : sparsity_(Sparsity::dense(x.size(), 1)),
      nonzeros_(std::vector<Scalar>(x.size())) {
        auto x_it = x.begin();
        for (auto&& d : nonzeros_) d = static_cast<Scalar>(*x_it++);
    }
#endif

    /** \brief Create a matrix from another matrix with a different entry type

     *  Assumes that the scalar conversion is valid.

        \identifier{18n} */
    template<typename A>
    Matrix(const Matrix<A>& x) : sparsity_(x.sparsity()), nonzeros_(std::vector<Scalar>(x.nnz())) {
      auto x_it = x->begin();
      for (auto&& d : nonzeros_) d = static_cast<Scalar>(*x_it++);
    }

#ifndef SWIG
    /// Construct from a vector
    Matrix(const std::vector<Scalar>& x);

    /// Construct from initializer list
    Matrix(std::initializer_list<Scalar> x) : Matrix<Scalar>(std::vector<Scalar>(x)) {}

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
    bool has_nz(casadi_int rr, casadi_int cc) const;

    /// Returns the truth value of a Matrix
    bool __nonzero__() const;

    ///@{
    /// Get a submatrix, single argument
    void get(Matrix<Scalar>& SWIG_OUTPUT(m), bool ind1, const Slice& rr) const;
    void get(Matrix<Scalar>& SWIG_OUTPUT(m), bool ind1, const Matrix<casadi_int>& rr) const;
    void get(Matrix<Scalar>& SWIG_OUTPUT(m), bool ind1, const Sparsity& sp) const;
    ///@}

    /// Get a submatrix, two arguments
    ///@{
    void get(Matrix<Scalar>& SWIG_OUTPUT(m), bool ind1,
                const Slice& rr, const Slice& cc) const;
    void get(Matrix<Scalar>& SWIG_OUTPUT(m), bool ind1,
                const Slice& rr, const Matrix<casadi_int>& cc) const;
    void get(Matrix<Scalar>& SWIG_OUTPUT(m), bool ind1,
                const Matrix<casadi_int>& rr, const Slice& cc) const;
    void get(Matrix<Scalar>& SWIG_OUTPUT(m), bool ind1,
                const Matrix<casadi_int>& rr, const Matrix<casadi_int>& cc) const;
    ///@}

    ///@{
    /// Set a submatrix, single argument
    void set(const Matrix<Scalar>& m, bool ind1, const Slice& rr);
    void set(const Matrix<Scalar>& m, bool ind1, const Matrix<casadi_int>& rr);
    void set(const Matrix<Scalar>& m, bool ind1, const Sparsity& sp);
    ///@}

    ///@{
    /// Set a submatrix, two arguments
    void set(const Matrix<Scalar>& m, bool ind1, const Slice& rr, const Slice& cc);
    void set(const Matrix<Scalar>& m, bool ind1, const Slice& rr, const Matrix<casadi_int>& cc);
    void set(const Matrix<Scalar>& m, bool ind1, const Matrix<casadi_int>& rr, const Slice& cc);
    void set(const Matrix<Scalar>& m, bool ind1, const Matrix<casadi_int>& rr,
                                                 const Matrix<casadi_int>& cc);
    ///@}

    ///@{
    /// Get a set of nonzeros
    void get_nz(Matrix<Scalar>& SWIG_OUTPUT(m), bool ind1, const Slice& k) const;
    void get_nz(Matrix<Scalar>& SWIG_OUTPUT(m), bool ind1, const Matrix<casadi_int>& k) const;
    ///@}

    ///@{
    /// Set a set of nonzeros
    void set_nz(const Matrix<Scalar>& m, bool ind1, const Slice& k);
    void set_nz(const Matrix<Scalar>& m, bool ind1, const Matrix<casadi_int>& k);
    ///@}

    Matrix<Scalar> operator+() const;
    Matrix<Scalar> operator-() const;

    /// \cond INTERNAL
    ///@{
    /** \brief  Create nodes by their ID

        \identifier{18o} */
    static Matrix<Scalar> binary(casadi_int op, const Matrix<Scalar> &x, const Matrix<Scalar> &y);
    static Matrix<Scalar> unary(casadi_int op, const Matrix<Scalar> &x);
    static Matrix<Scalar> scalar_matrix(casadi_int op,
                                          const Matrix<Scalar> &x, const Matrix<Scalar> &y);
    static Matrix<Scalar> matrix_scalar(casadi_int op,
                                          const Matrix<Scalar> &x, const Matrix<Scalar> &y);
    static Matrix<Scalar> matrix_matrix(casadi_int op,
                                          const Matrix<Scalar> &x, const Matrix<Scalar> &y);
    static std::vector< Matrix<Scalar> > call(const Function& f,
        const std::vector< Matrix<Scalar> > &x);
    static std::vector< Scalar > call(const Function& f, const std::vector< Scalar > &x);
    ///@}
    /// \endcond

#ifndef SWIG
    /// \cond CLUTTER
    ///@{
    /// Functions called by friend functions defined for GenericExpression
    static bool is_equal(const Matrix<Scalar> &x, const Matrix<Scalar> &y, casadi_int depth=0);
    static Matrix<Scalar> mmin(const Matrix<Scalar> &x);
    static Matrix<Scalar> mmax(const Matrix<Scalar> &x);
    ///@}

    ///@{
    /// Functions called by friend functions defined for GenericMatrix
    static Matrix<Scalar> simplify(const Matrix<Scalar> &x);
    static Matrix<Scalar> jacobian(const Matrix<Scalar> &f, const Matrix<Scalar> &x,
                                   const Dict& opts = Dict());
    static Sparsity jacobian_sparsity(const Matrix<Scalar> &f, const Matrix<Scalar> &x);
    static Matrix<Scalar> hessian(const Matrix<Scalar> &f, const Matrix<Scalar> &x,
      const Dict& opts = Dict());
    static Matrix<Scalar> hessian(const Matrix<Scalar> &f, const Matrix<Scalar> &x,
                                    Matrix<Scalar>& g, const Dict& opts = Dict());
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
    static Matrix<Scalar> expm_const(const Matrix<Scalar> &A, const Matrix<Scalar> &t);
    static Matrix<Scalar> expm(const Matrix<Scalar> &A);
    static Matrix<Scalar> solve(const Matrix<Scalar> &A, const Matrix<Scalar>& b);
    static Matrix<Scalar> solve(const Matrix<Scalar> &A, const Matrix<Scalar>& b,
                                  const std::string& lsolver, const Dict& opts);
    static Matrix<Scalar> inv(const Matrix<Scalar> &A);
    static Matrix<Scalar> inv(const Matrix<Scalar> &A,
                                  const std::string& lsolver, const Dict& opts);

    static casadi_int n_nodes(const Matrix<Scalar> &x);
    static std::string print_operator(const Matrix<Scalar> &x,
                                      const std::vector<std::string>& args);
    static void extract(std::vector<Matrix<Scalar>>& ex, std::vector<Matrix<Scalar>>& v,
      std::vector<Matrix<Scalar>>& vdef, const Dict& opts = Dict());
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
                                    bool short_circuit=false);
    static Matrix<Scalar> conditional(const Matrix<Scalar> &ind,
                                        const std::vector<Matrix<Scalar> > &x,
                                        const Matrix<Scalar> &x_default,
                                        bool short_circuit=false);
    static bool depends_on(const Matrix<Scalar> &x, const Matrix<Scalar> &arg);
    static bool contains(const std::vector<Matrix<Scalar> >& v, const Matrix<Scalar> &n);
    static bool contains_all(const std::vector<Matrix<Scalar> >& v,
        const std::vector<Matrix<Scalar> > &n);
    static bool contains_any(const std::vector<Matrix<Scalar> >& v,
        const std::vector<Matrix<Scalar> > &n);

    static Matrix<Scalar> mrdivide(const Matrix<Scalar> &x, const Matrix<Scalar> &y);
    static Matrix<Scalar> mldivide(const Matrix<Scalar> &x, const Matrix<Scalar> &y);
    static std::vector<Matrix<Scalar> > symvar(const Matrix<Scalar> &x);
    static Matrix<Scalar> det(const Matrix<Scalar> &x);
    static Matrix<Scalar> inv_minor(const Matrix<Scalar> &x);
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
    static Matrix<Scalar> einstein(const Matrix<Scalar>& A, const Matrix<Scalar>& B,
      const Matrix<Scalar>& C,
      const std::vector<casadi_int>& dim_a, const std::vector<casadi_int>& dim_b,
      const std::vector<casadi_int>& dim_c,
      const std::vector<casadi_int>& a, const std::vector<casadi_int>& b,
      const std::vector<casadi_int>& c);

    static Matrix<Scalar> einstein(const Matrix<Scalar>& A, const Matrix<Scalar>& B,
      const std::vector<casadi_int>& dim_a, const std::vector<casadi_int>& dim_b,
      const std::vector<casadi_int>& dim_c,
      const std::vector<casadi_int>& a, const std::vector<casadi_int>& b,
      const std::vector<casadi_int>& c);
    static Matrix<Scalar> cumsum(const Matrix<Scalar> &x, casadi_int axis=-1);
    static Matrix<Scalar> _logsumexp(const Matrix<Scalar>& x);
    static std::vector< Matrix<Scalar> > cse(const std::vector< Matrix<Scalar> >& e);
    ///@}

    ///@{
    /// Functions called by friend functions defined for SparsityInterface
    static Matrix<Scalar> blockcat(const std::vector< std::vector<Matrix<Scalar> > > &v);
    static Matrix<Scalar> horzcat(const std::vector<Matrix<Scalar> > &v);
    static std::vector<Matrix<Scalar> >
      horzsplit(const Matrix<Scalar>& x,
                const std::vector<casadi_int>& offset);
    static Matrix<Scalar> vertcat(const std::vector<Matrix<Scalar> > &v);
    static std::vector< Matrix<Scalar> >
      vertsplit(const Matrix<Scalar>& x,
                const std::vector<casadi_int>& offset);
    static std::vector< Matrix<Scalar> >
      diagsplit(const Matrix<Scalar>& x,
                const std::vector<casadi_int>& offset1,
                const std::vector<casadi_int>& offset2);
    static Matrix<Scalar> reshape(const Matrix<Scalar> &x, casadi_int nrow, casadi_int ncol);
    static Matrix<Scalar> reshape(const Matrix<Scalar> &x, const Sparsity& sp);
    static Matrix<Scalar> sparsity_cast(const Matrix<Scalar> &x, const Sparsity& sp);
    static Matrix<Scalar> kron(const Matrix<Scalar> &x, const Matrix<Scalar>& y);
    static Matrix<Scalar> mtimes(const Matrix<Scalar> &x, const Matrix<Scalar> &y);
    static Matrix<Scalar> mac(const Matrix<Scalar> &x,
                                const Matrix<Scalar> &y,
                                const Matrix<Scalar> &z);
    static void extract_parametric(const Matrix<Scalar> &expr,
        const Matrix<Scalar>& par,
        Matrix<Scalar>& expr_ret,
        std::vector<Matrix<Scalar> >& symbols,
        std::vector<Matrix<Scalar>>& parametric,
        const Dict& opts);
    static void separate_linear(const Matrix<Scalar> &expr,
      const Matrix<Scalar> &sym_lin, const Matrix<Scalar> &sym_const,
      Matrix<Scalar>& expr_const, Matrix<Scalar>& expr_lin, Matrix<Scalar>& expr_nonlin);
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
                                             const Matrix<Scalar> &b, casadi_int order=5);
    static Matrix<Scalar> gauss_quadrature(const Matrix<Scalar> &f,
                                             const Matrix<Scalar> &x, const Matrix<Scalar> &a,
                                             const Matrix<Scalar> &b, casadi_int order,
                                             const Matrix<Scalar>& w);
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
        casadi_int order=1, bool tr=false);
    static Matrix<Scalar> taylor(const Matrix<Scalar>& ex, const Matrix<Scalar>& x,
                                   const Matrix<Scalar>& a, casadi_int order);
    static Matrix<Scalar> mtaylor(const Matrix<Scalar>& ex, const Matrix<Scalar>& x,
                                    const Matrix<Scalar>& a, casadi_int order);
    static Matrix<Scalar> mtaylor(const Matrix<Scalar>& ex,
                                    const Matrix<Scalar>& x, const Matrix<Scalar>& a,
                                    casadi_int order,
                                    const std::vector<casadi_int>& order_contributions);
    static Matrix<Scalar> poly_coeff(const Matrix<Scalar>& ex, const Matrix<Scalar>&x);
    static Matrix<Scalar> poly_roots(const Matrix<Scalar>& p);
    static Matrix<Scalar> eig_symbolic(const Matrix<Scalar>& m);
    static Matrix<double> evalf(const Matrix<Scalar>& m);
    static void qr_sparse(const Matrix<Scalar>& A, Matrix<Scalar>& V, Matrix<Scalar>& R,
                          Matrix<Scalar>& beta, std::vector<casadi_int>& prinv,
                          std::vector<casadi_int>& pc, bool amd=true);
    static Matrix<Scalar> qr_solve(const Matrix<Scalar>& b, const Matrix<Scalar>& v,
                                   const Matrix<Scalar>& r, const Matrix<Scalar>& beta,
                                   const std::vector<casadi_int>& prinv,
                                   const std::vector<casadi_int>& pc, bool tr=false);
    static void qr(const Matrix<Scalar>& A, Matrix<Scalar>& Q, Matrix<Scalar>& R);
    static void ldl(const Matrix<Scalar>& A, Matrix<Scalar>& D, Matrix<Scalar>& LT,
                    std::vector<casadi_int>& p, bool amd=true);
    static Matrix<Scalar> ldl_solve(const Matrix<Scalar>& b, const Matrix<Scalar>& D,
                                    const Matrix<Scalar>& LT, const std::vector<casadi_int>& p);
    static Matrix<Scalar> all(const Matrix<Scalar>& x);
    static Matrix<Scalar> any(const Matrix<Scalar>& x);
    static Matrix<Scalar> adj(const Matrix<Scalar>& x);
    static Matrix<Scalar> minor(const Matrix<Scalar>& x, casadi_int i, casadi_int j);
    static Matrix<Scalar> cofactor(const Matrix<Scalar>& A, casadi_int i, casadi_int j);
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
\addtogroup expression_tools
@{
*/
    /** \brief Matrix adjoint

        \identifier{18p} */
    friend inline Matrix<Scalar> adj(const Matrix<Scalar>& A) {
      return Matrix<Scalar>::adj(A);
    }

    /** \brief Get the (i,j) minor matrix

        \identifier{18q} */
    friend inline Matrix<Scalar> minor(const Matrix<Scalar> &x, casadi_int i, casadi_int j) {
      return Matrix<Scalar>::minor(x, i, j);
    }

    /** \brief Get the (i,j) cofactor matrix

        \identifier{18r} */
    friend inline Matrix<Scalar> cofactor(const Matrix<Scalar> &x, casadi_int i, casadi_int j) {
      return Matrix<Scalar>::cofactor(x, i, j);
    }

    /** \brief  QR factorization using the modified Gram-Schmidt algorithm

     * More stable than the classical Gram-Schmidt, but may break down if the rows of A
     * are nearly linearly dependent
     * See J. Demmel: Applied Numerical Linear Algebra (algorithm 3.1.).
     * Note that in SWIG, Q and R are returned by value.

        \identifier{18s} */
    friend inline void qr(const Matrix<Scalar>& A, Matrix<Scalar>& Q, Matrix<Scalar>& R) {
      return Matrix<Scalar>::qr(A, Q, R);
    }

    /** \brief Sparse direct QR factorization

     * See T. Davis: Direct Methods for Sparse Linear Systems

        \identifier{18t} */
    friend inline void qr_sparse(const Matrix<Scalar>& A, Matrix<Scalar>& V, Matrix<Scalar>& R,
                                 Matrix<Scalar>& beta, std::vector<casadi_int>& prinv,
                                 std::vector<casadi_int>& pc, bool amd=true) {
      return Matrix<Scalar>::qr_sparse(A, V, R, beta, prinv, pc, amd);
    }

    /** \brief Solve using a sparse QR factorization

        \identifier{18u} */
    friend inline Matrix<Scalar>
    qr_solve(const Matrix<Scalar>& b, const Matrix<Scalar>& v,
             const Matrix<Scalar>& r, const Matrix<Scalar>& beta,
             const std::vector<casadi_int>& prinv,
             const std::vector<casadi_int>& pc, bool tr=false) {
        return Matrix<Scalar>::qr_solve(b, v, r, beta, prinv, pc, tr);
    }

    /** \brief Obtain a Cholesky factorisation of a matrix

     * Performs and LDL transformation [L,D] = ldl(A) and
     * returns diag(sqrt(D))*L'

        \identifier{18v} */
    friend inline Matrix<Scalar> chol(const Matrix<Scalar>& A) {
      return Matrix<Scalar>::chol(A);
    }

    /** \brief Sparse LDL^T factorization

     * Returns D and the strictly upper triangular entries of L^T
     * I.e. ones on the diagonal are ignored.
     * Only guarenteed to work for positive definite matrices.

        \identifier{18w} */
    friend inline void ldl(const Matrix<Scalar>& A, Matrix<Scalar>& D, Matrix<Scalar>& LT,
                           std::vector<casadi_int>& p, bool amd=true) {
      return Matrix<Scalar>::ldl(A, D, LT, p, amd);
    }

    /** \brief Solve using a sparse LDL^T factorization

        \identifier{18x} */
    friend inline Matrix<Scalar>
    ldl_solve(const Matrix<Scalar>& b, const Matrix<Scalar>& D, const Matrix<Scalar>& LT,
              const std::vector<casadi_int>& p) {
      return Matrix<Scalar>::ldl_solve(b, D, LT, p);
    }

    /** \brief Returns true only if any element in the matrix is true

        \identifier{18y} */
    friend inline Matrix<Scalar> any(const Matrix<Scalar> &x) {
      return Matrix<Scalar>::any(x);
    }

    /** \brief Returns true only if every element in the matrix is true

        \identifier{18z} */
    friend inline Matrix<Scalar> all(const Matrix<Scalar> &x) {
      return Matrix<Scalar>::all(x);
    }

    /** \brief Inf-norm of a Matrix-Matrix product

        \identifier{190} */
    friend inline Matrix<Scalar>
      norm_inf_mul(const Matrix<Scalar> &x, const Matrix<Scalar> &y) {
      return Matrix<Scalar>::norm_inf_mul(x, y);
    }

    /** \brief  Make a matrix sparse by removing numerical zeros

        \identifier{191} */
    friend inline Matrix<Scalar>
      sparsify(const Matrix<Scalar>& A, double tol=0) {
      return Matrix<Scalar>::sparsify(A, tol);
    }

    /** \brief  Expand the expression as a weighted sum (with constant weights)

        \identifier{192} */
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

        \identifier{193} */
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

                                                                                                                \identifier{194} */
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

         \identifier{195} */
    friend inline Matrix<Scalar> heaviside(const Matrix<Scalar> &x) {
      return Matrix<Scalar>::heaviside(x);
    }

    /** \brief rectangle function
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

        \identifier{23n} */
    friend inline Matrix<Scalar> rectangle(const Matrix<Scalar> &x) {
      return Matrix<Scalar>::rectangle(x);
    }

    /** \brief triangle function
     *
     * \f[
     * \begin {cases}
     * \Lambda(x) = 0 &    |x| >= 1  \\
     * \Lambda(x) = 1-|x| &  |x| < 1
     * \end {cases}
     * \f]
     *

        \identifier{23o} */
    friend inline Matrix<Scalar> triangle(const Matrix<Scalar> &x) {
      return Matrix<Scalar>::triangle(x);
    }

    /** \brief ramp function
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

        \identifier{23p} */
    friend inline Matrix<Scalar> ramp(const Matrix<Scalar> &x) {
      return Matrix<Scalar>::ramp(x);
    }

    ///@{
    /** \brief  Integrate f from a to b using Gaussian quadrature with n points

        \identifier{196} */
    friend inline Matrix<Scalar>
      gauss_quadrature(const Matrix<Scalar> &f, const Matrix<Scalar> &x,
                       const Matrix<Scalar> &a, const Matrix<Scalar> &b,
                       casadi_int order=5) {
      return Matrix<Scalar>::gauss_quadrature(f, x, a, b, order);
    }
    friend inline Matrix<Scalar>
      gauss_quadrature(const Matrix<Scalar> &f, const Matrix<Scalar> &x,
                       const Matrix<Scalar> &a, const Matrix<Scalar> &b,
                       casadi_int order, const Matrix<Scalar>& w) {
      return Matrix<Scalar>::gauss_quadrature(f, x, a, b, order, w);
    }
    ///@}

    ///@{
    /** \brief univariate Taylor series expansion
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
     *
     * \code
     * taylor(sin(x),x,x0) -> sin(x0)+cos(x0)*(x-x0)
     * \endcode
     *
     * \sa linearize, mtaylor, linear_coeff

        \identifier{23q} */
    friend inline Matrix<Scalar> taylor(const Matrix<Scalar>& ex, const Matrix<Scalar>& x,
                                          const Matrix<Scalar>& a, casadi_int order=1) {
      return Matrix<Scalar>::taylor(ex, x, a, order);
    }
    friend inline Matrix<Scalar> taylor(const Matrix<Scalar>& ex, const Matrix<Scalar>& x) {
      return Matrix<Scalar>::taylor(ex, x, 0, 1);
    }
    ///@}

    /** \brief multivariate Taylor series expansion
     *
     * Do Taylor expansions until the aggregated order of a term is equal to 'order'.
     * The aggregated order of \f$x^n y^m\f$ equals \f$n+m\f$.
     *
     * \sa taylor

        \identifier{23r} */
    friend inline Matrix<Scalar> mtaylor(const Matrix<Scalar>& ex, const Matrix<Scalar>& x,
                                           const Matrix<Scalar>& a, casadi_int order=1) {
      return Matrix<Scalar>::mtaylor(ex, x, a, order);
    }

    /** \brief multivariate Taylor series expansion
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
     * \sa taylor

        \identifier{23s} */
    friend inline Matrix<Scalar> mtaylor(const Matrix<Scalar>& ex, const Matrix<Scalar>& x,
                                           const Matrix<Scalar>& a, casadi_int order,
                                           const std::vector<casadi_int>& order_contributions) {
      return Matrix<Scalar>::mtaylor(ex, x, a, order, order_contributions);
    }

    /** \brief extracts polynomial coefficients from an expression
     *
     * \param ex Scalar expression that represents a polynomial
     * \param x  Scalar symbol that the polynomial is build up with

        \identifier{197} */
    friend inline Matrix<Scalar> poly_coeff(const Matrix<Scalar>& f,
                                              const Matrix<Scalar>& x) {
      return Matrix<Scalar>::poly_coeff(f, x);
    }

    /** \brief Attempts to find the roots of a polynomial
     *
     *  This will only work for polynomials up to order 3
     *  It is assumed that the roots are real.
     *
        \identifier{198} */
    friend inline Matrix<Scalar> poly_roots(const Matrix<Scalar>& p) {
      return Matrix<Scalar>::poly_roots(p);
    }

    /** \brief Attempts to find the eigenvalues of a symbolic matrix

     *  This will only work for up to 3x3 matrices

        \identifier{199} */
    friend inline Matrix<Scalar> eig_symbolic(const Matrix<Scalar>& m) {
      return Matrix<Scalar>::eig_symbolic(m);
    }


    /** \brief Evaluates the expression numerically
    *
    * An error is raised when the expression contains symbols

        \identifier{19a} */
    inline friend Matrix<double> evalf(const Matrix<Scalar>& expr) {
      return Matrix<Scalar>::evalf(expr);
    }
/** @} */
#endif

    /** \brief Set or reset the depth to which equalities are being checked for simplifications

        \identifier{19b} */
    static void set_max_depth(casadi_int eq_depth=1);

    /** \brief Get the depth to which equalities are being checked for simplifications

        \identifier{19c} */
    static casadi_int get_max_depth();

    /** \brief Get function input

        \identifier{19d} */
    static std::vector<Matrix<Scalar> > get_input(const Function& f);

    /** \brief Get free

        \identifier{19e} */
    static std::vector<Matrix<Scalar> > get_free(const Function& f);

    /// Get name of the class
    static std::string type_name();

    /// Get strings corresponding to the nonzeros and the interdependencies
    void print_split(std::vector<std::string>& SWIG_OUTPUT(nz),
                    std::vector<std::string>& SWIG_OUTPUT(inter)) const;

    /// Print a representation of the object
    void disp(std::ostream& stream, bool more=false) const;

    /// Get string representation
    std::string get_str(bool more=false) const;

    /// Print scalar
    void print_scalar(std::ostream &stream) const;

    /// Print vector-style
    void print_vector(std::ostream &stream, bool truncate=true) const;

    /// Print dense matrix-stype
    void print_dense(std::ostream &stream, bool truncate=true) const;

    /// Print sparse matrix style
    void print_sparse(std::ostream &stream, bool truncate=true) const;

#ifndef SWIG
    /// Print default style
    static void print_default(std::ostream &stream, const Sparsity& sp, const Scalar* nonzeros,
      bool truncate=true);

    /// Print scalar
    static void print_scalar(std::ostream &stream, const Scalar& e);

    /// Print vector-style
    static void print_vector(std::ostream &stream, const Sparsity& sp, const Scalar* nonzeros,
      bool truncate=true);

    /// Print canonical style
    static void print_canonical(std::ostream &stream, const Sparsity& sp, const Scalar* nonzeros,
        bool truncate=true);

    /// Print scalar
    static void print_sparse(std::ostream &stream, const Sparsity& sp, const Scalar* nonzeros,
      bool truncate=true);

    /// Get strings corresponding to the nonzeros and the interdependencies
    static void print_split(casadi_int nnz, const Scalar* nonzeros, std::vector<std::string>& nz,
                    std::vector<std::string>& inter);

    /// Print dense matrix-stype
    static void print_dense(std::ostream &stream,  const Sparsity& sp, const Scalar* nonzeros,
      bool truncate=true);
#endif

    void clear();
    void resize(casadi_int nrow, casadi_int ncol);
    void reserve(casadi_int nnz);
    void reserve(casadi_int nnz, casadi_int ncol);

    /** \brief Erase a submatrix (leaving structural zeros in its place)

        Erase rows and/or columns of a matrix

        \identifier{19f} */
    void erase(const std::vector<casadi_int>& rr, const std::vector<casadi_int>& cc,
               bool ind1=false);

    /** \brief Erase a submatrix (leaving structural zeros in its place)

        Erase elements of a matrix

        \identifier{19g} */
    void erase(const std::vector<casadi_int>& rr, bool ind1=false);

    /** \brief Remove columns and rows

        Remove/delete rows and/or columns of a matrix

        \identifier{19h} */
    void remove(const std::vector<casadi_int>& rr, const std::vector<casadi_int>& cc);

    /** \brief Enlarge matrix

        Make the matrix larger by inserting empty rows and columns,
        keeping the existing non-zeros

        \identifier{19i} */
    void enlarge(casadi_int nrow, casadi_int ncol,
                  const std::vector<casadi_int>& rr, const std::vector<casadi_int>& cc,
                  bool ind1=false);

#ifndef SWIG
    ///@{
    /// Access the non-zero elements
    std::vector<Scalar>& nonzeros();
    const std::vector<Scalar>& nonzeros() const;
    ///@}

    ///@{
    /// Get a pointer to the data
    Scalar* ptr();
    const Scalar* ptr() const;
    friend inline Scalar* get_ptr(Matrix<Scalar>& v) { return v.ptr(); }
    friend inline const Scalar* get_ptr(const Matrix<Scalar>& v) { return v.ptr(); }
    ///@}

    /// Const access the sparsity - reference to data member
    const Sparsity& sparsity() const;

#endif // SWIG

    /** \brief Get an owning reference to the sparsity pattern

        \identifier{19j} */
    Sparsity get_sparsity() const;

    /** \brief  Get an output

        \identifier{28i} */
    Matrix<Scalar> get_output(casadi_int oind) const;

    /** \brief Construct a sparse matrix from triplet form

     * Default matrix size is max(col) x max(row)

        \identifier{23t} */
    ///@{
    static Matrix<Scalar> triplet(const std::vector<casadi_int>& row,
                                  const std::vector<casadi_int>& col,
                                  const Matrix<Scalar>& d);
    static Matrix<Scalar> triplet(const std::vector<casadi_int>& row,
                                  const std::vector<casadi_int>& col,
                                  const Matrix<Scalar>& d, casadi_int nrow, casadi_int ncol);
    static Matrix<Scalar> triplet(const std::vector<casadi_int>& row,
                                  const std::vector<casadi_int>& col,
                                  const Matrix<Scalar>& d,
                                  const std::pair<casadi_int, casadi_int>& rc);
    ///@}

    ///@{
    /** \brief  create a matrix with all inf

        \identifier{19k} */
    static Matrix<Scalar> inf(const Sparsity& sp);
    static Matrix<Scalar> inf(casadi_int nrow=1, casadi_int ncol=1);
    static Matrix<Scalar> inf(const std::pair<casadi_int, casadi_int>& rc);
    ///@}

    ///@{
    /** \brief  create a matrix with all nan

        \identifier{19l} */
    static Matrix<Scalar> nan(const Sparsity& sp);
    static Matrix<Scalar> nan(casadi_int nrow=1, casadi_int ncol=1);
    static Matrix<Scalar> nan(const std::pair<casadi_int, casadi_int>& rc);
    ///@}

    /** \brief  create an n-by-n identity matrix

        \identifier{19m} */
    static Matrix<Scalar> eye(casadi_int n);

    /** \brief Returns a number that is unique for a given symbolic scalar
     *
     * Only defined if symbolic scalar.

        \identifier{19n} */
    casadi_int element_hash() const;

    /// Checks if expression does not contain NaN or Inf
    bool is_regular() const;

    /** \brief Check if smooth

        \identifier{19o} */
    bool is_smooth() const;

    /** \brief Check if SX is a leaf of the SX graph

        Only defined if symbolic scalar.

        \identifier{19p} */
    bool is_leaf() const;

    /** \brief Check whether a binary SX is commutative

        Only defined if symbolic scalar.

        \identifier{19q} */
    bool is_commutative() const;

    /** \brief Check if symbolic (Dense)

        Check if an expression is a pure symbol.
        Pure means that no operations should happen on the symbol
        (not even vec, transpose, index, concatenation, ...)
        By consequence, a slice of a vector-shaped MX symbol is not a symbol.
        However, the SX type is really a container format with scalar entries.
        Therefore, a slice of a vector-shaped SX symbol is still a symbol.

        Sparse matrices invariable return false

        \seealso is_valid_input

        \identifier{19r} */
    bool is_symbolic() const;

    /** \brief Check if matrix can be used to define function inputs.


        is_valid_input is more forgiving than is_symbolic.
        Some compositions are allowed: vec, vertcat.

        Sparse matrices can return true if all non-zero elements are symbolic

        \seealso is_symbolic

        \identifier{19s} */
    bool is_valid_input() const;

    /// \cond INTERNAL
    /** \brief Detect duplicate symbolic expressions

        If there are symbolic primitives appearing more than once, the function will return
        true and the names of the duplicate expressions will be passed to casadi_warning.
        Note: Will mark the node using SXElem::set_temp.
        Make sure to call reset_input() after usage.

        \identifier{19t} */
    bool has_duplicates() const;

    /** \brief Reset the marker for an input expression

        \identifier{19u} */
    void reset_input() const;
  /// \endcond

    /** \brief Check if the matrix is constant (note that false negative answers are possible)

        \identifier{19v} */
    bool is_constant() const;

    /** \brief Check if function call

        \identifier{28j} */
    bool is_call() const;

    /** \brief  Check if evaluation output

        \identifier{28k} */
    bool is_output() const;

    /** \brief  Check if a multiple output node

        \identifier{28l} */
    bool has_output() const;

    /** \brief Get the index of evaluation output - only valid when is_output() is true

        \identifier{28m} */
    casadi_int which_output() const;

    /** \brief Get function - only valid when is_call() is true

        \identifier{28n} */
    Function which_function() const;

    /** \brief Check if the matrix is integer-valued

     * (note that false negative answers are possible)

        \identifier{19w} */
    bool is_integer() const;

    /** \brief  check if the matrix is 0 (note that false negative answers are possible)

        \identifier{19x} */
    bool is_zero() const;

    /** \brief  check if the matrix is 1 (note that false negative answers are possible)

        \identifier{19y} */
    bool is_one() const;

    /** \brief  check if the matrix is -1 (note that false negative answers are possible)

        \identifier{19z} */
    bool is_minus_one() const;

    /** \brief  check if the matrix is an identity matrix (note that false negative answers

     * are possible)

        \identifier{1a0} */
    bool is_eye() const;

    /// Get operation type
    casadi_int op() const;

    /// Is it a certain operation
    bool is_op(casadi_int op) const;

    /** \brief  Check if the matrix has any zero entries which are not structural zeros

        \identifier{1a1} */
    bool has_zeros() const;

    /** \brief Get all nonzeros

        \identifier{1a2} */
    std::vector<Scalar> get_nonzeros() const;

    /** \brief Get all elements

        \identifier{1a3} */
    std::vector<Scalar> get_elements() const;

#ifndef SWIG
    /** \brief Get all nonzeros

        \identifier{1a4} */
    template<typename A>
    std::vector<A> get_nonzeros() const;
#endif // SWIG

    /** \brief Type conversion to double

        \identifier{1a5} */
    explicit operator double() const;

    /** \brief Type conversion to casadi_int

        \identifier{1a6} */
    explicit operator casadi_int() const;

#ifndef SWIG
    /** \brief Type conversion to a vector

        \identifier{1a7} */
    template<typename A>
    explicit operator std::vector<A>() const;
#endif // SWIG

    /** \brief Get name (only if symbolic scalar)

        \identifier{1a8} */
    std::string name() const;

    /** \brief Get expressions of the children of the expression

        Only defined if symbolic scalar.
        Wraps SXElem SXElem::dep(casadi_int ch=0) const.

        \identifier{1a9} */
    Matrix<Scalar> dep(casadi_int ch=0) const;

    /** \brief Get the number of dependencies of a binary SXElem

        Only defined if symbolic scalar.

        \identifier{1aa} */
    casadi_int n_dep() const;

    // @{
    /// Set the 'precision, width & scientific' used in printing and serializing to streams
    static void set_precision(casadi_int precision);
    static void set_width(casadi_int width);
    static void set_scientific(bool scientific);
    // @}

    /// Seed the random number generator
    static void rng(casadi_int seed);

    ///@{
    /** \brief Create a matrix with uniformly distributed random numbers

        \identifier{1ab} */
    static Matrix<Scalar> rand(  // NOLINT(runtime/threadsafe_fn)
        casadi_int nrow=1,
        casadi_int ncol=1);
    static Matrix<Scalar> rand(const Sparsity& sp);  // NOLINT(runtime/threadsafe_fn)
    static Matrix<Scalar> rand(  // NOLINT(runtime/threadsafe_fn)
        const std::pair<casadi_int, casadi_int>& rc);
    ///@}

    /** \brief Export matrix in specific language
     *
     * lang: only 'matlab' supported for now
     * \verbatim
     * options:
     *   inline: Indicates if you want everything on a single line (default: False)
     *   name: Name of exported variable (default: 'm')
     *   indent_level: Level of indentation (default: 0)
     *   spoof_zero: Replace numerical zero by a 1e-200 (default: false)
     *               might be needed for matlab sparse construct,
     *               which doesn't allow numerical zero
     * \endverbatim

        \identifier{1ac} */
    void export_code(const std::string& lang,
        std::ostream &stream=casadi::uout(), const Dict& options=Dict()) const;

    /** Obtain information about sparsity */
    Dict info() const;
    #ifndef SWIG
        /** \brief Serialize an object

            \identifier{1ad} */
        void serialize(std::ostream &stream) const;
    #endif

    /** \brief Serialize

        \identifier{1ae} */
    std::string serialize() const;

    /** \brief Build Sparsity from serialization

        \identifier{1af} */
    static Matrix<Scalar> deserialize(std::istream& stream);

    /** \brief Build Sparsity from serialization

        \identifier{1ag} */
    static Matrix<Scalar> deserialize(const std::string& s);

    /** \brief Serialize an object

        \identifier{1ah} */
    void serialize(SerializingStream& s) const;

    static Matrix<Scalar> deserialize(DeserializingStream& s);

    // @{
    /** Export numerical matrix to file
    *
    * Supported formats:
    *
    * \verbatim
    *   - .mtx   Matrix Market (sparse)
    *   - .txt   Ascii full precision representation (sparse)
    *            Whitespace separated, aligned.
    *            Comments with # % or /
    *            Uses C locale
    *            Structural zeros represented by 00
    *            Does not scale well for large sparse matrices
    * \endverbatim
    *
    */
    void to_file(const std::string& filename, const std::string& format="") const;
#ifndef SWIG
    static void to_file(const std::string& filename, const Sparsity& sp,
      const Scalar* nonzeros, const std::string& format="");
#endif

    static Matrix<double> from_file(const std::string& filename, const std::string& format_hint="");
    //@}

#ifndef SWIG
    /// Sparse matrix with a given sparsity with all values same
    Matrix(const Sparsity& sp, const Scalar& val, bool dummy);

    /// Sparse matrix with a given sparsity and non-zero elements.
    Matrix(const Sparsity& sp, const std::vector<Scalar>& d, bool dummy);

    /// \cond INTERNAL
    static Matrix<Scalar> _sym(const std::string& name, const Sparsity& sp);
    /// \endcond

#ifdef CASADI_WITH_THREADSAFE_SYMBOLICS
    static std::mutex& get_mutex_temp();
#endif //CASADI_WITH_THREADSAFE_SYMBOLICS

  private:
    /// Sparsity of the matrix in a compressed column storage (CCS) format
    Sparsity sparsity_;

    /// Nonzero elements
    std::vector<Scalar> nonzeros_;

    /// Precision used in streams
    static casadi_int stream_precision_;
    static casadi_int stream_width_;
    static bool stream_scientific_;

    /// Random number generator
    static std::default_random_engine rng_;

#endif // SWIG
  };

  /// Implementation of Matrix::get_nonzeros (in public API)
  template<typename Scalar>
  template<typename A>
  std::vector<A> Matrix<Scalar>::get_nonzeros() const {
    std::vector<A> ret(nnz());
    auto r = ret.begin();
    for (auto&& e : nonzeros()) *r++ = static_cast<A>(e);
    return ret;
  }

  /// Implementation of std::vector(Matrix) (in public API)
  template<typename Scalar>
  template<typename A>
  Matrix<Scalar>::operator std::vector<A>() const {
    // Get sparsity pattern
    casadi_int size1 = this->size1(), size2 = this->size2();
    const casadi_int *colind = this->colind(), *row = this->row();
    // Copy the nonzeros
    auto it = nonzeros().begin();
    std::vector<A> ret(numel(), 0);
    for (casadi_int cc=0; cc<size2; ++cc) {
      for (casadi_int el=colind[cc]; el<colind[cc+1]; ++el) {
        ret[row[el] + cc*size1] = static_cast<A>(*it++);
      }
    }
    return ret;
  }

} // namespace casadi

#endif // CASADI_MATRIX_DECL_HPP
