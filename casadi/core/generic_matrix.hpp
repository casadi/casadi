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
#include "calculus.hpp"
#include "sparsity_interface.hpp"
#include "generic_type.hpp"

namespace casadi {
  /** \brief Empty Base
      This class is extended in SWIG.
   */
  struct CASADI_EXPORT GenericMatrixCommon {};

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
  class CASADI_EXPORT GenericMatrix
    : public GenericMatrixCommon,
      public SparsityInterface<MatType> {
    using SparsityInterface<MatType>::self;
  public:

    /** \brief Get the number of (structural) non-zero elements */
    int nnz() const;

    /** \brief Get the number of non-zeros in the lower triangular half */
    int nnz_lower() const;

    /** \brief Get the number of non-zeros in the upper triangular half */
    int nnz_upper() const;

    /** \brief Get get the number of non-zeros on the diagonal */
    int nnz_diag() const;

    /** \brief Get the number of elements */
    int numel() const;

    /** \brief Get the number of elements in slice (cf. MATLAB) */
    int numel(int i) const { return 1;}

    /** \brief Get the first dimension (i.e. number of rows) */
    int size1() const;

    /** \brief Get the second dimension (i.e. number of columns) */
    int size2() const;

    /** \brief Get string representation of dimensions.
        The representation is (nrow x ncol = numel | size)
    */
    std::string dim() const;

    /** \brief  Get the shape */
    std::pair<int, int> size() const;

    /** \brief  Get the size along a particular dimensions */
    int size(int axis) const;

    /** \brief Check if the sparsity is empty, i.e. if one of the dimensions is zero
     * (or optionally both dimensions) */
    bool is_empty(bool both=false) const { return sparsity().is_empty(both);}

    /** \brief  Check if the matrix expression is dense */
    bool is_dense() const { return sparsity().is_dense();}

    /** \brief  Check if the matrix expression is scalar */
    bool is_scalar(bool scalar_and_dense=false) const;

    /** \brief  Check if the matrix expression is square */
    bool is_square() const { return sparsity().is_square();}

    /** \brief  Check if the matrix is a row or column vector */
    bool is_vector() const { return sparsity().is_vector();}

    /** \brief  Check if the matrix is a row vector (i.e. size1()==1) */
    bool is_row() const { return sparsity().is_row();}

    /** \brief  Check if the matrix is a column vector (i.e. size2()==1) */
    bool is_column() const { return sparsity().is_column();}

    /** \brief Check if the matrix is upper triangular */
    bool is_triu() const { return sparsity().is_triu();}

    /** \brief Check if the matrix is lower triangular */
    bool is_tril() const { return sparsity().is_tril();}

    ///@{
    /** \brief Get the sparsity pattern. See the Sparsity class for details. */
    std::vector<int> get_row() const { return sparsity().get_row(); }
    std::vector<int> get_colind() const { return sparsity().get_colind(); }
#ifndef SWIG
    const int* row() const { return sparsity().row(); }
    const int* colind() const { return sparsity().colind(); }
#endif
    int row(int el) const { return sparsity().row(el); }
    int colind(int col) const { return sparsity().colind(col); }
    ///@}

    /** \brief Get the location of all non-zero elements as they would appear in a Dense matrix
        A : DenseMatrix  4 x 3
        B : SparseMatrix 4 x 3 , 5 structural non-zeros

        k = A.find()
        A[k] will contain the elements of A that are non-zero in B
    */
    std::vector<int> find(bool ind1=SWIG_IND1) const { return sparsity().find(ind1);}

    /** \brief Get the sparsity pattern */
    SWIG_CONSTREF(Sparsity) sparsity() const;

#ifndef SWIG
    /// \cond INTERNAL
    /** \brief Access the sparsity, make a copy if there are multiple references to it */
    Sparsity& sparsityRef();
    /// \endcond


    /// \cond CLUTTER
    /**  @{  */
    /** \brief Accessed by friend functions */
    static int sprank(const MatType &x) { return Sparsity::sprank(x.sparsity());}
    static int norm_0_mul(const MatType &x, const MatType &y) {
      return Sparsity::norm_0_mul(x.sparsity(), y.sparsity());
    }
    static MatType tril(const MatType &x, bool includeDiagonal=true) {
      return project(x, Sparsity::tril(x.sparsity(), includeDiagonal));
    }
    static MatType triu(const MatType &x, bool includeDiagonal=true) {
      return project(x, Sparsity::triu(x.sparsity(), includeDiagonal));
    }
    static MatType sum_square(const MatType &x) { return dot(x, x);}
    static MatType linspace(const MatType &a, const MatType &b, int nsteps);
    static MatType cross(const MatType &a, const MatType &b, int dim=-1);
    static MatType skew(const MatType &a);
    static MatType inv_skew(const MatType &a);
    static MatType tril2symm(const MatType &x);
    static MatType triu2symm(const MatType &x);
    static MatType repsum(const MatType &x, int n, int m=1);
    /** @}  */
    /// \endcond

    /** \brief  Get vector nonzero or slice of nonzeros */
    template<typename K>
    const MatType operator[](const K& k) const {
      MatType ret;
      self().get_nz(ret, false, k);
      return ret;
    }

    /** \brief  Access vector nonzero or slice of nonzeros */
    template<typename K>
    NonZeros<MatType, K> operator[](const K& k) {
      return NonZeros<MatType, K>(self(), k);
    }

    /** \brief  Get vector element or slice */
    template<typename RR>
    const MatType operator()(const RR& rr) const {
      MatType ret;
      self().get(ret, false, rr);
      return ret;
    }

    /** \brief  Get Matrix element or slice */
    template<typename RR, typename CC>
    const MatType operator()(const RR& rr, const CC& cc) const {
      MatType ret;
      self().get(ret, false, rr, cc);
      return ret;
    }

    /** \brief Access Matrix elements (one argument) */
    template<typename RR>
    SubIndex<MatType, RR> operator()(const RR& rr) {
      return SubIndex<MatType, RR>(self(), rr);
    }

    /** \brief Access Matrix elements (two arguments) */
    template<typename RR, typename CC>
    SubMatrix<MatType, RR, CC> operator()(const RR& rr, const CC& cc) {
      return SubMatrix<MatType, RR, CC>(self(), rr, cc);
    }
#endif // SWIG

#if !defined(SWIG) || defined(DOXYGEN)
/**
\ingroup expression_tools
@{
*/
    /** \brief Matrix power x^n
     */
    inline friend MatType mpower(const MatType& x, const MatType& n) {
      return MatType::mpower(x, n);
    }

    /** \brief Matrix divide (cf. slash '/' in MATLAB)
     */
    inline friend MatType mrdivide(const MatType& x, const MatType& n) {
      return MatType::mrdivide(x, n);
    }

    /** \brief Matrix divide (cf. backslash '\' in MATLAB)
     */
    inline friend MatType mldivide(const MatType& x, const MatType& n) {
      return MatType::mldivide(x, n);
    }

    /** \brief Get all symbols contained in the supplied expression
     * Get all symbols on which the supplied expression depends
     * \see SXFunction::getFree(), MXFunction::getFree()
     */
    inline friend std::vector<MatType> symvar(const MatType& x) {
      return MatType::symvar(x);
    }

    ///@{
    /** \brief Calculate bilinear form x^T A y
     */
    inline friend MatType bilin(const MatType &A, const MatType &x, const MatType &y) {
      return MatType::bilin(A, x, y);
    }
    static MatType bilin(const MatType& A, const MatType& x, const MatType& y);
    ///@}

    ///@{
    /** \brief Make a rank-1 update to a matrix A
     * Calculates A + 1/2 * alpha * x*y'
     */
    inline friend MatType rank1(const MatType &A, const MatType &alpha,
                                const MatType &x, const MatType &y) {
      return MatType::rank1(A, alpha, x, y);
    }
    static MatType rank1(const MatType& A, const MatType& alpha,
                         const MatType& x, const MatType& y);
    ///@}

    /** \brief Calculate some of squares: sum_ij X_ij^2
     */
    inline friend MatType sum_square(const MatType &x) {
      return MatType::sum_square(x);
    }

    /** \brief Matlab's \c linspace command
     */
    inline friend MatType linspace(const MatType &a, const MatType &b, int nsteps) {
      return MatType::linspace(a, b, nsteps);
    }

    /** \brief Matlab's \c cross command
     */
    inline friend MatType cross(const MatType &a, const MatType &b, int dim = -1) {
      return MatType::cross(a, b, dim);
    }

    /** \brief Generate a skew symmetric matrix from a 3-vector
     */
    inline friend MatType skew(const MatType &a) {
      return MatType::skew(a);
    }

    /** \brief Generate the 3-vector progenitor of a skew symmetric matrix
     */
    inline friend MatType inv_skew(const MatType &a) {
      return MatType::inv_skew(a);
    }

    /** \brief Matrix determinant (experimental) */
    inline friend MatType det(const MatType& A) { return MatType::det(A);}

    /** \brief Matrix inverse (experimental) */
    inline friend MatType inv(const MatType& A) { return MatType::inv(A);}

    /** \brief Matrix trace */
    inline friend MatType trace(const MatType& x) { return MatType::trace(x);}

    /** \brief Convert a lower triangular matrix to a symmetric one
     */
    inline friend MatType tril2symm(const MatType &a) { return MatType::tril2symm(a);}

    /** \brief Convert a upper triangular matrix to a symmetric one
     */
    inline friend MatType triu2symm(const MatType &a) { return MatType::triu2symm(a);}

    /** \brief  Frobenius norm  */
    inline friend MatType norm_F(const MatType &x) { return MatType::norm_F(x);}

    /** \brief  2-norm  */
    inline friend MatType norm_2(const MatType &x) { return MatType::norm_2(x);}

    /** \brief 1-norm  */
    inline friend MatType norm_1(const MatType &x) { return MatType::norm_1(x);}

    /** \brief Infinity-norm */
    inline friend MatType norm_inf(const MatType &x) { return MatType::norm_inf(x);}

    /** \brief Return a row-wise summation of elements */
    inline friend MatType sum1(const MatType &x) { return MatType::sum1(x);}

    /** \brief Return a column-wise summation of elements */
    inline friend MatType sum2(const MatType &x) { return MatType::sum2(x);}

    /** \brief Inner product of two matrices
        with x and y matrices of the same dimension
    */
    inline friend MatType dot(const MatType &x, const MatType &y) {
      return MatType::dot(x, y);
    }

    /** \brief Computes the nullspace of a matrix A
     *
     * Finds Z m-by-(m-n) such that AZ = 0
     * with A n-by-m with m > n
     *
     * Assumes A is full rank
     *
     * Inspired by Numerical Methods in Scientific Computing by Ake Bjorck
     */
    inline friend MatType nullspace(const MatType& A) {
      return MatType::nullspace(A);
    }

    /** \brief  Evaluate a polynomial with coefficients p in x
    */
    inline friend MatType polyval(const MatType& p, const MatType& x) {
      return MatType::polyval(p, x);
    }

    /** \brief   Get the diagonal of a matrix or construct a diagonal
        When the input is square, the diagonal elements are returned.
        If the input is vector-like, a diagonal matrix is constructed with it. */
    inline friend MatType diag(const MatType &A) {
      return MatType::diag(A);
    }

    /** \brief  Unite two matrices no overlapping sparsity
     */
    inline friend MatType unite(const MatType& A, const MatType& B) {
      return MatType::unite(A, B);
    }

    /** \brief  Make the matrix dense if not already
     */
    inline friend MatType densify(const MatType& x) {
      return MatType::densify(x);
    }

    /** \brief  Make the matrix dense and assign nonzeros to a value
     */
    inline friend MatType densify(const MatType& x, const MatType& val) {
      return MatType::densify(x, val);
    }

    /** \brief Create a new matrix with a given sparsity pattern but with the
      * nonzeros taken from an existing matrix
      */
    inline friend MatType project(const MatType& A, const Sparsity& sp,
                                  bool intersect=false) {
      return MatType::project(A, sp, intersect);
    }

    /** \brief Branching on MX nodes
        Ternary operator, "cond ? if_true : if_false"
    */
    inline friend MatType if_else(const MatType &cond, const MatType &if_true,
                                  const MatType &if_false, bool short_circuit=true) {
      return MatType::if_else(cond, if_true, if_false, short_circuit);
    }

    /** \brief Create a switch
     *
     * If the condition \param ind evaluates to the integer k, where 0<=k<f.size(),
     * then x[k] will be returned, otherwise \param x_default will be returned.
     */
    inline friend MatType conditional(const MatType& ind, const std::vector<MatType> &x,
                                      const MatType &x_default, bool short_circuit=true) {
      return MatType::conditional(ind, x, x_default, short_circuit);
    }

    /** \brief Check if expression depends on the argument
        The argument must be symbolic
    */
    inline friend bool depends_on(const MatType& f, const MatType &arg) {
      return MatType::depends_on(f, arg);
    }

    /** \brief  Substitute variable v with expression vdef in an expression ex */
    friend inline MatType substitute(const MatType& ex, const MatType& v,
                                     const MatType& vdef) {
      return MatType::substitute(ex, v, vdef);
    }

    /** \brief  Substitute variable var with expression expr in multiple expressions */
    friend inline std::vector<MatType>
      substitute(const std::vector<MatType>& ex, const std::vector<MatType>& v,
                 const std::vector<MatType>& vdef) {
      return MatType::substitute(ex, v, vdef);
    }

    /** \brief Inplace substitution with piggyback expressions
     * Substitute variables v out of the expressions vdef sequentially,
     * as well as out of a number of other expressions piggyback */
    inline friend void
      substitute_inplace(const std::vector<MatType>& v,
                        std::vector<MatType>& inout_vdef,
                        std::vector<MatType>& inout_ex, bool reverse=false) {
      return MatType::substitute_inplace(v, inout_vdef, inout_ex, reverse);
    }

    /** \brief  Solve a system of equations: A*x = b
        The solve routine works similar to Matlab's backslash when A is square and nonsingular.
        The algorithm used is the following:
        1. A simple forward or backward substitution if A is upper or lower triangular
        2. If the linear system is at most 3-by-3, form the inverse via minor expansion and multiply
        3. Permute the variables and equations as to get a (structurally) nonzero diagonal,
        then perform a QR factorization without pivoting and solve the factorized system.

        Note 1: If there are entries of the linear system known to be zero, these will be removed.
        Elements that are very small, or will evaluate to be zero, can still cause numerical errors,
        due to the lack of pivoting (which is not possible since cannot compare the size of entries)

        Note 2: When permuting the linear system, a BLT (block lower triangular) transformation is
        formed. Only the permutation part of this is however used. An improvement would be to solve
        block-by-block if there are multiple BLT blocks.
    */
    friend inline MatType solve(const MatType& A, const MatType& b) {
      return MatType::solve(A, b);
    }

    /** \brief Solve a system of equations: A*x = b
     */
    friend inline MatType solve(const MatType& A, const MatType& b,
                                const std::string& lsolver,
                                const Dict& dict = Dict()) {
      return MatType::solve(A, b, lsolver, dict);
    }

    /** \brief Computes the Moore-Penrose pseudo-inverse
     *
     * If the matrix A is fat (size1<size2), mul(A, pinv(A)) is unity.
     *
     *  pinv(A)' = (AA')^(-1) A
     *
     *
     * If the matrix A is slender (size1>size2), mul(pinv(A), A) is unity.
     *
     *  pinv(A) = (A'A)^(-1) A'
     *
     */
    friend inline MatType pinv(const MatType& A) {
      return MatType::pinv(A);
    }

    /** \brief Computes the Moore-Penrose pseudo-inverse
     *
     * If the matrix A is fat (size1>size2), mul(A, pinv(A)) is unity.
     * If the matrix A is slender (size2<size1), mul(pinv(A), A) is unity.
     *
     */
    friend inline MatType pinv(const MatType& A, const std::string& lsolver,
                               const Dict& dict = Dict()) {
      return MatType::pinv(A, lsolver, dict);
    }

    ///@{
    /** \brief Calculate Jacobian
    */
    inline friend MatType jacobian(const MatType &ex, const MatType &arg, bool symmetric=false) {
      return MatType::jacobian(ex, arg, symmetric);
    }
    inline friend MatType gradient(const MatType &ex, const MatType &arg) {
      return MatType::gradient(ex, arg);
    }
    inline friend MatType tangent(const MatType &ex, const MatType &arg) {
      return MatType::tangent(ex, arg);
    }
    ///@}

    /** \brief Calculate the Jacobian and multiply by a vector from the right
        This is equivalent to <tt>mul(jacobian(ex, arg), v)</tt> or
        <tt>mul(jacobian(ex, arg).T, v)</tt> for
        tr set to false and true respectively. If contrast to these
        expressions, it will use directional derivatives which is typically (but
        not necessarily) more efficient if the complete Jacobian is not needed and v has few rows.
    */
    friend inline MatType jtimes(const MatType &ex, const MatType &arg,
                                 const MatType &v, bool tr=false) {
      return MatType::jtimes(ex, arg, v, tr);
    }

    ///@{
    // Hessian and (optionally) gradient
    inline friend MatType hessian(const MatType &ex, const MatType &arg) {
      return MatType::hessian(ex, arg);
    }
    inline friend MatType hessian(const MatType &ex, const MatType &arg, MatType& output_g) {
      return MatType::hessian(ex, arg, output_g);
    }
    ///@}

    ///@{
    /** \brief Find out which variables enter nonlinearly */
    inline friend std::vector<bool> nl_var(const MatType &expr, const MatType &var) {
      return MatType::nl_var(expr, var);
    }
    ///@}

    /** Count number of nodes */
    inline friend int n_nodes(const MatType& A) {
      return MatType::n_nodes(A);
    }

    /// Simplify an expression
    friend inline MatType simplify(const MatType &x) {
      return MatType::simplify(x);
    }

    /** \brief Get a string representation for a binary MatType, using custom arguments */
    inline friend std::string
      print_operator(const MatType& xb, const std::vector<std::string>& args) {
      return MatType::print_operator(xb, args);
    }

    /** \brief Extract shared subexpressions from an set of expressions */
    inline friend void shared(std::vector<MatType>& ex,
                                     std::vector<MatType>& v,
                                     std::vector<MatType>& vdef,
                                     const std::string& v_prefix="v_",
                                     const std::string& v_suffix="") {
      MatType::shared(ex, v, vdef, v_prefix, v_suffix);
    }

    /** \brief Extract shared subexpressions from an set of expressions */
    inline friend void shared(const std::vector<MatType>& ex,
                                     std::vector<MatType>& ex_output,
                                     std::vector<MatType>& v,
                                     std::vector<MatType>& vdef,
                                     const std::string& v_prefix="v_",
                                     const std::string& v_suffix="") {
      ex_output = ex;
      shared(ex_output, v, vdef, v_prefix, v_suffix);
    }

    /** \brief Given a repeated matrix, computes the sum of repeated parts
     */
    inline friend MatType repsum(const MatType &A, int n, int m=1) {
      return MatType::repsum(A, n, m);
    }


/** @} */
#endif // SWIG

    /** @name Construct symbolic primitives
        The "sym" function is intended to work in a similar way as "sym" used
        in the Symbolic Toolbox for Matlab but instead creating a
        CasADi symbolic primitive.
    */
    ///@{

    /** \brief Create an nrow-by-ncol symbolic primitive */
    static MatType sym(const std::string& name, int nrow=1, int ncol=1) {
      return sym(name, Sparsity::dense(nrow, ncol));
    }

    /** \brief  Construct a symbolic primitive with given dimensions */
    static MatType sym(const std::string& name, const std::pair<int, int> &rc) {
      return sym(name, rc.first, rc.second);
    }

    /** \brief Create symbolic primitive with a given sparsity pattern */
    static MatType sym(const std::string& name, const Sparsity& sp);

    /** \brief Create a vector of length p with with matrices
     * with symbolic primitives of given sparsity */
    static std::vector<MatType > sym(const std::string& name, const Sparsity& sp, int p);

    /** \brief Create a vector of length p with nrow-by-ncol symbolic primitives */
    static std::vector<MatType > sym(const std::string& name, int nrow, int ncol, int p) {
      return sym(name, Sparsity::dense(nrow, ncol), p);
    }

    /** \brief Create a vector of length r of vectors of length p with
     * symbolic primitives with given sparsity*/
    static std::vector<std::vector<MatType> >
      sym(const std::string& name, const Sparsity& sp, int p, int r);

    /** \brief Create a vector of length r of vectors of length p
     * with nrow-by-ncol symbolic primitives */
    static std::vector<std::vector<MatType> >
      sym(const std::string& name, int nrow, int ncol, int p, int r) {
      return sym(name, Sparsity::dense(nrow, ncol), p, r);
    }
    ///@}

    ///@{
    /** \brief Create a dense matrix or a matrix with specified sparsity with all entries zero */
    static MatType zeros(int nrow=1, int ncol=1) { return zeros(Sparsity::dense(nrow, ncol)); }
    static MatType zeros(const Sparsity& sp) { return MatType(sp, 0, false);}
    static MatType zeros(const std::pair<int, int>& rc) { return zeros(rc.first, rc.second);}
    ///@}

    ///@{
    /** \brief Create a dense matrix or a matrix with specified sparsity with all entries one */
    static MatType ones(int nrow=1, int ncol=1) { return ones(Sparsity::dense(nrow, ncol)); }
    static MatType ones(const Sparsity& sp) { return MatType(sp, 1, false);}
    static MatType ones(const std::pair<int, int>& rc) { return ones(rc.first, rc.second);}
    ///@}
  };

#ifndef SWIG
  // Implementations

  template<typename MatType>
  const Sparsity& GenericMatrix<MatType>::sparsity() const {
    return self().sparsity();
  }

  template<typename MatType>
  Sparsity& GenericMatrix<MatType>::sparsityRef() {
    return self().sparsityRef();
  }

  template<typename MatType>
  int GenericMatrix<MatType>::nnz() const {
    return sparsity().nnz();
  }

  template<typename MatType>
  int GenericMatrix<MatType>::nnz_lower() const {
    return sparsity().nnz_lower();
  }

  template<typename MatType>
  int GenericMatrix<MatType>::nnz_upper() const {
    return sparsity().nnz_upper();
  }

  template<typename MatType>
  int GenericMatrix<MatType>::nnz_diag() const {
    return sparsity().nnz_diag();
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
  std::pair<int, int> GenericMatrix<MatType>::size() const {
    return sparsity().size();
  }

  template<typename MatType>
  int GenericMatrix<MatType>::size(int axis) const {
    return sparsity().size(axis);
  }

  template<typename MatType>
  std::string GenericMatrix<MatType>::dim() const {
    return sparsity().dim();
  }

  template<typename MatType>
  bool GenericMatrix<MatType>::is_scalar(bool scalar_and_dense) const {
    return sparsity().is_scalar(scalar_and_dense);
  }

#endif // SWIG

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

  template<typename MatType>
  MatType GenericMatrix<MatType>::linspace(const MatType& a, const MatType& b, int nsteps) {
    std::vector<MatType> ret(nsteps);
    ret[0] = a;
    MatType step = (b-a)/(nsteps-1);

    for (int i=1; i<nsteps-1; ++i)
      ret[i] = ret[i-1] + step;

    ret[nsteps-1] = b;
    return vertcat(ret);
  }

  template<typename MatType>
  MatType GenericMatrix<MatType>::cross(const MatType& a, const MatType& b, int dim) {
    casadi_assert_message(a.size1()==b.size1() && a.size2()==b.size2(),
                          "cross(a, b): Inconsistent dimensions. Dimension of a ("
                          << a.dim() << " ) must equal that of b ("
                          << b.dim() << ").");

    casadi_assert_message(a.size1()==3 || a.size2()==3,
                          "cross(a, b): One of the dimensions of a should have length 3, but got "
                          << a.dim() << ".");
    casadi_assert_message(dim==-1 || dim==1 || dim==2,
                          "cross(a, b, dim): Dim must be 1, 2 or -1 (automatic).");

    std::vector<MatType> ret(3);

    bool t = a.size1()==3;

    if (dim==1) t = true;
    if (dim==2) t = false;

    MatType a1 = t ? a(0, Slice()) : a(Slice(), 0);
    MatType a2 = t ? a(1, Slice()) : a(Slice(), 1);
    MatType a3 = t ? a(2, Slice()) : a(Slice(), 2);

    MatType b1 = t ? b(0, Slice()) : b(Slice(), 0);
    MatType b2 = t ? b(1, Slice()) : b(Slice(), 1);
    MatType b3 = t ? b(2, Slice()) : b(Slice(), 2);

    ret[0] = a2*b3-a3*b2;
    ret[1] = a3*b1-a1*b3;
    ret[2] = a1*b2-a2*b1;

    return t ? vertcat(ret) : horzcat(ret);
  }

  template<typename MatType>
  MatType GenericMatrix<MatType>::skew(const MatType& a) {
    casadi_assert_message(a.is_vector() && (a.size1()==3 || a.size2()==3),
                          "skew(a): Expecting 3-vector, got " << a.dim() << ".");

    MatType x = a(0);
    MatType y = a(1);
    MatType z = a(2);
    return blockcat(std::vector< std::vector<MatType> >({{0, -z, y}, {z, 0, -x}, {-y, x, 0}}));
  }

  template<typename MatType>
  MatType GenericMatrix<MatType>::inv_skew(const MatType& a) {
    casadi_assert_message(a.size1()==3 && a.size2()==3,
                          "inv_skew(a): Expecting 3-by-3 matrix, got " << a.dim() << ".");

    return 0.5*vertcat(std::vector<MatType>({a(2, 1)-a(1, 2), a(0, 2)-a(2, 0), a(1, 0)-a(0, 1)}));
  }


  template<typename MatType>
  MatType GenericMatrix<MatType>::tril2symm(const MatType& x) {
    casadi_assert_message(x.is_square(),
                          "Shape error in tril2symm. Expecting square shape but got "
                          << x.dim());
    casadi_assert_message(x.nnz_upper()-x.nnz_diag()==0,
                          "Sparsity error in tril2symm. Found above-diagonal entries in argument: "
                          << x.dim());
    return x +  x.T() - diag(diag(x));
  }

  template<typename MatType>
  MatType GenericMatrix<MatType>::repsum(const MatType& x, int n, int m) {
    casadi_assert(x.size1() % n==0);
    casadi_assert(x.size2() % m==0);
    std::vector< std::vector< MatType> > s =
      blocksplit(x, x.size1()/n, x.size2()/m);
    MatType sum = 0;
    for (int i=0;i<s.size();++i) {
      for (int j=0;j<s[i].size();++j) {
        sum = sum + s[i][j];
      }
    }
    return sum;
  }

  template<typename MatType>
  MatType GenericMatrix<MatType>::triu2symm(const MatType& x) {
    casadi_assert_message(x.is_square(),
                          "Shape error in triu2symm. Expecting square shape but got "
                          << x.dim());
    casadi_assert_message(x.nnz_lower()-x.nnz_diag()==0,
                          "Sparsity error in triu2symm. Found below-diagonal entries in argument: "
                          << x.dim());
    return x + x.T() - diag(diag(x));
  }

  template<typename MatType>
  MatType GenericMatrix<MatType>::bilin(const MatType& A, const MatType& x,
                                        const MatType& y) {
    // Check/correct x
    casadi_assert(x.is_vector());
    if (!x.is_column()) return bilin(A, x.T(), y);
    if (!x.is_dense()) return bilin(A, densify(x), y);

    // Check/correct y
    casadi_assert(y.is_vector());
    if (!y.is_column()) return bilin(A, x, y.T());
    if (!y.is_dense()) return bilin(A, x, densify(y));

    // Assert dimensions
    casadi_assert_message(x.size1()==A.size1() && y.size1()==A.size2(),
                          "Dimension mismatch. Got x.size1() = " << x.size1()
                          << " and y.size1() = " << y.size1()
                          << " but A.size() = " << A.size());

    // Call the class specific method
    return MatType::_bilin(A, x, y);
  }

  template<typename MatType>
  MatType GenericMatrix<MatType>::rank1(const MatType& A, const MatType& alpha,
                                        const MatType& x, const MatType& y) {
    // Check/correct x
    casadi_assert(x.is_vector());
    if (!x.is_column()) return rank1(A, alpha, x.T(), y);
    if (!x.is_dense()) return rank1(A, alpha, densify(x), y);

    // Check/correct y
    casadi_assert(y.is_vector());
    if (!y.is_column()) return rank1(A, alpha, x, y.T());
    if (!y.is_dense()) return rank1(A, alpha, x, densify(y));

    // Check alpha, quick return
    casadi_assert(alpha.is_scalar());
    if (!alpha.is_dense()) return A;

    // Assert dimensions
    casadi_assert_message(x.size1()==A.size1() && y.size1()==A.size2(),
                          "Dimension mismatch. Got x.size1() = " << x.size1()
                          << " and y.size1() = " << y.size1()
                          << " but A.size() = " << A.size());

    // Call the class specific method
    return MatType::_rank1(A, alpha, x, y);
  }

} // namespace casadi

#endif // CASADI_GENERIC_MATRIX_HPP
