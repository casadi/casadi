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
#include "../casadi_math.hpp"
#include "sparsity_interface.hpp"
#include "../generic_type.hpp"

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

    /** \brief [DEPRECATED] Alias for nnz
        The meaning of this function will change in future revision.
        \see nnz()
    */
    int size() const { return nnz();}

    /** \brief Get the number of non-zeros in the lower triangular half */
    int sizeL() const;

    /** \brief Get the number of non-zeros in the upper triangular half */
    int sizeU() const;

    /** \brief Get get the number of non-zeros on the diagonal */
    int sizeD() const;

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
    std::string dimString() const;

    /** \brief  Get the shape */
    std::pair<int, int> shape() const;

    /** \brief  Get the size along a particular dimensions */
    int shape(int axis) const;

    /** \brief Check if the sparsity is empty, i.e. if one of the dimensions is zero
     * (or optionally both dimensions) */
    bool isempty(bool both=false) const { return sparsity().isempty(both);}

    /** \brief  Check if the matrix expression is dense */
    bool isdense() const { return sparsity().isdense();}

    /** \brief  Check if the matrix expression is scalar */
    bool isscalar(bool scalar_and_dense=false) const;

    /** \brief  Check if the matrix expression is square */
    bool issquare() const { return sparsity().issquare();}

    /** \brief  Check if the matrix is a row or column vector */
    bool isvector() const { return sparsity().isvector();}

    /** \brief  Check if the matrix is a row vector (i.e. size1()==1) */
    bool isrow() const { return sparsity().isrow();}

    /** \brief  Check if the matrix is a column vector (i.e. size2()==1) */
    bool iscolumn() const { return sparsity().iscolumn();}

    /** \brief Check if the matrix is upper triangular */
    bool istriu() const { return sparsity().istriu();}

    /** \brief Check if the matrix is lower triangular */
    bool istril() const { return sparsity().istril();}

    ///@{
    /** \brief Get the sparsity pattern. See the Sparsity class for details. */
    std::vector<int> getRow() const { return sparsity().getRow(); }
    std::vector<int> getColind() const { return sparsity().getColind(); }
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
    int zz_sprank() const { return sprank(sparsity());}
    int zz_norm_0_mul(const MatType &A) const { return norm_0_mul(sparsity(), A.sparsity());}
    MatType zz_tril(bool includeDiagonal=true) const {
      return project(self(), tril(sparsity(), includeDiagonal));
    }
    MatType zz_triu(bool includeDiagonal=true) const {
      return project(self(), triu(sparsity(), includeDiagonal));
    }
    MatType zz_quad_form(const MatType &A) const {
      casadi_assert(isvector());
      if (!iscolumn()) return quad_form(self().T(), A);
      return inner_prod(self(), mul(A, self()));
    }
    MatType zz_quad_form() const {
      casadi_assert(isvector());
      return inner_prod(self(), self());
    }
    MatType zz_sum_square() const { return inner_prod(self(), self());}
    MatType zz_linspace(const MatType &b, int nsteps) const;
    MatType zz_cross(const MatType &b, int dim=-1) const;
    MatType zz_tril2symm() const;
    MatType zz_triu2symm() const;
    MatType zz_densify() const {
      MatType ret = self();
      ret.makeDense();
      return ret;
    }
    MatType zz_repsum(int n, int m=1) const;
    /** @}  */
    /// \endcond

    /** \brief  Get vector nonzero or slice of nonzeros */
    template<typename K>
    const MatType operator[](const K& k) const {
      MatType ret;
      self().getNZ(ret, false, k);
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
      return x.zz_mpower(n);
    }

    /** \brief Matrix divide (cf. slash '/' in MATLAB)
     */
    inline friend MatType mrdivide(const MatType& x, const MatType& n) {
      return x.zz_mrdivide(n);
    }

    /** \brief Matrix divide (cf. backslash '\' in MATLAB)
     */
    inline friend MatType mldivide(const MatType& x, const MatType& n) {
      return x.zz_mldivide(n);
    }

    /** \brief Get all symbols contained in the supplied expression
     * Get all symbols on which the supplied expression depends
     * \see SXFunction::getFree(), MXFunction::getFree()
     */
    inline friend std::vector<MatType> symvar(const MatType& x) {
      return x.zz_symvar();
    }

    /** \brief Calculate quadratic form X^T A X
     */
    inline friend MatType quad_form(const MatType &X, const MatType &A) {
      return X.zz_quad_form(A);
    }

    /** \brief Calculate quadratic form X^T X
     */
    inline friend MatType quad_form(const MatType &X) {
      return X.zz_quad_form();
    }

    /** \brief Calculate some of squares: sum_ij X_ij^2 
     */
    inline friend MatType sum_square(const MatType &X) {
      return X.zz_sum_square();
    }

    /** \brief Matlab's \c linspace command
     */
    inline friend MatType linspace(const MatType &a, const MatType &b, int nsteps) {
      return a.zz_linspace(b, nsteps);
    }

    /** \brief Matlab's \c cross command
     */
    inline friend MatType cross(const MatType &a, const MatType &b, int dim = -1) {
      return a.zz_cross(b, dim);
    }

    /** \brief Matrix determinant (experimental) */
    inline friend MatType det(const MatType& A) { return A.zz_det();}

    /** \brief Matrix inverse (experimental) */
    inline friend MatType inv(const MatType& A) { return A.zz_inv();}

    /** \brief Matrix trace */
    inline friend MatType trace(const MatType& a) { return a.zz_trace();}

    /** \brief Convert a lower triangular matrix to a symmetric one
     */
    inline friend MatType tril2symm(const MatType &a) { return a.zz_tril2symm();}

    /** \brief Convert a upper triangular matrix to a symmetric one
     */
    inline friend MatType triu2symm(const MatType &a) { return a.zz_triu2symm();}

    /** \brief  Frobenius norm  */
    inline friend MatType norm_F(const MatType &x) { return x.zz_norm_F();}

    /** \brief  2-norm  */
    inline friend MatType norm_2(const MatType &x) { return x.zz_norm_2();}

    /** \brief 1-norm  */
    inline friend MatType norm_1(const MatType &x) { return x.zz_norm_1();}

    /** \brief Infinity-norm */
    inline friend MatType norm_inf(const MatType &x) { return x.zz_norm_inf();}

    /** \brief Return a col-wise summation of elements */
    inline friend MatType sumCols(const MatType &x) { return x.zz_sumCols();}

    /** \brief Return a row-wise summation of elements */
    inline friend MatType sumRows(const MatType &x) { return x.zz_sumRows();}

    /** \brief Inner product of two matrices
        with x and y matrices of the same dimension
    */
    inline friend MatType inner_prod(const MatType &x, const MatType &y) {
      return x.zz_inner_prod(y);
    }

    /** \brief  Take the outer product of two vectors
        Equals
        \code
        x*y.T()
        \endcode
        with x and y vectors
    */
    inline friend MatType outer_prod(const MatType &x, const MatType &y) {
      return x.zz_outer_prod(y);
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
      return A.zz_nullspace();
    }

    /** \brief  Evaluate a polynomial with coefficients p in x
    */
    inline friend MatType polyval(const MatType& p, const MatType& x) {
      return p.zz_polyval(x);
    }

    /** \brief   Get the diagonal of a matrix or construct a diagonal
        When the input is square, the diagonal elements are returned.
        If the input is vector-like, a diagonal matrix is constructed with it. */
    inline friend MatType diag(const MatType &A) {
      return A.zz_diag();
    }

    /** \brief  Unite two matrices no overlapping sparsity
     */
    inline friend MatType unite(const MatType& A, const MatType& B) {
      return A.zz_unite(B);
    }

    /** \brief  Make the matrix dense if not already
     */
    inline friend MatType densify(const MatType& x) {
      return x.zz_densify();
    }

    /** \brief Create a new matrix with a given sparsity pattern but with the
      * nonzeros taken from an existing matrix
      */
    inline friend MatType project(const MatType& A, const Sparsity& sp,
                                  bool intersect=false) {
      return A.zz_project(sp, intersect);
    }

    /** \brief Branching on MX nodes
        Ternary operator, "cond ? if_true : if_false"
    */
    inline friend MatType if_else(const MatType &cond, const MatType &if_true,
                                  const MatType &if_false, bool short_circuit=true) {
      return cond.zz_if_else(if_true, if_false, short_circuit);
    }

    /** \brief Create a switch
     *
     * If the condition \param ind evaluates to the integer k, where 0<=k<f.size(),
     * then x[k] will be returned, otherwise \param x_default will be returned.
     */
    inline friend MatType conditional(const MatType& ind, const std::vector<MatType> &x,
                                      const MatType &x_default, bool short_circuit=true) {
      return ind.zz_conditional(x, x_default, short_circuit);
    }

    /** \brief Check if expression depends on the argument
        The argument must be symbolic
    */
    inline friend bool dependsOn(const MatType& f, const MatType &arg) {
      return f.zz_dependsOn(arg);
    }

    /** \brief  Substitute variable v with expression vdef in an expression ex */
    friend inline MatType substitute(const MatType& ex, const MatType& v,
                                     const MatType& vdef) {
      return ex.zz_substitute(v, vdef);
    }

    /** \brief  Substitute variable var with expression expr in multiple expressions */
    friend inline std::vector<MatType>
      substitute(const std::vector<MatType>& ex, const std::vector<MatType>& v,
                 const std::vector<MatType>& vdef) {
      return MatType::zz_substitute(ex, v, vdef);
    }

    /** \brief Inplace substitution with piggyback expressions
     * Substitute variables v out of the expressions vdef sequentially,
     * as well as out of a number of other expressions piggyback */
    inline friend void
      substituteInPlace(const std::vector<MatType>& v,
                        std::vector<MatType>& inout_vdef,
                        std::vector<MatType>& inout_ex, bool reverse=false) {
      return MatType::zz_substituteInPlace(v, inout_vdef, inout_ex, reverse);
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
      return A.zz_solve(b);
    }

    /** \brief Solve a system of equations: A*x = b
     */
    friend inline MatType solve(const MatType& A, const MatType& b,
                                const std::string& lsolver,
                                const Dict& dict = Dict()) {
      return A.zz_solve(b, lsolver, dict);
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
      return A.zz_pinv();
    }

    /** \brief Computes the Moore-Penrose pseudo-inverse
     *
     * If the matrix A is fat (size1>size2), mul(A, pinv(A)) is unity.
     * If the matrix A is slender (size2<size1), mul(pinv(A), A) is unity.
     *
     */
    friend inline MatType pinv(const MatType& A, const std::string& lsolver,
                               const Dict& dict = Dict()) {
      return A.zz_pinv(lsolver, dict);
    }

    ///@{
    /** \brief Calculate jacobian via source code transformation
    */
    inline friend MatType jacobian(const MatType &ex, const MatType &arg) {
      return ex.zz_jacobian(arg);
    }
    inline friend MatType gradient(const MatType &ex, const MatType &arg) {
      return ex.zz_gradient(arg);
    }
    inline friend MatType tangent(const MatType &ex, const MatType &arg) {
      return ex.zz_tangent(arg);
    }
    ///@}

    ///@{
    // Hessian and (optionally) gradient
    inline friend MatType hessian(const MatType &ex, const MatType &arg) {
      return ex.zz_hessian(arg);
    }
    inline friend MatType hessian(const MatType &ex, const MatType &arg, MatType& output_g) {
      return ex.zz_hessian(arg, output_g);
    }
    ///@}

    /** Count number of nodes */
    inline friend int countNodes(const MatType& A) {
      return A.zz_countNodes();
    }

    /** \brief Get a string representation for a binary MatType, using custom arguments */
    inline friend std::string
      getOperatorRepresentation(const MatType& xb, const std::vector<std::string>& args) {
      return xb.zz_getOperatorRepresentation(args);
    }

    /** \brief Extract shared subexpressions from an set of expressions */
    inline friend void extractShared(std::vector<MatType>& ex,
                                     std::vector<MatType>& v,
                                     std::vector<MatType>& vdef,
                                     const std::string& v_prefix="v_",
                                     const std::string& v_suffix="") {
      MatType::zz_extractShared(ex, v, vdef, v_prefix, v_suffix);
    }

    /** \brief Extract shared subexpressions from an set of expressions */
    inline friend void extractShared(const std::vector<MatType>& ex,
                                     std::vector<MatType>& ex_output,
                                     std::vector<MatType>& v,
                                     std::vector<MatType>& vdef,
                                     const std::string& v_prefix="v_",
                                     const std::string& v_suffix="") {
      ex_output = ex;
      extractShared(ex_output, v, vdef, v_prefix, v_suffix);
    }

    /** \brief Given a repeated matrix, computes the sum of repeated parts
     */
    inline friend MatType repsum(const MatType &A, int n, int m=1) {
      return A.zz_repsum(n, m);
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

#if !defined(SWIG) || !defined(SWIGMATLAB)
    ///@{
    /** \brief [DEPRECATED] Create a sparse matrix with all zeros 
        Use MatType(nrow, ncol) instead **/
    static MatType sparse(int nrow=1, int ncol=1) { return MatType(nrow, ncol);}
    static MatType sparse(const std::pair<int, int>& rc) { return MatType(rc);}
    ///@}

    /** \brief [DEPRECATED] Create a sparse matrix with nonzeros given as a (dense) vector 
        Use MatType(Sparsity, nz) instead **/
    static MatType sparse(const Sparsity& sp, const MatType& nz) { return MatType(sp, nz); }
#endif // !defined(SWIG) || !defined(SWIGMATLAB)

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
#ifdef casadi_implementation
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
  int GenericMatrix<MatType>::shape(int axis) const {
    if (axis==1)
      return sparsity().size1();
    if (axis==2)
      return sparsity().size2();
    casadi_error("Axis must be 1 or 2.");
  }

  template<typename MatType>
  std::string GenericMatrix<MatType>::dimString() const {
    return sparsity().dimString();
  }

  template<typename MatType>
  bool GenericMatrix<MatType>::isscalar(bool scalar_and_dense) const {
    return sparsity().isscalar(scalar_and_dense);
  }

#endif
#endif // SWIG

#ifdef casadi_implementation
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
  MatType GenericMatrix<MatType>::zz_linspace(const MatType &b, int nsteps) const {
    std::vector<MatType> ret(nsteps);
    ret[0] = self();
    MatType step = (b-self())/(nsteps-1);

    for (int i=1; i<nsteps-1; ++i)
      ret[i] = ret[i-1] + step;

    ret[nsteps-1] = b;
    return vertcat(ret);
  }

  template<typename MatType>
  MatType GenericMatrix<MatType>::zz_cross(const MatType &b, int dim) const {
    const MatType &a = self();
    casadi_assert_message(a.size1()==b.size1() && a.size2()==b.size2(),
                          "cross(a, b): Inconsistent dimensions. Dimension of a ("
                          << a.dimString() << " ) must equal that of b ("
                          << b.dimString() << ").");

    casadi_assert_message(a.size1()==3 || a.size2()==3,
                          "cross(a, b): One of the dimensions of a should have length 3, but got "
                          << a.dimString() << ".");
    casadi_assert_message(dim==-1 || dim==1 || dim==2,
                          "cross(a, b, dim): Dim must be 1, 2 or -1 (automatic).");

    std::vector<MatType> ret(3);

    bool t = a.size1()==3;

    if (dim==1) t = true;
    if (dim==2) t = false;

    MatType a1 = t ? a(0, ALL) : a(ALL, 0);
    MatType a2 = t ? a(1, ALL) : a(ALL, 1);
    MatType a3 = t ? a(2, ALL) : a(ALL, 2);

    MatType b1 = t ? b(0, ALL) : b(ALL, 0);
    MatType b2 = t ? b(1, ALL) : b(ALL, 1);
    MatType b3 = t ? b(2, ALL) : b(ALL, 2);

    ret[0] = a2*b3-a3*b2;
    ret[1] = a3*b1-a1*b3;
    ret[2] = a1*b2-a2*b1;

    return t ? vertcat(ret) : horzcat(ret);
  }

  template<typename MatType>
  MatType GenericMatrix<MatType>::zz_tril2symm() const {
    casadi_assert_message(self().issquare(),
                          "Shape error in tril2symm. Expecting square shape but got "
                          << self().dimString());
    casadi_assert_message(self().sizeU()-self().sizeD()==0,
                          "Sparsity error in tril2symm. Found above-diagonal entries in argument: "
                          << self().dimString());
    return self() +  self().T() - diag(diag(self()));
  }

  template<typename MatType>
  MatType GenericMatrix<MatType>::zz_repsum(int n, int m) const {
    casadi_assert(self().size1() % n==0);
    casadi_assert(self().size2() % m==0);
    std::vector< std::vector< MatType> > s =
      blocksplit(self(), self().size1()/n, self().size2()/m);
    MatType sum = 0;
    for (int i=0;i<s.size();++i) {
      for (int j=0;j<s[i].size();++j) {
        sum = sum + s[i][j];
      }
    }
    return sum;
  }

  template<typename MatType>
  MatType GenericMatrix<MatType>::zz_triu2symm() const {
    casadi_assert_message(self().issquare(),
                          "Shape error in triu2symm. Expecting square shape but got "
                          << self().dimString());
    casadi_assert_message(self().sizeL()-self().sizeD()==0,
                          "Sparsity error in triu2symm. Found below-diagonal entries in argument: "
                          << self().dimString());
    return self() + self().T() - diag(diag(self()));
  }

#endif

} // namespace casadi

#endif // CASADI_GENERIC_MATRIX_HPP
