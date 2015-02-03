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


#ifndef CASADI_MX_HPP
#define CASADI_MX_HPP
#include "../shared_object.hpp"
#include "../matrix/matrix.hpp"
#include "../matrix/generic_expression.hpp"
#include "../generic_type.hpp"
#include <vector>
namespace casadi {

  /** \brief  Forward declaration */
  class MXNode;
  class Function;


  /** \brief MX - Matrix expression

      The MX class is used to build up trees made up from MXNodes. It is a more general
      graph representation than the scalar expression, SX, and much less efficient for small
      objects. On the other hand, the class allows much more general operations than does SX,
      in particular matrix valued operations and calls to arbitrary differentiable functions.

      The MX class is designed to have identical syntax with the Matrix<> template class,
      and uses Matrix<double> as its internal representation of the values at a node. By keeping
      the syntaxes identical, it is possible to switch from one class to the other,
      as well as inlining MX functions to SXElement functions.

      Note that an operation is always "lazy", making a matrix multiplication will create a
      matrix multiplication node, not perform the actual multiplication.

      \author Joel Andersson
      \date 2010-2011
  */
  class CASADI_EXPORT MX : public GenericExpression<MX>,
                           public GenericMatrix<MX>,
                           public SharedObject {
  public:

    /** \brief  Default constructor */
    MX();

    /** \brief Create a sparse matrix with all structural zeros */
    MX(int nrow, int ncol);

#ifndef SWIG
    /** \brief Create a sparse matrix with all structural zeros */
    explicit MX(const std::pair<int, int>& rc);
#endif // SWIG

    /** \brief Sparse matrix with a given sparsity and zero entries
        Same as MX::zeros(sparsity)
     */
    explicit MX(const Sparsity& sp);

    /** \brief Construct matrix with a given sparsity and nonzeros */
    MX(const Sparsity& sp, const MX& val);

    /** \brief  Create scalar constant (also implicit type conversion) */
    MX(double x);

    /** \brief  Copy constructor */
    MX(const MX& x);

    /** \brief  Create vector constant (also implicit type conversion) */
    MX(const std::vector<double> &x);

    /** \brief  Create sparse matrix constant (also implicit type conversion) */
    MX(const Matrix<double> &x);

/// \cond INTERNAL
    /** \brief  Destructor */
    ~MX();
/// \endcond

#ifndef SWIG
/// \cond INTERNAL
    /** \brief  Create from node */
    static MX create(MXNode* node);

    /** \brief  Create from node (multiple-outputs) */
    static std::vector<MX> createMultipleOutput(MXNode* node);
/// \endcond

    /// Get a non-zero element, with bounds checking
    const MX at(int k) const;

    /// Access a non-zero element, with bounds checking
    NonZeros<MX, int> at(int k);

#endif // SWIG

    /// Returns the truth value of an MX expression
    bool __nonzero__() const;

    /// \cond INTERNAL
    /// Scalar type
    typedef MX ScalarType;
    /// \endcond

    /** \brief Get the sparsity pattern */
    const Sparsity& sparsity() const;

    /// Access the sparsity, make a copy if there are multiple references to it
    Sparsity& sparsityRef();

    /** \brief Erase a submatrix (leaving structural zeros in its place)
        Erase rows and/or columns of a matrix */
    void erase(const std::vector<int>& rr, const std::vector<int>& cc, bool ind1=false);

    /** \brief Erase a submatrix (leaving structural zeros in its place)
        Erase elements of a matrix */
    void erase(const std::vector<int>& rr, bool ind1=false);

    /** \brief Enlarge matrix
        Make the matrix larger by inserting empty rows and columns, keeping the existing non-zeros */
    void enlarge(int nrow, int ncol,
                 const std::vector<int>& rr, const std::vector<int>& cc, bool ind1=false);

    MX operator-() const;

    /// Get the non-zero elements
    MX nonzeros() const;

#ifndef SWIG
    /// \cond INTERNAL
    ///@{
    /** \brief  Access a member of the node */
    MXNode* operator->();

    /** \brief  Const access a member of the node */
    const MXNode* operator->() const;
    ///@}
    /// \endcond
#endif // SWIG

    /** \brief Get the nth dependency as MX */
    MX getDep(int ch=0) const;

    /** \brief  Number of outputs */
    int getNumOutputs() const;

    /** \brief  Get an output */
    MX getOutput(int oind=0) const;

    /** \brief Get the number of dependencies of a binary SXElement */
    int getNdeps() const;

    /// Get the name.
    std::string getName() const;

    /// Get the value (only for scalar constant nodes)
    double getValue() const;

    /// Get the value (only for constant nodes)
    Matrix<double> getMatrixValue() const;

    /// Check if symbolic
    bool isSymbolic() const;

    /// Check if constant
    bool isConstant() const;

    /// Check if evaluation
    bool isEvaluation() const;

    /// Check if evaluation output
    bool isEvaluationOutput() const;

    /// Get the index of evaluation output - only valid when isEvaluationoutput() is true
    int getEvaluationOutput() const;

    /// Is it a certain operation
    bool isOperation(int op) const;

    /// Check if multiplication
    bool isMultiplication() const;

    /// Check if commutative operation
    bool isCommutative() const;

    /// Check if norm
    bool isNorm() const;

    /** \brief  check if all nonzeros are symbolic
     * (this function is currently identical to isSymbolic) */
    bool isSymbolicSparse() const;

    /** \brief  check if identity */
    bool isIdentity() const;

    /** \brief  check if zero (note that false negative answers are possible) */
    bool isZero() const;

    /** \brief  check if zero (note that false negative answers are possible) */
    bool isOne() const;

    /** \brief  check if zero (note that false negative answers are possible) */
    bool isMinusOne() const;

    /** \brief  Is the expression a transpose? */
    bool isTranspose() const;

    /// Checks if expression does not contain NaN or Inf
    bool isRegular() const;

    /// Get function
    Function getFunction();

    /// Is binary operation
    bool isBinary() const;

    /// Is unary operation
    bool isUnary() const;

    /// Get operation type
    int getOp() const;

    /** \brief Returns a number that is unique for a given MXNode.
     * If the MX does not point to any node, 0 is returned.
     */
    long __hash__() const;

    /// \cond INTERNAL
    /// Get the temporary variable
    int getTemp() const;

    /// Set the temporary variable
    void setTemp(int t);
    /// \endcond

    ///@{
    /** \brief  Create nodes by their ID */
    static MX binary(int op, const MX &x, const MX &y);
    static MX unary(int op, const MX &x);
    ///@}

    ///@{
    /** \brief  create a matrix with all inf */
    static MX inf(const Sparsity& sp);
    static MX inf(int nrow=1, int ncol=1);
    static MX inf(const std::pair<int, int>& rc);
    ///@}

    ///@{
    /** \brief  create a matrix with all nan */
    static MX nan(const Sparsity& sp);
    static MX nan(int nrow=1, int ncol=1);
    static MX nan(const std::pair<int, int>& rc);
    ///@}

    /** \brief  Identity matrix */
    static MX eye(int ncol);

    ///@{
    /// Get a submatrix, single argument
    const MX getSub(bool ind1, const Slice& rr) const;
    const MX getSub(bool ind1, const Matrix<int>& rr) const;
    const MX getSub(bool ind1, const Sparsity& sp) const;
    ///@}

    /// Get a submatrix, two arguments
    ///@{
    const MX getSub(bool ind1, const Slice& rr, const Slice& cc) const;
    const MX getSub(bool ind1, const Slice& rr, const Matrix<int>& cc) const;
    const MX getSub(bool ind1, const Matrix<int>& rr, const Slice& cc) const;
    const MX getSub(bool ind1, const Matrix<int>& rr, const Matrix<int>& cc) const;
    ///@}

    ///@{
    /// Set a submatrix, single argument
    void setSub(const MX& m, bool ind1, const Slice& rr);
    void setSub(const MX& m, bool ind1, const Matrix<int>& rr);
    void setSub(const MX& m, bool ind1, const Sparsity& sp);
    ///@}

    ///@{
    /// Set a submatrix, two arguments
    ///@}
    void setSub(const MX& m, bool ind1, const Slice& rr, const Slice& cc);
    void setSub(const MX& m, bool ind1, const Slice& rr, const Matrix<int>& cc);
    void setSub(const MX& m, bool ind1, const Matrix<int>& rr, const Slice& cc);
    void setSub(const MX& m, bool ind1, const Matrix<int>& rr, const Matrix<int>& cc);
    ///@}

    ///@{
    /// Get a set of nonzeros
    MX getNZ(bool ind1, const Slice& kk) const;
    MX getNZ(bool ind1, const Matrix<int>& kk) const;
    ///@}

    ///@{
    /// Set a set of nonzeros
    void setNZ(const MX& m, bool ind1, const Slice& kk);
    void setNZ(const MX& m, bool ind1, const Matrix<int>& kk);
    ///@}

    /** \brief Append a matrix vertically (NOTE: only efficient if vector) */
    void append(const MX& y);

    /** \brief Append a matrix horizontally */
    void appendColumns(const MX& y);

    /// \cond CLUTTER
    /// all binary operations
    MX zz_plus(const MX& y) const;
    MX zz_minus(const MX& y) const;
    MX zz_times(const MX& y) const;
    MX zz_rdivide(const MX& y) const;
    MX zz_lt(const MX& y) const;
    MX zz_le(const MX& y) const;
    MX zz_eq(const MX& y) const;
    MX zz_ne(const MX& y) const;
    MX zz_power(const MX& b) const;
    MX zz_mpower(const MX& b) const;
    MX zz_inner_prod(const MX& y) const;
    MX zz_outer_prod(const MX& y) const;
    MX zz_min(const MX& y) const;
    MX zz_max(const MX& y) const;
    MX zz_mod(const MX& y) const;
    MX zz_atan2(const MX& y) const;
    MX zz_and(const MX& y) const;
    MX zz_or(const MX& y) const;
    MX zz_if_else_zero(const MX& y) const;
    /// \endcond

    /// \cond CLUTTER
    MX __truediv__(const MX& y) const { return zz_rdivide(y);}
    MX __constpow__(const MX& b) const;
    MX __mrdivide__(const MX& b) const;
    MX __copysign__(const MX& y) const;
    MX constpow(const MX& y) const;
    /// \endcond
    MX printme(const MX& y) const;

    /// \cond CLUTTER

    // all unary operations
    MX zz_exp() const;
    MX zz_log() const;
    MX zz_log10() const;
    MX zz_sqrt() const;
    MX zz_sin() const;
    MX zz_cos() const;
    MX zz_tan() const;
    MX zz_asin() const;
    MX zz_acos() const;
    MX zz_atan() const;
    MX zz_floor() const;
    MX zz_ceil() const;
    MX zz_abs() const;
    MX zz_sign() const;
    MX zz_erfinv() const;
    MX zz_erf() const;
    MX zz_sinh() const;
    MX zz_cosh() const;
    MX zz_tanh() const;
    MX zz_asinh() const;
    MX zz_acosh() const;
    MX zz_atanh() const;
    MX zz_not() const;
    static MX zz_horzcat(const std::vector<MX>& x);
    static MX zz_diagcat(const std::vector<MX>& x);
    static MX zz_vertcat(const std::vector<MX>& x);
    std::vector<MX> zz_horzsplit(const std::vector<int>& offset) const;
    std::vector<MX> zz_diagsplit(const std::vector<int>& offset1,
                                 const std::vector<int>& offset2) const;
    std::vector<MX> zz_vertsplit(const std::vector<int>& offset) const;
    static MX zz_blockcat(const std::vector< std::vector<MX > > &v);
    MX zz_norm_2() const;
    MX zz_norm_F() const;
    MX zz_norm_1() const;
    MX zz_norm_inf() const;
    MX zz_mtimes(const MX& y) const;
    MX zz_mtimes(const MX& y, const MX& z) const;
    MX zz_simplify() const;
    MX zz_reshape(int nrow, int ncol) const;
    MX zz_reshape(const Sparsity& sp) const;
    MX zz_vecNZ() const;
    MX zz_if_else(const MX &if_true, const MX &if_false) const;
    MX zz_unite(const MX& B) const;
    MX zz_trace() const;
    MX zz_repmat(const Sparsity& sp) const;
    MX zz_repmat(int n, int m=1) const;
    static MX zz_createParent(std::vector<MX> &deps);
    static MX zz_createParent(const std::vector<Sparsity> &deps, std::vector<MX>& children);
    static MX zz_createParent(const std::vector<MX> &deps, std::vector<MX>& children);
    MX zz_diag() const;
    int zz_countNodes() const;
    MX zz_sumCols() const;
    MX zz_sumRows() const;
    MX zz_sumAll() const;
    MX zz_polyval(const MX& x) const;
    std::string zz_getOperatorRepresentation(const std::vector<std::string>& args) const;
    static void zz_substituteInPlace(const std::vector<MX>& v,
                                     std::vector<MX>& vdef, bool reverse=false);
    static void zz_substituteInPlace(const std::vector<MX>& v, std::vector<MX>& vdef,
                              std::vector<MX>& ex, bool reverse=false);
    MX zz_substitute(const MX& v, const MX& vdef) const;
    static std::vector<MX> zz_substitute(const std::vector<MX> &ex,
                                         const std::vector<MX> &v,
                                         const std::vector<MX> &vdef);
    MX zz_graph_substitute(const std::vector<MX> &v, const std::vector<MX> &vdef) const;
    static std::vector<MX> zz_graph_substitute(const std::vector<MX> &ex,
                                               const std::vector<MX> &expr,
                                               const std::vector<MX> &exprs);
    static void zz_extractShared(std::vector<MX>& ex, std::vector<MX>& v, std::vector<MX>& vdef,
                                 const std::string& v_prefix="v_", const std::string& v_suffix="");
    void zz_printCompact(std::ostream &stream=CASADI_COUT) const;
    MX zz_jacobian(const MX &arg) const;
    MX zz_gradient(const MX &arg) const;
    MX zz_tangent(const MX &arg) const;
    MX zz_det() const;
    MX zz_inv() const;
    std::vector<MX> zz_getSymbols() const;
    static std::vector<MX> zz_getSymbols(const std::vector<MX>& e);
    bool zz_dependsOn(const std::vector<MX> &arg) const;
    static MX zz_matrix_expand(const MX& e,
                               const std::vector<MX> &boundary = std::vector<MX>());
    static std::vector<MX> zz_matrix_expand(const std::vector<MX>& e,
                                            const std::vector<MX>& boundary = std::vector<MX>());
    MX zz_kron(const MX& b) const;
    MX zz_solve(const MX& b, const std::string& lsolver,
                const Dictionary& dict = Dictionary()) const;
    MX zz_pinv(const std::string& lsolver, const Dictionary& dict = Dictionary()) const;
    MX zz_nullspace() const;
    bool zz_isEqual(const MX& y, int depth=0) const;
#ifndef SWIG
    bool zz_isEqual(const MXNode* y, int depth=0) const;
#endif // SWIG

    /// \endcond

    /** \brief returns itself, but with an assertion attached
    *
    *  If y does not evaluate to 1, a runtime error is raised
    */
    MX attachAssert(const MX& y, const std::string &fail_message="") const;

    /** \brief Set sparse */
    MX setSparse(const Sparsity& sp, bool intersect=false) const;

    /// Make the matrix dense
    void makeDense(const MX& val = 0);

    /// Lift an expression
    void lift(const MX& x_guess);

    /// Add an expression to the expression if the expression is non-empty, otherwise assign
    void addToSum(const MX& x);

    /// Transpose the matrix
    MX T() const;

    /** \brief Get an IMatrix representation of a GetNonzeros or SetNonzeros node */
    Matrix<int> mapping() const;

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

    /// Check if a particular cast is allowed
    static bool testCast(const SharedObjectNode* ptr);

#ifndef SWIG
    /// Construct constant matrix with a given sparsity and all
    MX(const Sparsity& sp, int val, bool dummy);
    MX(const Sparsity& sp, double val, bool dummy);
  private:

    /// Create an expression from a node: extra dummy arguments to avoid ambiguity for 0/NULL
    MX(MXNode* node, bool dummy1, bool dummy2, bool dummy3, bool dummy4);

    // Maximum number of calls
    static long max_num_calls_in_print_;

    // Depth when checking equalities
    static int eq_depth_;

#endif // SWIG
  };

  // Create matrix symbolic primitive
  template<>
  MX GenericMatrix<MX>::sym(const std::string& name, const Sparsity& sp);

  ///@{
  /// Some typedefs
  typedef std::vector<MX> MXVector;
  typedef std::vector< std::vector<MX> > MXVectorVector;
  /// \cond INTERNAL
  typedef MX* MXPtr;
  typedef std::vector<MXPtr> MXPtrV;
  typedef std::vector<MXPtrV> MXPtrVV;
  /// \endcond
  ///@}

} // namespace casadi

#endif // CASADI_MX_HPP
