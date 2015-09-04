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
  class CASADI_EXPORT MX :
    public GenericExpression<MX>,
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

    /** \brief Create a sparse matrix from a sparsity pattern.
        Same as MX::ones(sparsity)
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

    /** \brief Get the sparsity pattern */
    const Sparsity& sparsity() const;

    /// Access the sparsity, make a copy if there are multiple references to it
    Sparsity& sparsityRef();

    /// \cond INTERNAL
    /// Scalar type
    typedef MX ScalarType;
    /// \endcond

    /// Base class
    typedef GenericMatrix<MX> B;

    /// Expose base class functions
    using B::zz_horzsplit;
    using B::zz_diagsplit;
    using B::zz_vertsplit;
#endif // SWIG

    /// Returns the truth value of an MX expression
    bool __nonzero__() const;

    /** \brief Get an owning reference to the sparsity pattern */
    Sparsity getSparsity() const { return sparsity();}

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
    int nOut() const;

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

    /** \brief Check if matrix can be used to define function inputs.
        Valid inputs for MXFunctions are combinations of Reshape, concatenations and SymbolicMX
    */
    bool isValidInput() const;

    /** \brief Get the number of symbolic primitive
        Assumes isValidInput() returns true.
    */
    int numPrimitives() const;

    /** \brief Get symbolic primitives */
    std::vector<MX> getPrimitives() const;

    /** \brief Split up an expression along symbolic primitives */
    std::vector<MX> splitPrimitives(const MX& x) const;

    /** \brief Join an expression along symbolic primitives */
    MX joinPrimitives(std::vector<MX>& v) const;

    /// \cond INTERNAL
    /** \brief Detect duplicate symbolic expressions
        If there are symbolic primitives appearing more than once, the function will return
        true and the names of the duplicate expressions will be printed to userOut<true, PL_WARN>().
        Note: Will mark the node using MX::setTemp.
        Make sure to call resetInput() after usage.
    */
    bool hasDuplicates();

    /** \brief Reset the marker for an input expression */
    void resetInput();
  /// \endcond

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

    /** \brief  Number of functions */
    int numFunctions() const;

    /// Get function
    Function getFunction(int i=0);

    /// Is binary operation
    bool isBinary() const;

    /// Is unary operation
    bool isUnary() const;

    /// Get operation type
    int getOp() const;

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

#ifndef SWIG
    /** \brief Avoid shadowing SharedObject::get() */
    using SharedObject::get;
#endif // SWIG

    ///@{
    /// Get a submatrix, single argument
    void get(MX& SWIG_OUTPUT(m), bool ind1, const Slice& rr) const;
    void get(MX& SWIG_OUTPUT(m), bool ind1, const Matrix<int>& rr) const;
    void get(MX& SWIG_OUTPUT(m), bool ind1, const Sparsity& sp) const;
    ///@}

    /// Get a submatrix, two arguments
    ///@{
    void get(MX& SWIG_OUTPUT(m), bool ind1, const Slice& rr, const Slice& cc) const;
    void get(MX& SWIG_OUTPUT(m), bool ind1, const Slice& rr, const Matrix<int>& cc) const;
    void get(MX& SWIG_OUTPUT(m), bool ind1, const Matrix<int>& rr, const Slice& cc) const;
    void get(MX& SWIG_OUTPUT(m), bool ind1, const Matrix<int>& rr, const Matrix<int>& cc) const;
    ///@}

    ///@{
    /// Set a submatrix, single argument
    void set(const MX& m, bool ind1, const Slice& rr);
    void set(const MX& m, bool ind1, const Matrix<int>& rr);
    void set(const MX& m, bool ind1, const Sparsity& sp);
    ///@}

    ///@{
    /// Set a submatrix, two arguments
    ///@}
    void set(const MX& m, bool ind1, const Slice& rr, const Slice& cc);
    void set(const MX& m, bool ind1, const Slice& rr, const Matrix<int>& cc);
    void set(const MX& m, bool ind1, const Matrix<int>& rr, const Slice& cc);
    void set(const MX& m, bool ind1, const Matrix<int>& rr, const Matrix<int>& cc);
    ///@}

    ///@{
    /// Get a set of nonzeros
    void getNZ(MX& SWIG_OUTPUT(m), bool ind1, const Slice& kk) const;
    void getNZ(MX& SWIG_OUTPUT(m), bool ind1, const Matrix<int>& kk) const;
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

#ifndef SWIG
    /// \cond CLUTTER
    ///@{
    /// Functions called by friend functions defined for GenericExpression
    MX zz_plus(const MX& y) const;
    MX zz_minus(const MX& y) const;
    MX zz_times(const MX& y) const;
    MX zz_rdivide(const MX& y) const;
    MX zz_lt(const MX& y) const;
    MX zz_le(const MX& y) const;
    MX zz_eq(const MX& y) const;
    MX zz_ne(const MX& y) const;
    MX zz_atan2(const MX& y) const;
    MX zz_min(const MX& y) const;
    MX zz_max(const MX& y) const;
    MX zz_and(const MX& y) const;
    MX zz_or(const MX& y) const;
    MX zz_abs() const;
    MX zz_sqrt() const;
    MX zz_sin() const;
    MX zz_cos() const;
    MX zz_tan() const;
    MX zz_asin() const;
    MX zz_acos() const;
    MX zz_atan() const;
    MX zz_sinh() const;
    MX zz_cosh() const;
    MX zz_tanh() const;
    MX zz_asinh() const;
    MX zz_acosh() const;
    MX zz_atanh() const;
    MX zz_exp() const;
    MX zz_log() const;
    MX zz_log10() const;
    MX zz_floor() const;
    MX zz_ceil() const;
    MX zz_erf() const;
    MX zz_erfinv() const;
    MX zz_sign() const;
    MX zz_power(const MX& b) const;
    MX zz_mod(const MX& y) const;
    MX zz_simplify() const;
    bool zz_isEqual(const MX& y, int depth) const;
    bool zz_isEqual(const MXNode* y, int depth) const;
    MX zz_copysign(const MX& y) const;
    MX zz_constpow(const MX& y) const;
    ///@}

    ///@{
    /// Functions called by friend functions defined for GenericExpression
    MX zz_jacobian(const MX& arg) const;
    MX zz_gradient(const MX& arg) const;
    MX zz_tangent(const MX& arg) const;
    MX zz_hessian(const MX& arg) const;
    MX zz_hessian(const MX& arg, MX& g) const;
    MX zz_substitute(const MX& v, const MX& vdef) const;
    static std::vector<MX> zz_substitute(const std::vector<MX> &ex,
                                         const std::vector<MX> &v,
                                         const std::vector<MX> &vdef);
    static void zz_substituteInPlace(const std::vector<MX>& v,
                                     std::vector<MX>& vdef,
                                     std::vector<MX>& ex, bool reverse);
    MX zz_solve(const MX& b, const std::string& lsolver="symbolicqr",
                const Dict& dict = Dict()) const;
    MX zz_pinv(const std::string& lsolver="symbolicqr",
               const Dict& dict = Dict()) const;
    int zz_countNodes() const;
    std::string zz_getOperatorRepresentation(const std::vector<std::string>& args) const;
    static void zz_extractShared(std::vector<MX>& ex, std::vector<MX>& v,
                                 std::vector<MX>& vdef, const std::string& v_prefix,
                                 const std::string& v_suffix);
    MX zz_if_else(const MX& if_true, const MX& if_false, bool short_circuit=true) const;
    MX zz_conditional(const std::vector<MX> &x, const MX& x_default,
                      bool short_circuit=true) const;
    bool zz_dependsOn(const MX& arg) const;
    MX zz_not() const;
    ///@}

    ///@{
    /// Functions called by friend functions defined for SparsityInterface
    static MX zz_horzcat(const std::vector<MX>& x);
    static MX zz_diagcat(const std::vector<MX>& x);
    static MX zz_vertcat(const std::vector<MX>& x);
    std::vector<MX> zz_horzsplit(const std::vector<int>& offset) const;
    std::vector<MX> zz_diagsplit(const std::vector<int>& offset1,
                                 const std::vector<int>& offset2) const;
    std::vector<MX> zz_vertsplit(const std::vector<int>& offset) const;
    static MX zz_blockcat(const std::vector< std::vector<MX > > &v);
    MX zz_mtimes(const MX& y) const;
    MX zz_mac(const MX& y, const MX& z) const;
    MX zz_reshape(int nrow, int ncol) const;
    MX zz_reshape(const Sparsity& sp) const;
    MX zz_vecNZ() const;
    MX zz_kron(const MX& b) const;
    MX zz_repmat(int n, int m=1) const;
    ///@}

    ///@{
    /// Functions called by friend functions defined for GenericMatrix
    MX zz_mpower(const MX& b) const;
    MX zz_inner_prod(const MX& y) const;
    MX zz_outer_prod(const MX& y) const;
    MX zz_mrdivide(const MX& b) const;
    MX zz_mldivide(const MX& b) const;
    MX zz_if_else_zero(const MX& y) const;
    MX zz_norm_2() const;
    MX zz_norm_F() const;
    MX zz_norm_1() const;
    MX zz_norm_inf() const;
    MX zz_unite(const MX& B) const;
    MX zz_trace() const;
    MX zz_diag() const;
    MX zz_sumCols() const;
    MX zz_sumRows() const;
    MX zz_polyval(const MX& x) const;
    MX zz_det() const;
    MX zz_inv() const;
    std::vector<MX> zz_symvar() const;
    MX zz_nullspace() const;
    MX zz_repsum(int n, int m=1) const;
    ///@}

    ///@{
    /// Functions called by friend functions defined for this class
    MX zz_find() const;

    MX zz_graph_substitute(const std::vector<MX> &v, const std::vector<MX> &vdef) const;
    static std::vector<MX> zz_graph_substitute(const std::vector<MX> &ex,
                                               const std::vector<MX> &expr,
                                               const std::vector<MX> &exprs);
    static MX zz_matrix_expand(const MX& e, const std::vector<MX> &boundary);
    static std::vector<MX> zz_matrix_expand(const std::vector<MX>& e,
                                            const std::vector<MX>& boundary);
    ///@}
    /// \endcond

#endif // SWIG

    MX printme(const MX& y) const;

#if !defined(SWIG) || defined(DOXYGEN)
/**
\ingroup expression_tools
@{
*/
    /** \brief Find first nonzero
     * If failed, returns the number of rows
     */
    inline friend MX find(const MX& x) {
      return x.zz_find();
    }

    /** \brief Substitute single expression in graph
     * Substitute variable v with expression vdef in an expression ex, preserving nodes
     */
    inline friend MX graph_substitute(const MX& ex, const std::vector<MX> &v,
                                      const std::vector<MX> &vdef) {
      return ex.zz_graph_substitute(v, vdef);
    }

    /** \brief Substitute multiple expressions in graph
     * Substitute variable var with expression expr in
     * multiple expressions, preserving nodes
     */
    inline friend std::vector<MX>
      graph_substitute(const std::vector<MX> &ex,
                       const std::vector<MX> &v,
                       const std::vector<MX> &vdef) {
      return MX::zz_graph_substitute(ex, v, vdef);
    }

    /** \brief Expand MX graph to SXFunction call
     *
     *  Expand the given expression e, optionally
     *  supplying expressions contained in it at which expansion should stop.
     *
     */
    inline friend MX
      matrix_expand(const MX& e, const std::vector<MX> &boundary = std::vector<MX>()) {
      return MX::zz_matrix_expand(e, boundary);
    }

    /** \brief Expand MX graph to SXFunction call
     *
     *  Expand the given expression e, optionally
     *  supplying expressions contained in it at which expansion should stop.
     *
     */
    inline friend std::vector<MX>
      matrix_expand(const std::vector<MX>& e,
                    const std::vector<MX> &boundary = std::vector<MX>()) {
      return MX::zz_matrix_expand(e, boundary);
    }
/** @} */
#endif // SWIG

    /** \brief returns itself, but with an assertion attached
    *
    *  If y does not evaluate to 1, a runtime error is raised
    */
    MX attachAssert(const MX& y, const std::string& fail_message="") const;

    /** \brief Monitor an expression
    * Returns itself, but with the side effect of printing the nonzeros along with a comment
    */
    MX monitor(const std::string& comment) const;

#if !defined(SWIG) || !defined(SWIGMATLAB)

    /** \brief Set sparse */
    MX zz_project(const Sparsity& sp, bool intersect=false) const;

#endif // !defined(SWIG) || !defined(SWIGMATLAB)

    /// Make the matrix dense
    void makeDense(const MX& val = 0);

    /// Lift an expression
    void lift(const MX& x_guess);

    /// Transpose the matrix
    MX T() const;

    /** \brief Get an IMatrix representation of a GetNonzeros or SetNonzeros node */
    Matrix<int> mapping() const;

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

    // Depth when checking equalities
    static int eq_depth_;

#endif // SWIG
  };

  // Create matrix symbolic primitive
  template<>
  MX GenericMatrix<MX>::sym(const std::string& name, const Sparsity& sp);

  ///@{
  /// Readability typedefs
  typedef std::vector<MX> MXVector;
  typedef std::vector<MXVector> MXVectorVector;
  typedef std::map<std::string, MX> MXDict;
  ///@}

} // namespace casadi

#endif // CASADI_MX_HPP
