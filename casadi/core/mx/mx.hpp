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
  class CASADI_CORE_EXPORT MX : public GenericExpression<MX>,
                                public GenericMatrix<MX>,
                                public SharedObject {
  public:

    /** \brief  Default constructor */
    MX();

    ///@{
    /** \brief Construct constant matrix with a given sparsity */
    MX(const Sparsity& sp, int val=0);
    MX(const Sparsity& sp, double val);
    MX(const Sparsity& sp, const MX& val);
    ///@}

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

    /// \cond CLUTTER
    ///@{
    /// Indexing for interfaced languages

    /// get a non-zero
    const MX nz_indexed_one_based(int k) const { return at(k-1);}
    const MX nz_indexed_zero_based(int k) const { return at(k);}
    const MX nz_indexed(const IndexList &k) const {
      return (*this)[k.getAll(size())];
    }
    const MX nz_indexed(const Slice &k) const {
      return (*this)[k.getAll(size())];
    }
    const MX nz_indexed(const Matrix<int> &k) const {
      return (*this)[k];
    }

    /// get a matrix element
    const MX indexed_one_based(int rr, int cc) const { return (*this)(rr-1, cc-1);}
    const MX indexed_zero_based(int rr, int cc) const { return (*this)(rr, cc);}
    const MX indexed(const IndexList &rr, const IndexList &cc) const {
      return (*this)(rr.getAll(size1()), cc.getAll(size2()));
    }
    const MX indexed(const Slice &rr, const Slice &cc) const {
      return (*this)(rr.getAll(size1()), cc.getAll(size2()));
    }
    const MX indexed(const Matrix<int> &k) const {
      return (*this)(k);
    }
    const MX indexed(const Sparsity &sp) const {
      return (*this)(sp);
    }
    const MX indexed(const Slice &rr, const Matrix<int>& cc) const { return (*this)(rr, cc); }
    const MX indexed(const Matrix<int>& rr, const IndexList& cc) const {
      return (*this)(rr, cc.getAll(size2()));
    }
    const MX indexed(const Matrix<int>& rr, const Slice& cc) const { return (*this)(rr, cc); }
    const MX indexed(const IndexList& rr, const Matrix<int>& cc) const {
      return (*this)(rr.getAll(size1()), cc);
    }
    const MX indexed(const Matrix<int>& rr, const Matrix<int>& cc) const {
      return (*this)(rr, cc);
    }

    // get a vector element
    const MX indexed_one_based(int rr) const {
      casadi_assert_message(isDense() && isVector(),
                            "Matrix must be a dense vector, but got " << dimString() << ".");
      return (*this)(rr-1);
    }
    const MX indexed_zero_based(int rr) const {
      casadi_assert_message(isDense() && isVector(),
                            "Matrix must be a dense vector, but got " << dimString() << ".");
      return (*this)(rr);
    }
    const MX indexed(const IndexList &rr) const {
      casadi_assert_message(isDense() && isVector(), "Matrix must be a dense vector, but got "
                            << dimString() << ".");
      return (*this)(rr.getAll(size1()));
    }
    const MX indexed(const Slice &rr) const {
      casadi_assert_message(isDense() && isVector(), "Matrix must be a dense vector, but got "
                            << dimString() << ".");
      return (*this)(rr.getAll(size1()));
    }

    /// set a non-zero
    void nz_indexed_one_based_assignment(int k, const MX &m) { at(k-1) = m(0, 0);}
    void nz_indexed_zero_based_assignment(int k, const MX &m) { at(k) = m(0, 0);}
    void nz_indexed_assignment(const IndexList &k, const MX &m) {
      (*this)[k.getAll(size())] = m;
    }
    void nz_indexed_assignment(const Slice &k, const MX &m) {
      (*this)[k.getAll(size())] = m;
    }
    void nz_indexed_assignment(const Matrix<int> &k, const MX &m) {
      (*this)[k] = m;
    }

    /// set a matrix element
    void indexed_one_based_assignment(int rr, int cc, const MX &m) { (*this)(rr-1, cc-1) = m;}
    void indexed_zero_based_assignment(int rr, int cc, const MX &m) { (*this)(rr, cc) = m;}
    void indexed_assignment(const IndexList &rr, const IndexList &cc, const MX &m) {
      setSub(m, rr.getAll(size1()), cc.getAll(size2()));
    }

    void indexed_assignment(const Slice &rr, const Slice &cc, const MX &m) {
      (*this)(rr.getAll(size1()), cc.getAll(size2())) = m;
    }

    void indexed_zero_based_assignment(const Matrix<int>& k, const MX &m) {
      (*this)[k] = m;
    }
    void indexed_assignment(const Sparsity& sp, const MX &m) {
      (*this)(sp) = m;
    }
    void indexed_assignment(const Matrix<int> &rr, const Slice& cc, const MX& m) {
      (*this)(rr, cc.getAll(size2())) = m;
    }
    void indexed_assignment(const Slice& rr, const Matrix<int>& cc, const MX& m) {
      (*this)(rr.getAll(size1()), cc) = m;
    }
    void indexed_assignment(const Matrix<int>& rr, const IndexList& cc, const MX& m) {
      (*this)(rr, cc.getAll(size2())) = m;
    }
    void indexed_assignment(const IndexList& rr, const Matrix<int>& cc, const MX& m) {
      (*this)(rr.getAll(size1()), cc) = m;
    }
    void indexed_assignment(const Matrix<int>& rr, const Matrix<int>& cc, const MX& m) {
      (*this)(rr, cc) = m;
    }
    ///@}

    // set a vector element
    void indexed_one_based_assignment(int rr, const MX &m) {
      casadi_assert_message(isDense() && isVector(),
                            "Matrix must be a dense vector, but got " << dimString() << ".");
      (*this)(rr-1) = m;
    }
    void indexed_zero_based_assignment(int rr, const MX &m) {
      casadi_assert_message(isDense() && isVector(),
                            "Matrix must be a dense vector, but got " << dimString() << ".");
      (*this)(rr) = m;
    }
    void indexed_assignment(const IndexList &rr, const MX &m) {
      casadi_assert_message(isDense() && isVector(),
                            "Matrix must be a dense vector, but got " << dimString() << ".");
      (*this)(rr.getAll(size1())) = m;
    }

    void indexed_assignment(const Slice &rr, const MX &m) {
      (*this)(rr.getAll(size1())) = m;
    }

    /// \endcond

    /// \cond INTERNAL
    /// Scalar type
    typedef MX ScalarType;
    /// \endcond

    /** \brief Get the sparsity pattern */
    const Sparsity& sparsity() const;

    /// Access the sparsity, make a copy if there are multiple references to it
    Sparsity& sparsityRef();

    /** \brief Erase a submatrix */
    void erase(const std::vector<int>& rr, const std::vector<int>& cc);

    /** \brief Enlarge matrix
        Make the matrix larger by inserting empty rows and columns, keeping the existing non-zeros */
    void enlarge(int nrow, int ncol, const std::vector<int>& rr, const std::vector<int>& cc);

    MX operator-() const;

    /// \cond INTERNAL
    ///@{
    /** \brief  Access a member of the node */
    MXNode* operator->();

    /** \brief  Const access a member of the node */
    const MXNode* operator->() const;
    ///@}
    /// \endcond

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

    /** \brief Check if two nodes are equivalent up to a given depth.
     *  Depth=0 checks if the expressions are identical, i.e. points to the same node.
     *
     *  a = x*x
     *  b = x*x
     *
     *  a.isEqual(b, 0)  will return false, but a.isEqual(b, 1) will return true
     */
    bool isEqual(const MX& y, int depth=0) const;
#ifndef SWIG
    bool isEqual(const MXNode* y, int depth=0) const;
#endif // SWIG

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

    ///@{
    /** \brief  create a matrix by repeating an existing matrix */
    static MX repmat(const MX& x, const Sparsity& sp);
    static MX repmat(const MX& x, int nrow, int ncol=1);
    static MX repmat(const MX& x, const std::pair<int, int> &rc);
    ///@}

    /** \brief  Identity matrix */
    static MX eye(int ncol);


    /// \cond INTERNAL
    const MX sub(int rr, int cc) const;
    const MX sub(const std::vector<int>& rr, int cc) const;
    const MX sub(int rr, const std::vector<int>& cc) const;
    const MX sub(const std::vector<int>& rr, const std::vector<int>& cc) const;
    const MX sub(const Sparsity& sp, int dummy=0) const;
    const MX sub(const Matrix<int>& rr, int cc) const { return sub(rr, Slice(cc));}
    const MX sub(const Matrix<int>& rr, const std::vector<int>& cc) const;
    const MX sub(const std::vector<int>& rr, const Matrix<int>& cc) const;
    const MX sub(int rr, const Slice& cc) const {return sub(rr, cc.getAll(size2()));}
    const MX sub(const Slice& rr, int cc) const {return sub(rr.getAll(size1()), cc);}
    const MX sub(const Slice& rr, const Slice& cc) const
    { return sub(rr.getAll(size1()), cc.getAll(size2())); }
    const MX sub(const Matrix<int>& rr, const Slice& cc) const
    { return sub(rr, cc.getAll(size2())); }
    const MX sub(const Slice& rr, const Matrix<int>& cc) const
    { return sub(rr.getAll(size1()), cc); }
    const MX sub(const Matrix<int>& rr, const Matrix<int>& cc) const;
    const MX sub(int rr, const Matrix<int>& cc) const { return sub(Slice(rr), cc);}

    void setSub(const MX& m, int rr, int cc);
    void setSub(const MX& m, const std::vector<int>& rr, int cc);
    void setSub(const MX& m, int rr, const std::vector<int>& cc);
    void setSub(const MX& m, const std::vector<int>& rr, const std::vector<int>& cc);
    void setSub(const MX& m, const std::vector<int>& rr, Slice cc)
    { setSub(m, rr, cc.getAll(size2()));}
    void setSub(const MX& m, const Matrix<int>& rr, const std::vector<int>& cc);
    void setSub(const MX& m, const Matrix<int>& rr, int cc)
    { setSub(m, rr, std::vector<int>(1, cc)); }
    void setSub(const MX& m, const std::vector<int>& rr, const Matrix<int>& cc);
    void setSub(const MX& m, int rr, const Matrix<int>& cc)
    { setSub(m, std::vector<int>(1, rr), cc); }
    void setSub(const MX& m, const Slice& rr, const Slice& cc);
    void setSub(const MX& m, const Matrix<int>& rr, const Matrix<int>& cc);
    void setSub(const MX& m, const Sparsity& sp, int dummy);

    MX getNZ(int k) const;
    MX getNZ(const std::vector<int>& k) const;
    MX getNZ(const Slice& k) const { return getNZ(k.getAll(size()));}
    MX getNZ(const Matrix<int>& k) const;
    void setNZ(int k, const MX& el);
    void setNZ(const std::vector<int>& k, const MX& el);
    void setNZ(const Slice& k, const MX& m) { setNZ(k.getAll(size()), m);}
    void setNZ(const Matrix<int>& k, const MX& m);
    /// \endcond


    /** \brief Append a matrix vertically (NOTE: only efficient if vector) */
    void append(const MX& y);

    /** \brief Append a matrix horizontally */
    void appendColumns(const MX& y);

    /// \cond SWIGINTERNAL
    /// all binary operations
    MX __add__(const MX& y) const;
    MX __sub__(const MX& y) const;
    MX __mul__(const MX& y) const;
    MX __div__(const MX& y) const;
    MX __lt__(const MX& y) const;
    MX __le__(const MX& y) const;
    MX __eq__(const MX& y) const;
    MX __ne__(const MX& y) const;
    MX __truediv__(const MX& y) const { return __div__(y);}
    MX __pow__(const MX& b) const;
    MX __constpow__(const MX& b) const;
    MX __mrdivide__(const MX& b) const;
    MX __mpower__(const MX& b) const;
    /// \endcond

    MX mul(const MX& y, const Sparsity &sp_z=Sparsity()) const;
    MX mul_full(const MX& y, const Sparsity &sp_z=Sparsity()) const;
    MX inner_prod(const MX& y) const;
    MX outer_prod(const MX& y) const;
    MX constpow(const MX& y) const;
    MX fmin(const MX& y) const;
    MX fmax(const MX& y) const;
    MX fmod(const MX& y) const;
    MX printme(const MX& y) const;
    MX arctan2(const MX& y) const;
    MX logic_and(const MX& y) const;
    MX logic_or(const MX& y) const;
    MX if_else_zero(const MX& y) const;
    MX __copysign__(const MX& y) const;

    // all unary operations
    MX exp() const;
    MX log() const;
    MX log10() const;
    MX sqrt() const;
    MX sin() const;
    MX cos() const;
    MX tan() const;
    MX arcsin() const;
    MX arccos() const;
    MX arctan() const;
    MX floor() const;
    MX ceil() const;
    MX fabs() const;
    MX sign() const;
    MX erfinv() const;
    MX erf() const;
    MX sinh() const;
    MX cosh() const;
    MX tanh() const;
    MX arcsinh() const;
    MX arccosh() const;
    MX arctanh() const;
    MX logic_not() const;

    /** \brief returns itself, but with an assertion attached
    *
    *  If y does not evaluate to 1, a runtime error is raised
    */
    MX attachAssert(const MX& y, const std::string &fail_message="") const;

    /** \brief Set sparse */
    MX setSparse(const Sparsity& sp, bool intersect=false) const;

    /// Make the matrix dense
    void densify(const MX& val = 0);

    /// Lift an expression
    void lift(const MX& x_guess);

    /// Add an expression to the expression if the expression is non-empty, otherwise assign
    void addToSum(const MX& x);

    /// Transpose the matrix
    MX trans() const;

#ifndef SWIG
    /// Transpose the matrix (shorthand)
    MX T() const { return trans();}
#endif

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
