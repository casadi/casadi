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


#ifndef CASADI_MX_HPP
#define CASADI_MX_HPP
#include "shared_object.hpp"
#include "matrix_fwd.hpp"
#include "sx_fwd.hpp"
#include "dm.hpp"
#include "generic_matrix.hpp"
#include "generic_expression.hpp"
#include "generic_type.hpp"
#include "printable.hpp"
#include <vector>
#ifdef CASADI_WITH_THREAD
#ifdef CASADI_WITH_THREAD_MINGW
#include <mingw.mutex.h>
#else // CASADI_WITH_THREAD_MINGW
#include <mutex>
#endif // CASADI_WITH_THREAD_MINGW
#endif //CASADI_WITH_THREAD

namespace casadi {

  /** \brief  Forward declaration

      \identifier{px} */
  class MXNode;
  class Function;
  class SerializingStream;
  class DeserializingStream;

#ifndef SWIG
  struct ConvexifyData {
    std::vector<casadi_int> scc_offset, scc_mapping;
    Sparsity scc_sp;
    Sparsity Hrsp;
    Sparsity Hsp;
    casadi_convexify_config<double> config;
    casadi_int sz_iw;
    casadi_int sz_w;
    casadi_int *iw;
    double* w;
  };
#endif

  /** \brief MX - Matrix expression

      The MX class is used to build up trees made up from MXNodes. It is a more general
      graph representation than the scalar expression, SX, and much less efficient for small
      objects. On the other hand, the class allows much more general operations than does SX,
      in particular matrix valued operations and calls to arbitrary differentiable functions.

      The MX class is designed to have identical syntax with the Matrix<> template class,
      and uses DM (i.e. Matrix<double>) as its internal representation of the values at a node. By keeping
      the syntaxes identical, it is possible to switch from one class to the other,
      as well as inlining MX functions to SXElem functions.

      Note that an operation is always "lazy", making a matrix multiplication will create a
      matrix multiplication node, not perform the actual multiplication.

      \author Joel Andersson
      \date 2010-2011

      \identifier{py} */
  class CASADI_EXPORT MX :
    public SWIG_IF_ELSE(GenericExpressionCommon, GenericExpression<MX>),
    public SWIG_IF_ELSE(PrintableCommon, Printable<MX>),
    public GenericMatrix<MX>,
    public SharedObject {
  public:
    /** \brief Get type name

        \identifier{pz} */
    static std::string type_name() {return "MX";}

    /** \brief  Default constructor

        \identifier{q0} */
    MX();

    /** \brief Create a sparse matrix with all structural zeros

        \identifier{q1} */
    MX(casadi_int nrow, casadi_int ncol);

#ifndef SWIG
    /** \brief Create a sparse matrix with all structural zeros

        \identifier{q2} */
    explicit MX(const std::pair<casadi_int, casadi_int>& rc);
#endif // SWIG

    /** \brief Create a sparse matrix from a sparsity pattern.

        Same as MX::ones(sparsity)

        \identifier{q3} */
    explicit MX(const Sparsity& sp);

    /** \brief Construct matrix with a given sparsity and nonzeros

        \identifier{q4} */
    MX(const Sparsity& sp, const MX& val);

    /** \brief Construct matrix with a given sparsity and a file with nonzeros

        \identifier{q5} */
    MX(const Sparsity& sp, const std::string& fname);


    /** \brief Construct matrix with a given sparsity and nonzeros,

     * configurable in codegen via a pool

        \identifier{2aa} */
    MX(const Matrix<double>& val, const std::string& name);

    /** \brief  Create scalar constant (also implicit type conversion)

        \identifier{q6} */
    MX(double x);

#ifndef SWIG
    /** \brief  Create vector constant (also implicit type conversion)

        \identifier{q7} */
    MX(const std::vector<double> &x);
#endif

    /** \brief  Create sparse matrix constant (also implicit type conversion)

        \identifier{q8} */
    MX(const Matrix<double> &x);

/// \cond INTERNAL
    /** \brief  Destructor

        \identifier{q9} */
    ~MX();
/// \endcond

#ifndef SWIG
/// \cond INTERNAL
    /** \brief  Create from node

        \identifier{qa} */
    static MX create(MXNode* node);

    /** \brief  Create from node (multiple-outputs)

        \identifier{qb} */
    static std::vector<MX> createMultipleOutput(MXNode* node);
/// \endcond

    /** \brief Get the sparsity pattern

        \identifier{qc} */
    const Sparsity& sparsity() const;

    /// \cond INTERNAL
    /// Scalar type
    typedef MX ScalarType;
    /// \endcond

    /// Base class
    typedef GenericMatrix<MX> B;

    /// Expose base class functions
    using B::horzsplit;
    using B::diagsplit;
    using B::vertsplit;
    using B::mtimes;
    using B::repmat;
#endif // SWIG

    /// Returns the truth value of an MX expression
    bool __nonzero__() const;

    /** \brief Get an owning reference to the sparsity pattern

        \identifier{qd} */
    Sparsity get_sparsity() const { return sparsity();}

    /** \brief Get nonzeros as list of scalar MXes
    *
    * Since MX is not a containter, the scalar MXes may be complex
    * When the expression satisfies is_valid_input, the results may be simple
    *
    * For example: 
    * vertcat(x,y).nonzeros()
    * will return {x,y}
    *
    * \sa expr.nz[:]
    *

        \identifier{2bh} */
    std::vector<MX> get_nonzeros() const;

    /** \brief Erase a submatrix (leaving structural zeros in its place)

        Erase rows and/or columns of a matrix

        \identifier{qe} */
    void erase(const std::vector<casadi_int>& rr, const std::vector<casadi_int>& cc,
                bool ind1=false);

    /** \brief Erase a submatrix (leaving structural zeros in its place)

        Erase elements of a matrix

        \identifier{qf} */
    void erase(const std::vector<casadi_int>& rr, bool ind1=false);

    /** \brief Enlarge matrix

        Make the matrix larger by inserting empty rows and columns, keeping the existing non-zeros

        \identifier{qg} */
    void enlarge(casadi_int nrow, casadi_int ncol,
                  const std::vector<casadi_int>& rr, const std::vector<casadi_int>& cc,
                  bool ind1=false);

    MX operator-() const;

#ifndef SWIG
    /// \cond INTERNAL
    ///@{
    /** \brief  Access a member of the node

        \identifier{qh} */
    MXNode* operator->();

    /** \brief  Const access a member of the node

        \identifier{qi} */
    const MXNode* operator->() const;
    ///@}
    /// \endcond
#endif // SWIG

    /** \brief Get the nth dependency as MX

        \identifier{qj} */
    MX dep(casadi_int ch=0) const;

    /** \brief  Number of outputs

        \identifier{qk} */
    casadi_int n_out() const;

    /** \brief  Get an output

        \identifier{ql} */
    MX get_output(casadi_int oind) const;

    /** \brief Get the number of dependencies of a binary SXElem

        \identifier{qm} */
    casadi_int n_dep() const;

    /// Get the name.
    std::string name() const;

    /// Get the value (only for scalar constant nodes)
    explicit operator double() const;

    /// Get the value (only for constant nodes)
    explicit operator Matrix<double>() const;

    /// Check if symbolic
    bool is_symbolic() const;

    /// Check if constant
    bool is_constant() const;

    /// Check if evaluation
    bool is_call() const;

    /// Get function - only valid when is_call() is true
    Function which_function() const;

    /// Check if evaluation output
    bool is_output() const;

    /** \brief  Check if a multiple output node

        \identifier{284} */
    bool has_output() const;

    /// Get the index of evaluation output - only valid when is_output() is true
    casadi_int which_output() const;

    /// Is it a certain operation
    bool is_op(casadi_int op) const;

    /// Check if multiplication
    bool is_multiplication() const;

    /// Check if commutative operation
    bool is_commutative() const;

    /// Check if norm
    bool is_norm() const;

    /** \brief Check if matrix can be used to define function inputs.

        Valid inputs for MXFunctions are combinations of Reshape, concatenations and SymbolicMX

        \identifier{qn} */
    bool is_valid_input() const;

    /** \brief Get the number of primitives for MXFunction inputs/outputs

        \identifier{qo} */
    casadi_int n_primitives() const;

    /** \brief Get primitives

        \identifier{qp} */
    std::vector<MX> primitives() const;

    /// @{
    /** \brief Split up an expression along symbolic primitives

        \identifier{qq} */
    std::vector<MX> split_primitives(const MX& x) const;
    std::vector<SX> split_primitives(const SX& x) const;
    std::vector<DM> split_primitives(const DM& x) const;
    /// @}

    /// @{
    /** \brief Join an expression along symbolic primitives

        \identifier{qr} */
    MX join_primitives(const std::vector<MX>& v) const;
    SX join_primitives(const std::vector<SX>& v) const;
    DM join_primitives(const std::vector<DM>& v) const;
    /// @}

    /// \cond INTERNAL
    /** \brief Detect duplicate symbolic expressions

        If there are symbolic primitives appearing more than once, the function will return
        true and the names of the duplicate expressions will be passed to casadi_warning.
        Note: Will mark the node using MX::set_temp.
        Make sure to call reset_input() after usage.

        \identifier{qs} */
    bool has_duplicates() const;

    /** \brief Reset the marker for an input expression

        \identifier{qt} */
    void reset_input() const;
  /// \endcond

    /** \brief  check if identity

        \identifier{qu} */
    bool is_eye() const;

    /** \brief  check if zero (note that false negative answers are possible)

        \identifier{qv} */
    bool is_zero() const;

    /** \brief  check if zero (note that false negative answers are possible)

        \identifier{qw} */
    bool is_one() const;

    /** \brief  check if zero (note that false negative answers are possible)

        \identifier{qx} */
    bool is_minus_one() const;

    /** \brief  Is the expression a transpose?

        \identifier{qy} */
    bool is_transpose() const;

    /// Checks if expression does not contain NaN or Inf
    bool is_regular() const;

    /// Is binary operation
    bool is_binary() const;

    /// Is unary operation
    bool is_unary() const;

    /// Get operation type
    casadi_int op() const;

    /** Obtain information about node */
    Dict info() const;

    /** \brief Serialize an object

        \identifier{qz} */
    void serialize(SerializingStream& s) const;

    /** \brief Deserialize with type disambiguation

        \identifier{r0} */
    static MX deserialize(DeserializingStream& s);

    /// \cond INTERNAL
    /// Get the temporary variable
    casadi_int get_temp() const;

    /// Set the temporary variable
    void set_temp(casadi_int t) const;
    /// \endcond

    ///@{
    /** \brief  Create nodes by their ID

        \identifier{r1} */
    static MX binary(casadi_int op, const MX &x, const MX &y);
    static MX unary(casadi_int op, const MX &x);
    ///@}

    ///@{
    /** \brief  create a matrix with all inf

        \identifier{r2} */
    static MX inf(const Sparsity& sp);
    static MX inf(casadi_int nrow=1, casadi_int ncol=1);
    static MX inf(const std::pair<casadi_int, casadi_int>& rc);
    ///@}

    ///@{
    /** \brief  create a matrix with all nan

        \identifier{r3} */
    static MX nan(const Sparsity& sp);
    static MX nan(casadi_int nrow=1, casadi_int ncol=1);
    static MX nan(const std::pair<casadi_int, casadi_int>& rc);
    ///@}

    /** \brief  Identity matrix

        \identifier{r4} */
    static MX eye(casadi_int n);

#ifndef SWIG
    /// Get a const pointer to the node
    MXNode* get() const;
#endif // SWIG

    ///@{
    /// Get a submatrix, single argument
    void get(MX& SWIG_OUTPUT(m), bool ind1, const Slice& rr) const;
    void get(MX& SWIG_OUTPUT(m), bool ind1, const Matrix<casadi_int>& rr) const;
    void get(MX& SWIG_OUTPUT(m), bool ind1, const Sparsity& sp) const;
    void get(MX& SWIG_OUTPUT(m), bool ind1, const MX& rr) const;
    void get(MX& SWIG_OUTPUT(m), bool ind1, const casadi_int rr) const {
      get(m, ind1, Matrix<casadi_int>(rr));
    }
    ///@}

    /// Get a submatrix, two arguments
    ///@{
    void get(MX& SWIG_OUTPUT(m), bool ind1, const Slice& rr, const Slice& cc) const;
    void get(MX& SWIG_OUTPUT(m), bool ind1, const Slice& rr, const Matrix<casadi_int>& cc) const;
    void get(MX& SWIG_OUTPUT(m), bool ind1, const Slice& rr, casadi_int cc) const {
      get(m, ind1, rr, Matrix<casadi_int>(cc));
    }
    void get(MX& SWIG_OUTPUT(m), bool ind1, const Matrix<casadi_int>& rr, const Slice& cc) const;
    void get(MX& SWIG_OUTPUT(m), bool ind1, casadi_int rr, const Slice& cc) const {
      get(m, ind1, Matrix<casadi_int>(rr), cc);
    }
    void get(MX& SWIG_OUTPUT(m), bool ind1, const Matrix<casadi_int>& rr,
                                            const Matrix<casadi_int>& cc) const;
    void get(MX& SWIG_OUTPUT(m), bool ind1, casadi_int rr,
                                            casadi_int cc) const {
      get(m, ind1, Matrix<casadi_int>(rr), Matrix<casadi_int>(cc));
    }
    void get(MX& SWIG_OUTPUT(m), bool ind1, const MX& rr, const Slice& cc) const;
    void get(MX& SWIG_OUTPUT(m), bool ind1, const Slice& rr, const MX& cc) const;
    void get(MX& SWIG_OUTPUT(m), bool ind1, const MX& rr, const MX& cc) const;
    ///@}

    ///@{
    /// Set a submatrix, single argument
    void set(const MX& m, bool ind1, const Slice& rr);
    void set(const MX& m, bool ind1, const Matrix<casadi_int>& rr);
    void set(const MX& m, bool ind1, const Sparsity& sp);
    ///@}

    ///@{
    /// Set a submatrix, two arguments
    void set(const MX& m, bool ind1, const Slice& rr, const Slice& cc);
    void set(const MX& m, bool ind1, const Slice& rr, const Matrix<casadi_int>& cc);
    void set(const MX& m, bool ind1, const Matrix<casadi_int>& rr, const Slice& cc);
    void set(const MX& m, bool ind1, const Matrix<casadi_int>& rr, const Matrix<casadi_int>& cc);
    ///@}

    ///@{
    /// Get a set of nonzeros
    void get_nz(MX& SWIG_OUTPUT(m), bool ind1, const Slice& kk) const;
    void get_nz(MX& SWIG_OUTPUT(m), bool ind1, const Matrix<casadi_int>& kk) const;
    void get_nz(MX& SWIG_OUTPUT(m), bool ind1, const MX& kk) const;
    void get_nz(MX& SWIG_OUTPUT(m), bool ind1, casadi_int kk) const {
      get_nz(m, ind1, Matrix<casadi_int>(kk));
    }
    void get_nz(MX& SWIG_OUTPUT(m), bool ind1, const MX& inner, const Slice& outer) const;
    void get_nz(MX& SWIG_OUTPUT(m), bool ind1, const Slice& inner, const MX& outer) const;
    void get_nz(MX& SWIG_OUTPUT(m), bool ind1, const MX& inner, const MX& outer) const;
    ///@}

    ///@{
    /// Set a set of nonzeros
    void set_nz(const MX& m, bool ind1, const Slice& kk);
    void set_nz(const MX& m, bool ind1, const Matrix<casadi_int>& kk);
    void set_nz(const MX& m, bool ind1, const MX& kk);
    void set_nz(const MX& m, bool ind1, casadi_int kk) { set_nz(m, ind1, Matrix<casadi_int>(kk)); }
    ///@}

    ///@{
    /** \brief Computes an einstein dense tensor contraction

        Computes the product:
        C_c = A_a + B_b
          where a b c are index/einstein notation in an encoded form

        For example, an matrix-matrix product may be written as:
        C_ij = A_ik B_kj

        The encoded form uses strictly negative numbers to indicate labels.
        For the above example, we would have:
        a {-1, -3} b {-3, -2} c {-1 -2}

        \identifier{r5} */
    static MX einstein(const MX& A, const MX& B, const MX& C,
      const std::vector<casadi_int>& dim_a, const std::vector<casadi_int>& dim_b,
      const std::vector<casadi_int>& dim_c,
      const std::vector<casadi_int>& a, const std::vector<casadi_int>& b,
      const std::vector<casadi_int>& c);

    static MX einstein(const MX& A, const MX& B,
      const std::vector<casadi_int>& dim_a, const std::vector<casadi_int>& dim_b,
      const std::vector<casadi_int>& dim_c,
      const std::vector<casadi_int>& a, const std::vector<casadi_int>& b,
      const std::vector<casadi_int>& c);
    ///@}

#ifndef SWIG
    /// \cond CLUTTER
    ///@{
    /// Functions called by friend functions defined for GenericExpression
    static bool is_equal(const MX& x, const MX& y, casadi_int depth=0);
    static MX mmin(const MX &x);
    static MX mmax(const MX &x);
    ///@}

    ///@{
    /// Functions called by friend functions defined for SparsityInterface
    static MX horzcat(const std::vector<MX>& x);
    static MX diagcat(const std::vector<MX>& x);
    static MX vertcat(const std::vector<MX>& x);
    static std::vector<MX> horzsplit(const MX& x, const std::vector<casadi_int>& offset);
    static std::vector<MX> diagsplit(const MX& x, const std::vector<casadi_int>& offset1,
                                     const std::vector<casadi_int>& offset2);
    static std::vector<MX> vertsplit(const MX& x, const std::vector<casadi_int>& offset);
    static MX blockcat(const std::vector< std::vector<MX > > &v);
    static MX mtimes(const MX& x, const MX& y);
    static MX mac(const MX& x, const MX& y, const MX& z);
    static MX reshape(const MX& x, casadi_int nrow, casadi_int ncol);
    static MX reshape(const MX& x, const Sparsity& sp);
    static MX sparsity_cast(const MX& x, const Sparsity& sp);
    static MX kron(const MX& x, const MX& b);
    static MX repmat(const MX& x, casadi_int n, casadi_int m=1);
    ///@}

    ///@{
    /// Functions called by friend functions defined for GenericMatrix
    static MX jacobian(const MX& f, const MX& x, const Dict& opts = Dict());
    static MX hessian(const MX& f, const MX& x, const Dict& opts = Dict());
    static MX hessian(const MX& f, const MX& x, MX& g, const Dict& opts = Dict());
    static std::vector<std::vector<MX> >
    forward(const std::vector<MX> &ex,
            const std::vector<MX> &arg,
            const std::vector<std::vector<MX> > &v,
            const Dict& opts = Dict());
    static std::vector<std::vector<MX> >
    reverse(const std::vector<MX> &ex,
            const std::vector<MX> &arg,
            const std::vector<std::vector<MX> > &v,
            const Dict& opts = Dict());
    static std::vector<bool> which_depends(const MX &expr, const MX &var,
        casadi_int order=1, bool tr=false);
    static Sparsity jacobian_sparsity(const MX& f, const MX& x);
    static MX substitute(const MX& ex, const MX& v, const MX& vdef);
    static std::vector<MX> substitute(const std::vector<MX> &ex,
                                         const std::vector<MX> &v,
                                         const std::vector<MX> &vdef);
    static void substitute_inplace(const std::vector<MX>& v,
                                  std::vector<MX>& vdef,
                                  std::vector<MX>& ex, bool reverse);
    static MX solve(const MX& a, const MX& b);
    static MX solve(const MX& a, const MX& b, const std::string& lsolver,
                    const Dict& dict = Dict());
    static MX inv_minor(const MX& A);
    static MX inv_node(const MX& A);
    static MX inv(const MX& A, const std::string& lsolver="qr", const Dict& dict = Dict());
    static MX pinv(const MX& A, const std::string& lsolver="qr",
               const Dict& dict = Dict());
    static MX expm_const(const MX& A, const MX& t);
    static MX expm(const MX& A);
    static casadi_int n_nodes(const MX& x);
    static std::string print_operator(const MX& x, const std::vector<std::string>& args);
    static void extract(std::vector<MX>& ex, std::vector<MX>& v,
      std::vector<MX>& vdef, const Dict& opts = Dict());
    static void shared(std::vector<MX>& ex, std::vector<MX>& v,
      std::vector<MX>& vdef, const std::string& v_prefix, const std::string& v_suffix);
    static MX if_else(const MX& cond, const MX& if_true,
                      const MX& if_false, bool short_circuit=false);
    static MX conditional(const MX& ind, const std::vector<MX> &x, const MX& x_default,
                          bool short_circuit=false);
    static bool depends_on(const MX& x, const MX& arg);
    static bool contains_all(const std::vector<MX>& v, const std::vector<MX> &n);
    static bool contains_any(const std::vector<MX>& v, const std::vector<MX> &n);
    static MX simplify(const MX& x);
    static MX dot(const MX& x, const MX& y);
    static MX mrdivide(const MX& a, const MX& b);
    static MX mldivide(const MX& a, const MX& b);
    static MX norm_2(const MX& x);
    static MX norm_fro(const MX& x);
    static MX norm_1(const MX& x);
    static MX norm_inf(const MX& x);
    static MX unite(const MX& A, const MX& B);
    static MX trace(const MX& x);
    static MX diag(const MX& x);
    static MX sum2(const MX& x);
    static MX sum1(const MX& x);
    static MX polyval(const MX& p, const MX& x);
    static MX det(const MX& x);
    static std::vector<MX> symvar(const MX& x);
    static MX nullspace(const MX& A);
    static MX repsum(const MX& x, casadi_int n, casadi_int m=1);
    static MX densify(const MX& x, const MX& val=0);
    static MX _bilin(const MX& A, const MX& x, const MX& y);
    static MX _rank1(const MX& A, const MX& alpha, const MX& x, const MX& y);
    static MX project(const MX& x, const Sparsity& sp, bool intersect=false);
    static MX cumsum(const MX &x, casadi_int axis=-1);
    static MX _logsumexp(const MX& x);
    static std::vector<MX> cse(const std::vector<MX>& e);
    static void extract_parametric(const MX &expr, const MX& par,
            MX& expr_ret, std::vector<MX>& symbols, std::vector<MX>& parametric, const Dict& opts);
    static void separate_linear(const MX &expr,
      const MX &sym_lin, const MX &sym_const,
      MX& expr_const, MX& expr_lin, MX& expr_nonlin);
    ///@}

    ///@{
    /// Functions called by friend functions defined for this class
    static MX find(const MX& x);
    static MX low(const MX& v, const MX& p, const Dict& options = Dict());
    static MX graph_substitute(const MX& x, const std::vector<MX> &v,
                               const std::vector<MX> &vdef);
    static MX graph_substitute(const MX& x, const std::vector<MX> &v,
                               const std::vector<MX> &vdef, bool& updated);
    static std::vector<MX> graph_substitute(const std::vector<MX> &ex,
                                            const std::vector<MX> &v,
                                            const std::vector<MX> &vdef);
    static std::vector<MX> graph_substitute(const std::vector<MX> &ex,
                                            const std::vector<MX> &v,
                                            const std::vector<MX> &vdef,
                                            bool& updated);
    static MX matrix_expand(const MX& e, const std::vector<MX> &boundary,
                            const Dict& options);
    static std::vector<MX> matrix_expand(const std::vector<MX>& e,
                                         const std::vector<MX>& boundary,
                                         const Dict& options);
    static MX lift(const MX& x, const MX& x_guess);
    static DM evalf(const MX& m);
    static MX bspline(const MX& x,
            const DM& coeffs,
            const std::vector< std::vector<double> >& knots,
            const std::vector<casadi_int>& degree,
            casadi_int m,
            const Dict& opts = Dict());
    static MX bspline(const MX& x, const MX& coeffs,
            const std::vector< std::vector<double> >& knots,
            const std::vector<casadi_int>& degree,
            casadi_int m,
            const Dict& opts = Dict());
    static MX convexify(const MX& H, const Dict& opts = Dict());
    static MX stop_diff(const MX& expr, casadi_int order);
    static MX stop_diff(const MX& expr, const MX& var, casadi_int order);
    static std::vector<MX> difference(const std::vector<MX>& a, const std::vector<MX>& b);
    ///@}
    /// \endcond

#endif // SWIG

    static DM bspline_dual(const std::vector<double>& x,
            const std::vector< std::vector<double> >& knots,
            const std::vector<casadi_int>& degree,
            const Dict& opts = Dict());

    /** \brief Low-level access to inlined linear interpolation
     *
     * Usually, you want to be using 'interpolant' instead.
     *
     * Accepts lookup_mode option.

        \identifier{r6} */
    static MX interpn_linear(const std::vector<MX>& x, const MX& v, const std::vector<MX>& xq,
      const Dict& opts=Dict());

    MX printme(const MX& b) const;

#if !defined(SWIG) || defined(DOXYGEN)
/**
\addtogroup expression_tools
@{
*/
    /** \brief Find first nonzero, returned as row index

     * If failed, returns the number of rows

        \identifier{r7} */
    inline friend MX find(const MX& x) {
      return MX::find(x);
    }

    /** \brief Find first nonzero

     * If failed, returns the number of rows

        \identifier{r8} */
    inline friend MX low(const MX& v, const MX& p, const Dict& options=Dict()) {
      return MX::low(v, p, options);
    }

    /** \brief Substitute single expression in graph

     * Substitute variable v with expression vdef in an expression ex, preserving nodes

        \identifier{r9} */
    inline friend MX graph_substitute(const MX& ex, const std::vector<MX> &v,
                                      const std::vector<MX> &vdef) {
      return MX::graph_substitute(ex, v, vdef);
    }

    inline friend MX graph_substitute(const MX& ex, const std::vector<MX> &v,
                                      const std::vector<MX> &vdef, bool& updated) {
      return MX::graph_substitute(ex, v, vdef, updated);
    }

    /** \brief Substitute multiple expressions in graph

     * Substitute variable var with expression expr in
     * multiple expressions, preserving nodes

        \identifier{ra} */
    inline friend std::vector<MX>
      graph_substitute(const std::vector<MX> &ex,
                       const std::vector<MX> &v,
                       const std::vector<MX> &vdef) {
      return MX::graph_substitute(ex, v, vdef);
    }

    inline friend std::vector<MX>
      graph_substitute(const std::vector<MX> &ex,
                       const std::vector<MX> &v,
                       const std::vector<MX> &vdef,
                       bool& updated) {
      return MX::graph_substitute(ex, v, vdef, updated);
    }

    /** \brief Expand MX graph to SXFunction call
     *
     *  Expand the given expression e, optionally
     *  supplying expressions contained in it at which expansion should stop.
     *
        \identifier{rb} */
    inline friend MX
      matrix_expand(const MX& e, const std::vector<MX> &boundary = std::vector<MX>(),
        const Dict& options = Dict()) {
      return MX::matrix_expand(e, boundary, options);
    }

    /** \brief Expand MX graph to SXFunction call
     *
     *  Expand the given expression e, optionally
     *  supplying expressions contained in it at which expansion should stop.
     *
        \identifier{rc} */
    inline friend std::vector<MX>
      matrix_expand(const std::vector<MX>& e,
                    const std::vector<MX> &boundary = std::vector<MX>(),
                    const Dict& options = Dict()) {
      return MX::matrix_expand(e, boundary, options);
    }


    inline friend MX bspline(const MX& x,
            const DM& coeffs,
            const std::vector< std::vector<double> >& knots,
            const std::vector<casadi_int>& degree,
            casadi_int m,
            const Dict& opts = Dict()) {
      return MX::bspline(x, coeffs, knots, degree, m, opts);
    }

    inline friend MX bspline(const MX& x, const MX& coeffs,
            const std::vector< std::vector<double> >& knots,
            const std::vector<casadi_int>& degree,
            casadi_int m,
            const Dict& opts = Dict()) {
      return MX::bspline(x, coeffs, knots, degree, m, opts);
    }

    inline friend DM bspline_dual(const std::vector<double>& x,
            const std::vector< std::vector<double> >& knots,
            const std::vector<casadi_int>& degree,
            const Dict& opts = Dict()) {
      return MX::bspline_dual(x, knots, degree, opts);
    }

    inline friend MX convexify(const MX& H,
            const Dict& opts = Dict()) {
      return MX::convexify(H, opts);
    }

    /** \brief Lift the expression

     * Experimental feature
     *
        \identifier{rd} */
    inline friend MX lift(const MX& x, const MX& x_guess) {
      return MX::lift(x, x_guess);
    }

    /** \brief Inverse node
     *
        \identifier{re} */
    inline friend MX inv_node(const MX& x) {
      return MX::inv_node(x);
    }

    /** \brief Evaluates the expression numerically
    *
    * An error is raised when the expression contains symbols

        \identifier{rf} */
    inline friend DM evalf(const MX& expr) {
      return MX::evalf(expr);
    }

    /** \brief Stop derivatives of an expression wrt to all its symbolic variables

        \identifier{25l} */
    inline friend MX stop_diff(const MX& expr, casadi_int order) {
      return MX::stop_diff(expr, order);
    }

    /** \brief Stop first derivatives of an expression wrt to all its symbolic variables
     * 
     * \seealso stop_diff

        \identifier{25m} */
    inline friend MX no_grad(const MX& expr) {
      return MX::stop_diff(expr, 1);
    }

    /** \brief Stop second derivatives of an expression wrt to all its symbolic variables
     * 
     * \seealso stop_diff

        \identifier{25n} */
    inline friend MX no_hess(const MX& expr) {
      return MX::stop_diff(expr, 2);
    }


    /** \brief Stop derivatives of an expression wrt to a select set of symbolic variables

        \identifier{25o} */
    inline friend MX stop_diff(const MX& expr, const MX& var, casadi_int order) {
      return MX::stop_diff(expr, var, order);
    }

    /** \bried Return all elements of a that do not occur in b, preserving order */
    inline friend std::vector<MX> difference(const std::vector<MX>& a, const std::vector<MX>& b) {
      return MX::difference(a, b);
    }

/** @} */
#endif // SWIG

    /** \brief returns itself, but with an assertion attached
    *
    *  If y does not evaluate to 1, a runtime error is raised

        \identifier{rg} */
    MX attachAssert(const MX& y, const std::string& fail_message="") const;

    /** \brief Monitor an expression

    * Returns itself, but with the side effect of printing the nonzeros along with a comment

        \identifier{rh} */
    MX monitor(const std::string& comment) const;

    /// Transpose the matrix
    MX T() const;

    /** \brief Get an IM representation of a GetNonzeros or SetNonzeros node

        \identifier{ri} */
    Matrix<casadi_int> mapping() const;

    /** \brief Set or reset the depth to which equalities are being checked for simplifications

        \identifier{rj} */
    static void set_max_depth(casadi_int eq_depth=1);

    /** \brief Get the depth to which equalities are being checked for simplifications

        \identifier{rk} */
    static casadi_int get_max_depth();

    /// Check if a particular cast is allowed
    static bool test_cast(const SharedObjectInternal* ptr);

    /** \brief Get function inputs

        \identifier{rl} */
    static std::vector<MX> get_input(const Function& f);

    /** \brief Get free variables

        \identifier{rm} */
    static std::vector<MX> get_free(const Function& f);

    /// Readability typedef
    typedef std::map<std::string, MX> MXDict;

    /** \brief Evaluate the MX node with new symbolic dependencies

        \identifier{rn} */
    void eval_mx(const std::vector<MX>& arg, std::vector<MX>& SWIG_OUTPUT(res)) const;

#ifndef SWIG
    ///@{
    /** \brief Called from MXFunction

        \identifier{ro} */
    void ad_forward(const std::vector<std::vector<MX> >& fseed,
                    std::vector<std::vector<MX> >& fsens) const;
    void ad_reverse(const std::vector<std::vector<MX> >& aseed,
                    std::vector<std::vector<MX> >& asens) const;
    ///@}

    /// Construct constant matrix with a given sparsity and values
    MX(const Sparsity& sp, double val, bool dummy);

    // Create matrix symbolic primitive
    static MX _sym(const std::string& name, const Sparsity& sp);

#ifdef CASADI_WITH_THREADSAFE_SYMBOLICS
    static std::mutex& get_mutex_temp() { return mutex_temp; }
    static std::mutex mutex_temp;
#endif //CASADI_WITH_THREADSAFE_SYMBOLICS
  private:

    /// Create an expression from a node: extra dummy arguments to avoid ambiguity for 0/NULL
    MX(MXNode* node, bool dummy1, bool dummy2, bool dummy3, bool dummy4);

    // Depth when checking equalities
    static casadi_int eq_depth_;

#endif // SWIG
  };


  ///@{
  /// Readability typedefs
  typedef std::vector<MX> MXVector;
  typedef std::initializer_list<MX> MXIList;
  typedef std::vector<MXVector> MXVectorVector;
  typedef std::map<std::string, MX> MXDict;
  ///@}

} // namespace casadi

#endif // CASADI_MX_HPP
