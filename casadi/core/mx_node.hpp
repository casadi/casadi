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


#ifndef CASADI_MX_NODE_HPP
#define CASADI_MX_NODE_HPP

#include "mx.hpp"
#include "shared_object_internal.hpp"
#include "sx_elem.hpp"
#include "calculus.hpp"
#include "code_generator.hpp"
#include "linsol.hpp"
#include <vector>
#include <stack>

namespace casadi {

  class SerializingStream;
  class DeserializingStream;

  /** \brief Node class for MX objects

      \author Joel Andersson
      \date 2010
      Internal class.

      \identifier{1qb} */
  class CASADI_EXPORT MXNode : public SharedObjectInternal {
    friend class MX;

  public:
    /// Constructor
    MXNode();

    /** \brief  Destructor

        \identifier{1qc} */
    ~MXNode() override=0;

    /** \brief Check the truth value of this node

        \identifier{1qd} */
    virtual bool __nonzero__() const;

    /** \brief Check if identically zero

        \identifier{1qe} */
    virtual bool is_zero() const { return false;}

    /** \brief Check if identically one

        \identifier{1qf} */
    virtual bool is_one() const { return false;}

    /** \brief Check if identically  minus one

        \identifier{1qg} */
    virtual bool is_minus_one() const { return false;}

    /** \brief Check if a certain value

        \identifier{1qh} */
    virtual bool is_value(double val) const { return false;}

    /** \brief Check if identity matrix

        \identifier{1qi} */
    virtual bool is_eye() const { return false;}

    /** \brief Check if unary operation

        \identifier{1qj} */
    virtual bool is_unary() const { return false;}

    /** \brief Check if binary operation

        \identifier{1qk} */
    virtual bool is_binary() const { return false;}

    /** \brief Find out which nodes can be inlined

        \identifier{1ql} */
    void can_inline(std::map<const MXNode*, casadi_int>& nodeind) const;

    /** \brief Print compact

        \identifier{1qm} */
    std::string print_compact(std::map<const MXNode*, casadi_int>& nodeind,
                             std::vector<std::string>& intermed) const;

    /** \brief  Print expression

        \identifier{1qn} */
    virtual std::string disp(const std::vector<std::string>& arg) const = 0;

    /** \brief Add a dependent function

        \identifier{1qo} */
    virtual void add_dependency(CodeGenerator& g) const {}

    /** \brief Is reference counting needed in codegen?

        \identifier{1qp} */
    virtual bool has_refcount() const { return false;}

    /** \brief Codegen incref

        \identifier{1qq} */
    virtual void codegen_incref(CodeGenerator& g, std::set<void*>& added) const {}

    /** \brief Codegen decref

        \identifier{1qr} */
    virtual void codegen_decref(CodeGenerator& g, std::set<void*>& added) const {}

    /** \brief Generate code for the operation

        \identifier{1qs} */
    virtual void generate(CodeGenerator& g,
                          const std::vector<casadi_int>& arg,
                          const std::vector<casadi_int>& res) const;

    /** \brief  Evaluate numerically

        \identifier{1qt} */
    virtual int eval(const double** arg, double** res, casadi_int* iw, double* w) const;

    /** \brief  Evaluate symbolically (SX)

        \identifier{1qu} */
    virtual int eval_sx(const SXElem** arg, SXElem** res, casadi_int* iw, SXElem* w) const;

    /** \brief  Evaluate symbolically (MX)

        \identifier{1qv} */
    virtual void eval_mx(const std::vector<MX>& arg, std::vector<MX>& res) const;

    /** \brief Calculate forward mode directional derivatives

        \identifier{1qw} */
    virtual void ad_forward(const std::vector<std::vector<MX> >& fseed,
                         std::vector<std::vector<MX> >& fsens) const;

    /** \brief Calculate reverse mode directional derivatives

        \identifier{1qx} */
    virtual void ad_reverse(const std::vector<std::vector<MX> >& aseed,
                         std::vector<std::vector<MX> >& asens) const;

    /** \brief  Propagate sparsity forward

        \identifier{1qy} */
    virtual int sp_forward(const bvec_t** arg, bvec_t** res, casadi_int* iw, bvec_t* w) const;

    /** \brief  Propagate sparsity backwards

        \identifier{1qz} */
    virtual int sp_reverse(bvec_t** arg, bvec_t** res, casadi_int* iw, bvec_t* w) const;

    /** \brief  Get the name

        \identifier{1r0} */
    virtual const std::string& name() const;

    /** \brief Get name of public class

        \identifier{1r1} */
    std::string class_name() const override;

    /** \brief  Print a description

        \identifier{1r2} */
    void disp(std::ostream& stream, bool more) const override;

    /** \brief  Check if valid function input

        \identifier{1r3} */
    virtual bool is_valid_input() const { return false;}

    /** \brief Get the number of symbolic primitives

        \identifier{1r4} */
    virtual casadi_int n_primitives() const;

    /** \brief Get symbolic primitives

        \identifier{1r5} */
    virtual void primitives(std::vector<MX>::iterator& it) const;

    /** \brief Split up an expression along symbolic primitives

        \identifier{1r6} */
    virtual void split_primitives(const MX& x, std::vector<MX>::iterator& it) const;

    /** \brief Join an expression along symbolic primitives

        \identifier{1r7} */
    virtual MX join_primitives(std::vector<MX>::const_iterator& it) const;

    /** \brief Detect duplicate symbolic expressions

        \identifier{1r8} */
    virtual bool has_duplicates() const;

    /** \brief Reset the marker for an input expression

        \identifier{1r9} */
    virtual void reset_input() const;

    /** \brief  Check if evaluation output

        \identifier{1ra} */
    virtual bool is_output() const {return false;}

    /** \brief  Check if a multiple output node

        \identifier{1rb} */
    virtual bool has_output() const {return false;}

    /** \brief  Get function output

        \identifier{1rc} */
    virtual casadi_int which_output() const;

    /** \brief  Get called function

        \identifier{1rd} */
    virtual const Function& which_function() const;

    /** \brief Get the operation

        \identifier{1re} */
    virtual casadi_int op() const = 0;

    /** Obtain information about node */
    virtual Dict info() const;

    /** \brief Serialize an object

        \identifier{1rf} */
    void serialize(SerializingStream& s) const;

    /** \brief Serialize an object without type information

        \identifier{1rg} */
    virtual void serialize_body(SerializingStream& s) const;

    /** \brief Serialize type information
     *
     * Information needed to unambiguously find the (lowest level sub)class,
     * such that its deserializing constructor can be called.

        \identifier{1rh} */
    virtual void serialize_type(SerializingStream& s) const;

    /** \brief Deserialize with type disambiguation
     *
     * Uses the information encoded with serialize_type to unambiguously find the (lowest level sub)class,
     * and calls its deserializing constructor.

        \identifier{1ri} */
    static MXNode* deserialize(DeserializingStream& s);

    /** \brief Check if two nodes are equivalent up to a given depth

        \identifier{1rj} */
    static bool is_equal(const MXNode* x, const MXNode* y, casadi_int depth);
    virtual bool is_equal(const MXNode* node, casadi_int depth) const { return false;}

    /** \brief Get equality checking depth

        \identifier{1rk} */
    inline static bool maxDepth() { return MX::get_max_depth();}

    /** \brief Checks if two nodes have the same operation and have

     * equivalent dependencies up to a given depth

        \identifier{1rl} */
    bool sameOpAndDeps(const MXNode* node, casadi_int depth) const;

    /** \brief  dependencies - functions that have to be evaluated before this one

        \identifier{1rm} */
    const MX& dep(casadi_int ind=0) const { return dep_.at(ind);}

    /** \brief  Number of dependencies

        \identifier{1rn} */
    casadi_int n_dep() const;

    /** \brief  Number of outputs

        \identifier{1ro} */
    virtual casadi_int nout() const { return 1;}

    /** \brief  Get an output

        \identifier{1rp} */
    virtual MX get_output(casadi_int oind) const;

    /// Get the sparsity
    const Sparsity& sparsity() const { return sparsity_;}

    /// Get the sparsity of output oind
    virtual const Sparsity& sparsity(casadi_int oind) const;

    template<class T>
    bool matches_sparsity(const std::vector<T>& arg) const {
        for (casadi_int i=0;i<dep_.size();++i) {
            if (dep_[i].sparsity()!=arg[i].sparsity()) {
                return false;
            }
        }
        return true;
    }

    /// Get shape
    casadi_int numel() const { return sparsity().numel(); }
    casadi_int nnz(casadi_int i=0) const { return sparsity(i).nnz(); }
    casadi_int size1() const { return sparsity().size1(); }
    casadi_int size2() const { return sparsity().size2(); }
    std::pair<casadi_int, casadi_int> size() const { return sparsity().size();}

    // Get IO index
    virtual casadi_int ind() const;

    // Get IO segment
    virtual casadi_int segment() const;

    // Get IO offset
    virtual casadi_int offset() const;

    /// Set the sparsity
    void set_sparsity(const Sparsity& sparsity);

    /** \brief Get required length of arg field

        \identifier{1rq} */
    virtual size_t sz_arg() const { return n_dep();}

    /** \brief Get required length of res field

        \identifier{1rr} */
    virtual size_t sz_res() const { return nout();}

    /** \brief Get required length of iw field

        \identifier{1rs} */
    virtual size_t sz_iw() const { return 0;}

    /** \brief Get required length of w field

        \identifier{1rt} */
    virtual size_t sz_w() const { return 0;}

    /// Set unary dependency
    void set_dep(const MX& dep);

    /// Set binary dependencies
    void set_dep(const MX& dep1, const MX& dep2);

    /// Set ternary dependencies
    void set_dep(const MX& dep1, const MX& dep2, const MX& dep3);

    /// Set multiple dependencies
    void set_dep(const std::vector<MX>& dep);

    /// Convert scalar to matrix
    inline static MX to_matrix(const MX& x, const Sparsity& sp) {
      if (x.size()==sp.size()) {
        return x;
      } else {
        return MX(sp, x);
      }
    }

    /// Get the value (only for scalar constant nodes)
    virtual double to_double() const;

    /// Get the value (only for constant nodes)
    virtual DM get_DM() const;

    /// Can the operation be performed inplace (i.e. overwrite the result)
    virtual casadi_int n_inplace() const { return 0;}

    /// Get an IM representation of a GetNonzeros or SetNonzeros node
    virtual Matrix<casadi_int> mapping() const;

    /// Create a horizontal concatenation node
    virtual MX get_horzcat(const std::vector<MX>& x) const;

    /// Create a horizontal split node
    virtual std::vector<MX> get_horzsplit(const std::vector<casadi_int>& output_offset) const;

    /// Create a repeated matrix node
    virtual MX get_repmat(casadi_int m, casadi_int n) const;

    /// Create a repeated sum node
    virtual MX get_repsum(casadi_int m, casadi_int n) const;

    /// Create a vertical concatenation node (vectors only)
    virtual MX get_vertcat(const std::vector<MX>& x) const;

    /// Create a vertical split node (vectors only)
    virtual std::vector<MX> get_vertsplit(const std::vector<casadi_int>& output_offset) const;

    /// Create a diagonal concatenation node
    virtual MX get_diagcat(const std::vector<MX>& x) const;

    /// Create a diagonal split node
    virtual std::vector<MX> get_diagsplit(const std::vector<casadi_int>& offset1,
                                         const std::vector<casadi_int>& offset2) const;

    /// Transpose
    virtual MX get_transpose() const;

    /// Reshape
    virtual MX get_reshape(const Sparsity& sp) const;

    /// Sparsity cast
    virtual MX get_sparsity_cast(const Sparsity& sp) const;

    /** \brief Matrix multiplication and addition

        \identifier{1ru} */
    virtual MX get_mac(const MX& y, const MX& z) const;

    /** \brief Einstein product and addition

        \identifier{1rv} */
    virtual MX get_einstein(const MX& A, const MX& B,
      const std::vector<casadi_int>& dim_c, const std::vector<casadi_int>& dim_a,
      const std::vector<casadi_int>& dim_b,
      const std::vector<casadi_int>& c, const std::vector<casadi_int>& a,
      const std::vector<casadi_int>& b) const;

    /** \brief Bilinear form

        \identifier{1rw} */
    virtual MX get_bilin(const MX& x, const MX& y) const;

    /** \brief Bilinear form

        \identifier{1rx} */
    virtual MX get_rank1(const MX& alpha, const MX& x, const MX& y) const;

    /** \brief Logsumexp

        \identifier{1ry} */
    virtual MX get_logsumexp() const;

    /** \brief Solve a system of linear equations
    *
    *      For system Ax = b:
    *
    *      A->get_solve(b)
    *
        \identifier{1rz} */
    virtual MX get_solve(const MX& r, bool tr, const Linsol& linear_solver) const;

    /** \brief Solve a system of linear equations, upper triangular A
    *
    *      For system Ax = b:
    *
    *      A->get_solve_triu(b)
    *
        \identifier{1s0} */
    virtual MX get_solve_triu(const MX& r, bool tr) const;

    /** \brief Solve a system of linear equations, lower triangular A
    *
    *      For system Ax = b:
    *
    *      A->get_solve_tril(b)
    *
        \identifier{1s1} */
    virtual MX get_solve_tril(const MX& r, bool tr) const;

    /** \brief Solve a system of linear equations, upper triangular A, unity diagonal
    *
    *      For system Ax = b:
    *
    *      A->get_solve_triu(b)
    *
        \identifier{1s2} */
    virtual MX get_solve_triu_unity(const MX& r, bool tr) const;

    /** \brief Solve a system of linear equations, lower triangular A, unity diagnal
    *
    *      For system Ax = b:
    *
    *      A->get_solve_tril(b)
    *
        \identifier{1s3} */
    virtual MX get_solve_tril_unity(const MX& r, bool tr) const;

    /** \brief Get the nonzeros of matrix
    *
    *   a->get_nzref(sp,nz)
    *
    *   returns Matrix(sp,a[nz])

        \identifier{1s4} */
    virtual MX get_nzref(const Sparsity& sp, const std::vector<casadi_int>& nz) const;

    /** \brief Get the nonzeros of matrix, parametrically
    *
        \identifier{1s5} */
    virtual MX get_nz_ref(const MX& nz) const;

    /** \brief Get the nonzeros of matrix, parametrically
    *
        \identifier{1s6} */
    virtual MX get_nz_ref(const MX& inner, const Slice& outer) const;

    /** \brief Get the nonzeros of matrix, parametrically
    *
        \identifier{1s7} */
    virtual MX get_nz_ref(const Slice& inner, const MX& outer) const;

    /** \brief Get the nonzeros of matrix, parametrically
    *
        \identifier{1s8} */
    virtual MX get_nz_ref(const MX& inner, const MX& outer) const;

    /** \brief Assign the nonzeros of a matrix to another matrix
    *
    *   a->get_nzassign(b,nz)
    *   returns b with b[nz]=a

        \identifier{1s9} */
    virtual MX get_nzassign(const MX& y, const std::vector<casadi_int>& nz) const;

    /** \brief Add the nonzeros of a matrix to another matrix
    *
    *   a->get_nzadd(b,nz)
    *   returns b with b[nz]+=a

        \identifier{1sa} */
    virtual MX get_nzadd(const MX& y, const std::vector<casadi_int>& nz) const;

    /** \brief Assign the nonzeros of a matrix to another matrix, parametrically
    *
    *   a->get_nzassign(b,nz)
    *   returns b with b[nz]=a

        \identifier{1sb} */
    virtual MX get_nzassign(const MX& y, const MX& nz) const;

    /** \brief Assign the nonzeros of a matrix to another matrix, parametrically
    *
    *   a->get_nzassign(b,nz)
    *   returns b with b[nz]=a

        \identifier{1sc} */
    virtual MX get_nzassign(const MX& y, const MX& inner, const Slice& outer) const;

    /** \brief Assign the nonzeros of a matrix to another matrix, parametrically
    *
    *   a->get_nzassign(b,nz)
    *   returns b with b[nz]=a

        \identifier{1sd} */
    virtual MX get_nzassign(const MX& y, const Slice& inner, const MX& outer) const;

    /** \brief Assign the nonzeros of a matrix to another matrix, parametrically
    *
    *   a->get_nzassign(b,nz)
    *   returns b with b[nz]=a

        \identifier{1se} */
    virtual MX get_nzassign(const MX& y, const MX& inner, const MX& outer) const;

    /** \brief Add the nonzeros of a matrix to another matrix, parametrically
    *
    *   a->get_nzadd(b,nz)
    *   returns b with b[nz]+=a

        \identifier{1sf} */
    virtual MX get_nzadd(const MX& y, const MX& nz) const;

    /** \brief Add the nonzeros of a matrix to another matrix, parametrically
    *
    *   a->get_nzadd(b,nz)
    *   returns b with b[nz]+=a

        \identifier{1sg} */
    virtual MX get_nzadd(const MX& y, const MX& inner, const Slice& outer) const;

    /** \brief Add the nonzeros of a matrix to another matrix, parametrically
    *
    *   a->get_nzadd(b,nz)
    *   returns b with b[nz]+=a

        \identifier{1sh} */
    virtual MX get_nzadd(const MX& y, const Slice& inner, const MX& outer) const;

    /** \brief Add the nonzeros of a matrix to another matrix, parametrically
    *
    *   a->get_nzadd(b,nz)
    *   returns b with b[nz]+=a

        \identifier{1si} */
    virtual MX get_nzadd(const MX& y, const MX& inner, const MX& outer) const;

    /// Get submatrix reference
    virtual MX get_subref(const Slice& i, const Slice& j) const;

    /// Get submatrix assignment
    virtual MX get_subassign(const MX& y, const Slice& i, const Slice& j) const;

    /// Create set sparse
    virtual MX get_project(const Sparsity& sp) const;

    /// Get a unary operation
    virtual MX get_unary(casadi_int op) const;

    /// Get a binary operation operation
    MX get_binary(casadi_int op, const MX& y) const;

    /// Get a binary operation operation (matrix-matrix)
    virtual MX _get_binary(casadi_int op, const MX& y, bool scX, bool scY) const;

    /// Determinant
    virtual MX get_det() const;

    /// Inverse
    virtual MX get_inv() const;

    /// Inner product
    virtual MX get_dot(const MX& y) const;

    /// Frobenius norm
    virtual MX get_norm_fro() const;

    /// Spectral norm
    virtual MX get_norm_2() const;

    /// Infinity norm
    virtual MX get_norm_inf() const;

    /// 1-norm
    virtual MX get_norm_1() const;

    /// Min
    virtual MX get_mmin() const;

    /// Max
    virtual MX get_mmax() const;

    /// Assertion
    MX get_assert(const MX& y, const std::string& fail_message) const;

    /// Monitor
    MX get_monitor(const std::string& comment) const;

    /// Find
    MX get_find() const;

    /// Find
    MX get_low(const MX& v, const Dict& options) const;

    /// BSpline
    MX get_bspline(const std::vector<double>& knots,
            const std::vector<casadi_int>& offset,
            const std::vector<double>& coeffs,
            const std::vector<casadi_int>& degree,
            casadi_int m,
            const std::vector<casadi_int>& lookup_mode) const;
    /// BSpline
    MX get_bspline(const MX& coeffs, const std::vector<double>& knots,
            const std::vector<casadi_int>& offset,
            const std::vector<casadi_int>& degree,
            casadi_int m,
            const std::vector<casadi_int>& lookup_mode) const;

    /// Convexify
    MX get_convexify(const Dict& opts) const;

    /** Temporary variables to be used in user algorithms like sorting,
        the user is responsible of making sure that use is thread-safe
        The variable is initialized to zero
    */
    mutable casadi_int temp;

    /** \brief  dependencies - functions that have to be evaluated before this one

        \identifier{1sj} */
    std::vector<MX> dep_;

    /** \brief  The sparsity pattern

        \identifier{1sk} */
    Sparsity sparsity_;

    /** \brief Propagate sparsities forward through a copy operation

        \identifier{1sl} */
    static void copy_fwd(const bvec_t* arg, bvec_t* res, casadi_int len);

    /** \brief Propagate sparsities backwards through a copy operation

        \identifier{1sm} */
    static void copy_rev(bvec_t* arg, bvec_t* res, casadi_int len);

    static std::map<casadi_int, MXNode* (*)(DeserializingStream&)> deserialize_map;

  protected:
    /** \brief Deserializing constructor

        \identifier{1sn} */
    explicit MXNode(DeserializingStream& s);
  };

  /// \endcond
} // namespace casadi

#endif // CASADI_MX_NODE_HPP
