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


#ifndef CASADI_IPMC_INTERFACE_HPP
#define CASADI_IPMC_INTERFACE_HPP

#include <casadi/interfaces/ipmc/casadi_nlpsol_ipmc_export.h>
#include "casadi/core/nlpsol_impl.hpp"
#include "casadi/core/im.hpp"
#include <ipmc/ipmc.h>

namespace casadi {
  #include "ipmc_runtime.hpp"
}

/** \pluginsection{Nlpsol,ipmc}

Ipmc is a reverse-communication block-structure exploiting nonlinear interior point method
written in C99.

Inspired by
  Vanroye, L., De Schutter, J., & Decré, W. (2024).
  Efficient Numerical Algorithms for Nonlinear Optimal Control.


With structure_detection = 'none' (default),
it will behave as a general-purpose dense nonlinear program solver.

With structure_detection = 'manual', you can specify a block structure.

Let's say you perform multiply shooting with a system

x_k+1 = A_k x_k + B_k u_k


Suppose your constraint Jacobian looks like:

     nx0  nu0  nx1  nu1  nx2  nu2
     -----------------------------
nx1  |A0  B0   I0
ng1  |C0  D0
nx2  |         A1   B1   I1
ng2  |         C1   D1
ng3  |                   C2   D2

with n* capturing the number of states, inputs, and constraints in each block.

You can then specify this structure with:

N = 2
nx = [nx0 ,nx1, nx2]
nu = [nu0, nu1, nu2]
ng = [ng1, ng2, ng3]

With structure_detection = 'auto', the block-defining parameters
nx, nu, ng, and N are automatically detected from the sparsity pattern.

*/

/// \cond INTERNAL
namespace casadi {

  class IpmcInterface;

  struct CASADI_NLPSOL_IPMC_EXPORT IpmcMemory : public NlpsolMemory {
    // Runtime data
    casadi_ipmc_data<double> d;

    // Sizes/work of the lifted problem; d.nlp points here iff n_lift_>0
    casadi_nlpsol_data<double> d_nlp_lift;

    // The block d.solver is built in, ipmc_memsize_ bytes
    std::vector<double> ipmc_block;
  };

  /** \brief \pluginbrief{Nlpsol,ipmc}

      @copydoc Nlpsol_doc
      @copydoc plugin_Nlpsol_ipmc

      init() fixes the structure of the problem ipmc solves and compiles the caller's NLP
      into the constant tables of casadi_ipmc_prob. The runtime makes no structural
      decision: it hands over numbers and moves values through the tables, identically in
      solve() and in generated C.

      Stage partition. N_, nxs_, nus_ and ngs_ come from the options ('manual'), from the
      Jacobian staircase ('auto', detect_structure) or are a single stage ('none').
      AB_blocks_ and CD_blocks_ locate the dynamics and path blocks in the Jacobian.

      Native slacks. Every column j of S is one slack variable s_j >= 0, priced
      z_j s_j + 1/2 Z_j s_j^2 and bounded by ubs_j. f_s must therefore be a separable
      quadratic; settle_slack_penalty verifies that and registers nlp_fs, which evaluates
      z and Z at each solve's p. slack_lo_/slack_up_ say which column relaxes which side
      of which row; each side takes at most one column (slack_maps). A column whose rows lie in one stage becomes an
      ipmc slack of that stage (soft_lo/soft_up), ordered stage after stage in slack_perm_.

      Lift. A column whose rows span several stages, an L-infinity budget, has no
      stage-local representation. lift() rewrites the problem: the column becomes one helper
      state m, appended to every stage with dynamics m_{k+1} = m_k, bounded 0 <= m_0 <= ubs
      and priced at stage 0; ubs = 0 makes ipmc pin it, hardening the sides it relaxes. A
      row with a lifted side splits into a lower half, carrying the lower bound, and an
      upper half, carrying the upper bound; a lifted lower side enters its half as c + m,
      a lifted upper side as c - m, and a side that is not lifted keeps its stage-local
      slack or stays hard. A lifted simple bound becomes such a pair of path rows and frees
      its variable. The helpers are trailing constant states, so ipmc sees nxc_ + n_lift_
      of them, the last n_lift_ declared slack helpers. The rewrite is linear with +-1
      coefficients, held as constant sparse maps (casadi_ipmc_rewrite_prob) and rebuilt from
      S and the stage partition, never stored. The oracle stays the caller's: at solve time
      the runtime maps bounds in, oracle outputs over and the solution back. Lifting the
      oracle symbolically costs more than it saves, since the wrapping MX call node defeats
      dead-code elimination and every single-output function would evaluate the whole
      oracle.

      ipmc's rows (build_rows). Every path row and every variable of the handed-over
      problem is an ipmc inequality row, per stage the path rows first and then the
      variables as simple bounds. ipmc holds a row without a finite bound inert and a hard
      row with equal bounds as an equality, deciding both anew at every solve, so 'equality'
      only serves structure detection. Two halves of a split row that both carry a helper
      are declared twins. With the slack and helper declarations this is the IpmcProblem the
      solver is built from, once per memory object, inside memory casadi owns:
      IpmcMemory::ipmc_block, or a static array in generated C. At every solve
      casadi_ipmc_hand_over checks the bounds of the gap-closing rows and hands over bounds,
      x0, penalty, ubs and s0.

      Pack tables (build_pack_tables). For every entry of ipmc's stage blocks (BAt,
      Gt_ineq, RSQ, RSQ_slack) a table states which nonzero of the caller's Jacobian or
      Hessian it reads, obtained by slicing integer matrices whose entries are nonzero
      codes. Simple bounds have no stored column. Structure that is promised but not handed
      over, the gap-closing identity and the constant states, becomes runtime checks.

      Solve loop. ipmc requests an evaluation; the interface calls the caller's oracle,
      lifts the result and packs it into ipmc's layout. The loop exists twice, in solve()
      and as emitted C in codegen_body().
  */
  class CASADI_NLPSOL_IPMC_EXPORT IpmcInterface : public Nlpsol {
  public:
    Sparsity jacg_sp_;
    Sparsity hesslag_sp_;

    explicit IpmcInterface(const std::string& name, const Function& nlp);
    ~IpmcInterface() override;

    // Get name of the plugin
    const char* plugin_name() const override { return "ipmc";}

    // Get name of the class
    std::string class_name() const override { return "IpmcInterface";}

    /** \brief  Create a new NLP Solver */
    static Nlpsol* creator(const std::string& name, const Function& nlp) {
      return new IpmcInterface(name, nlp);
    }

    ///@{
    /** \brief Options */
    static const Options options_;
    const Options& get_options() const override { return options_;}
    ///@}

    // Initialize the solver
    void init(const Dict& opts) override;

    /** \brief Create memory block */
    void* alloc_mem() const override { return new IpmcMemory();}

    /** \brief Initialize memory block */
    int init_mem(void* mem) const override;

    /** \brief Free memory block */
    void free_mem(void* mem) const override { delete static_cast<IpmcMemory*>(mem);}

    /// Get all statistics
    Dict get_stats(void* mem) const override;

    /** \brief Set the (persistent) work vectors */
    void set_work(void* mem, const double**& arg, double**& res,
                  casadi_int*& iw, double*& w) const override;

    // Serve ipmc's requests; codegen_body emits the same loop
    int solve(void* mem) const override;

    /// Exact Hessian?
    bool exact_hessian_;

    /// All IPMC options
    Dict opts_;

    /// A documentation string
    static const std::string meta_doc;

    void set_ipmc_prob();
    void set_ipmc_prob(CodeGenerator& g) const;
    // The IpmcProblem the solver is built from, in generated C
    std::string codegen_desc(CodeGenerator& g) const;
    // The options of the 'ipmc' dict onto a fresh solver
    void push_options(IpmcSolver* solver) const;
    // The constants of casadi_ipmc_prob, for both set_ipmc_prob
    template<typename V>
    void prob_fields(casadi_ipmc_prob<double>& p, V& v) const;

    /** \brief Generate code for the function body */
    void codegen_body(CodeGenerator& g) const override;

    /** \brief Generate code for the declarations of the C function */
    void codegen_declarations(CodeGenerator& g) const override;

    /** \brief Codegen alloc_mem */
    void codegen_init_mem(CodeGenerator& g) const override;

    /** \brief Thread-local memory object type */
    std::string codegen_mem_type() const override { return "struct casadi_ipmc_data"; }

    /** \brief Is thread-local memory object needed? */
    bool codegen_needs_mem() const override { return true; }

    /** \brief Serialize an object without type information */
    void serialize_body(SerializingStream& s) const override;

    /** \brief Deserialize into MX */
    static ProtoFunction* deserialize(DeserializingStream& s) { return new IpmcInterface(s); }

  protected:
    /** \brief Deserializing constructor */
    explicit IpmcInterface(DeserializingStream& s);

  private:
    // Constant tables read by the runtime
    casadi_ipmc_prob<double> p_;
    // Sizes of the lifted problem; p_.nlp points here iff n_lift_>0
    casadi_nlpsol_prob<double> p_nlp_lift_;

    // The caller's stage partition, k=0..N_; nxs_/ngs_ shadow Nlpsol's slack counts
    casadi_int N_;
    std::vector<casadi_int> nxs_;  // [N_+1] states
    std::vector<casadi_int> nus_;  // [N_+1] controls
    std::vector<casadi_int> ngs_;  // [N_+1] path rows
    casadi_int nxc_;               // trailing constant states per stage, option 'nxc'
    // The handed-over partition: nxs_ plus the n_lift_ helpers, ngs_ with split rows
    std::vector<casadi_int> nxh_, ngh_;
    // Dynamics [A B] and path [C D] blocks of the handed-over Jacobian
    std::vector<casadi_ocp_block> AB_blocks_, CD_blocks_;

    static Sparsity blocksparsity(casadi_int rows, casadi_int cols,
                                   const std::vector<casadi_ocp_block>& blocks, bool eye=false);
    Sparsity identity_sparsity() const;

    // Native slacks
    bool slacks_;                         // slack_native_ and ns>0
    // Derived by slack_maps() from slack_lo_/slack_up_ and the caller's partition
    // [ng_] / [nx_] column relaxing the lower / upper side of a row or variable, -1 if hard
    std::vector<casadi_int> slack_g_lo_, slack_g_up_, slack_x_lo_, slack_x_up_;
    std::vector<casadi_int> slack_perm_;  // [n_soft] column of each stage-local ipmc slack
    std::vector<casadi_int> slack_idx_;   // [ns] ipmc slack of a column, -1 if lifted
    std::vector<casadi_int> slack_ns_;    // [N_+1] stage-local slacks per stage
    std::vector<casadi_int> lift_col_;    // [n_lift_] column of each helper
    std::vector<casadi_int> lift_ent_;    // [ns] helper of a column, -1 if stage-local
    void slack_maps();

    // Lift of cross-stage slack columns
    casadi_int n_lift_;                   // helper states, one per lifted column
    casadi_int nxt_, nat_;                // nx, ng of the handed-over problem
    // Per handed-over z entry: the caller's z entry it carries (-1: none, a helper or the
    // helper dynamics), which of its bounds (LIFT_BOTH, _LO, _UP, _NONE) and its helper
    // term (code 4*e + neg, -1 for none)
    std::vector<casadi_int> lift_src_, lift_side_, lift_hlp_;
    enum { LIFT_BOTH, LIFT_LO, LIFT_UP, LIFT_NONE };
    // Runtime maps, see casadi_ipmc_rewrite_prob
    DM lift_Px_, lift_H_, lift_Pg_, lift_C_, lift_B_, lift_Blo_, lift_Bup_, lift_Gz_;
    std::vector<double> lift_lbz0_, lift_ubz0_;
    std::vector<casadi_int> lift_m0_;
    void lift();

    // ipmc's rows
    std::vector<casadi_int> ineq_z_;      // z-space index of every row, stage after stage
    std::vector<casadi_int> ineq_lo_, ineq_up_;  // [n_ineq] ipmc slack of a side, -1 if none
    // The IpmcProblem, built once per memory object, and its arrays
    std::vector<ipmc_int> d_nu_, d_nx_, d_ng_, d_ng_ineq_, d_ns_, d_soft_lo_, d_soft_up_;
    std::vector<ipmc_int> d_slack_helper_, d_nb_, d_idxb_;
    IpmcProblem desc_;
    size_t ipmc_memsize_;                 // bytes of the block the solver lives in
    void build_rows();
    // The caller's entry, x[i] or g[i], that z-space entry z of the handed-over problem carries
    std::string caller_entry(casadi_int z) const;

    // Pack tables
    // Jacobian / Hessian of the handed-over problem with nonzero codes as entries
    IM jac_codes() const;
    IM hess_codes() const;
    void build_pack_tables();
    Sparsity bat_sp_, rsq_sp_, rsqs_sp_, gi_sp_;
    std::vector<casadi_int> bat_code_, bat_blk_, bat_col_;
    std::vector<casadi_int> rsq_code_, rsq_blk_, rsq_col_;
    std::vector<casadi_int> rsqs_code_;
    std::vector<casadi_int> gi_code_, gi_blk_, gi_col_;
    std::vector<casadi_int> ichk_, cchk_;

    // The phases of init(), in order
    void settle_slack_penalty();
    void detect_structure(std::set<casadi_int>& errors);
    // Everything that follows from the caller's partition and S; also run on deserialization
    void build();

    // Create nlp_f, nlp_g, nlp_grad_f, nlp_jac_g, nlp_hess_l; set jacg_sp_, hesslag_sp_
    void create_ipmc_functions();

    // Who found the stage partition, the 'structure_detection' option
    enum StructureDetection {
      STRUCTURE_NONE,
      STRUCTURE_AUTO,
      STRUCTURE_MANUAL
    };
  };

} // namespace casadi
/// \endcond

#endif // CASADI_IPMC_INTERFACE_HPP
