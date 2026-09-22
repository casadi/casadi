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
#include "casadi/core/timing.hpp"
#include <ipmc/ipmc.h>

namespace casadi {
  #include "ipmc_runtime.hpp"
}

/** */

/** \pluginsection{Nlpsol,ipmc}

Ipmc is a reverse-communication block-structure exploiting nonlinear interior point method.

Inspired by 
  Vanroye, L., De Schutter, J., & Decré, W. (2024). Efficient Numerical Algorithms for Nonlinear Optimal Control.


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


**/

/// \cond INTERNAL
namespace casadi {

  class IpmcInterface;

  struct CASADI_NLPSOL_IPMC_EXPORT IpmcMemory : public NlpsolMemory {
    // Problem data structure
    casadi_ipmc_data<double> d;

    /* The REWRITTEN problem's z-space, used only when the interface lifted a
       cross-stage slack column (n_lift_>0).  d.nlp points here and d.rewrite.user
       at NlpsolMemory::d_nlp, and the two are the same pointer otherwise. */
    casadi_nlpsol_data<double> d_nlp_lift;

    // The interface this memory belongs to.
    const IpmcInterface* self = nullptr;
  };

  /** \brief \pluginbrief{Nlpsol,ipmc}

      @copydoc Nlpsol_doc
      @copydoc plugin_Nlpsol_ipmc
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

    /** \brief Initalize memory block */
    int init_mem(void* mem) const override;

    /** \brief Free memory block */
    void free_mem(void *mem) const override;

    /// Get all statistics
    Dict get_stats(void* mem) const override;

    /** \brief Trampoline from the solve loop into Nlpsol::callback
     *
     * Installed on the memory object as a plain function pointer, so that the
     * loop that reports progress is the same code the generated C runs (with
     * the pointer null there). */

    /** \brief Set the (persistent) work vectors */
    void set_work(void* mem, const double**& arg, double**& res,
                          casadi_int*& iw, double*& w) const override;

    // Solve the NLP
    int solve(void* mem) const override;

    /// Exact Hessian?
    bool exact_hessian_;

    /// All IPMC options
    Dict opts_;

    /// A documentation string
    static const std::string meta_doc;

    // Options

    /// Data for convexification
    ConvexifyData convexify_data_;

    /// convexify?
    bool convexify_;

    void set_ipmc_prob();
    void set_ipmc_prob(CodeGenerator& g) const;

    /** \brief Generate code for the function body */
    void codegen_body(CodeGenerator& g) const override;

    /** \brief Generate code for the declarations of the C function */
    void codegen_declarations(CodeGenerator& g) const override;

    /** \brief Codegen alloc_mem */
    void codegen_init_mem(CodeGenerator& g) const override;

    /** \brief Codegen free_mem */
    void codegen_free_mem(CodeGenerator& g) const override;

    /** \brief Thread-local memory object type */
    std::string codegen_mem_type() const override { return "struct casadi_ipmc_data"; }

    /** \brief Is thread-local memory object needed? */
    bool codegen_needs_mem() const override { return true; }

    /** \brief Serialize an object without type information */
    void serialize_body(SerializingStream &s) const override;

    /** \brief Deserialize into MX */
    static ProtoFunction* deserialize(DeserializingStream& s) { return new IpmcInterface(s); }

  protected:
    /** \brief Deserializing constructor */
    explicit IpmcInterface(DeserializingStream& s);

  private:
    // Memory structure
    casadi_ipmc_prob<double> p_;

    static Sparsity blocksparsity(casadi_int rows, casadi_int cols,
                                   const std::vector<casadi_ocp_block>& blocks, bool eye=false);
    static void blockptr(std::vector<double *>& vs, std::vector<double>& v,
      const std::vector<casadi_ocp_block>& blocks, bool eye=false);
    Sparsity Isp_, ABsp_, CDsp_, RSQsp_;

    std::vector< casadi_ocp_block > AB_blocks_, CD_blocks_, RSQ_blocks_, I_blocks_;

    std::vector<casadi_int> nxs_;
    std::vector<casadi_int> nus_;
    std::vector<casadi_int> ngs_;
    /* how many TRAILING states of every stage are constant (x_{k+1}=x_k).
       One number for the whole horizon; 0 when the user asked for none. */
    casadi_int nxc_;
    casadi_int N_;

    // ---- native soft constraints (slacks), PHASE 1: stage-local only -----
    // NOTE: nxs_/nus_/ngs_ above SHADOW Nlpsol::nxs_/ngs_ (which are slack
    // counts, not stagewise vectors).  Every base-class access below is
    // qualified as Nlpsol::.
    //
    // The slack -> stage map.  Built once in init() from Nlpsol::slack_S_ and
    // the CD block partition; a column whose rows span several stages is
    // rejected there (Phase 2 / DESIGN_linf_lifting.md lifts those, and only
    // has to add lifted columns to slack_stage_/slack_perm_/slack_offs_ --
    // slack_g_/slack_x_ keep their meaning unchanged).
    bool slacks_;                         // native slack mode is active
    std::vector<casadi_int> slack_g_;     // [ngu_] g row  -> casadi column, -1 hard
    std::vector<casadi_int> slack_x_;     // [nxu_] x index-> casadi column, -1 hard
    std::vector<casadi_int> slack_stage_; // [ns_]  casadi column -> stage
    std::vector<casadi_int> slack_perm_;  // [ns_]  ipmc slack index -> column
    std::vector<casadi_int> slack_ns_;    // [N_+1] columns local to stage k
    std::vector<casadi_int> slack_offs_;  // [N_+2] prefix sums of the per-stage count
    // The penalty, settled in init(): f_s is validated to be a separable
    // quadratic in s whose coefficients do not depend on p, so z = grad:f_s:s
    // at s=0 and Z = diag hess:f_s:s:s are plain numbers.  Dense, length
    // 2*ns_, laid out [lower; upper] like casadi's slack vector.  Nothing
    // evaluates f_s after this -- not the rewrite, not the runtime, not
    // generated C.
    std::vector<double> fs_z_;
    std::vector<double> fs_Z_;

    // ---- PHASE 2: cross-stage (L-infinity) columns, lifted at RUNTIME ----
    // A column whose softened rows span several stages cannot be expressed
    // through ipmc's slack pairs at all -- a pair belongs to ONE stage, and
    // ipmc's soft_idx is a stage-local index.  Such a column is therefore
    // removed from the problem before ipmc ever sees it: init() REWRITES the
    // problem STRUCTURE into an equivalent one that has no cross-stage
    // column left, and hands that to the ordinary un-lifted path.  The
    // ORACLE stays the caller's: the rewritten quantities are produced
    // numerically in the solve loop by the runtime's lift_* routines,
    // driven by the maps below (a symbolic lift_oracle was tried and cost
    // more than it saved: the wrapped call node defeats dead-code
    // elimination, see git history).
    //
    // The rewrite, one helper state component ("entity") per softened SIDE:
    //   * m is inserted at the END of every stage's state block, so the
    //     decision vector becomes [x_0, m, u_0, x_1, m, u_1, ..., x_N, m];
    //   * m_{k+1} - m_k = 0 is appended to the gap-closing rows of every
    //     transition, after the nx[k+1] casadi ones, so the dynamics row
    //     order still matches the state order;
    //   * every softened row becomes an ordinary HARD stage inequality
    //     c + m_l >= lb / c - m_u <= ub, so a two-sided row SPLITS in two
    //     (the two sides need different helper states); a softened SIMPLE
    //     bound becomes a path row of the stage owning the variable, and its
    //     own lbx/ubx is opened to +-inf;
    //   * 0 <= m <= ubs is an ordinary simple bound on the stage-0 copy, and
    //     ubs == 0 (or a side no row uses) makes it a fixed variable, which
    //     the existing equality/inequality split already handles;
    //   * the penalty on the helpers is sum_e (z_e m_e + 1/2 Z_e m_e^2),
    //     added to f, with z_e/Z_e the constants fs_z_/fs_Z_ hold for that
    //     entity's slack slot -- which is exactly f_s(s_subst,p) - f_s(0,p)
    //     with s_subst carrying m_0 in the lifted columns' slots, f_s being
    //     a separable quadratic that does not depend on p.
    // The helpers are trailing-aligned with identity dynamics, so they are
    // declared to ipmc as nxc_ + n_lift_ constant states and cost the Riccati
    // recursion nothing.  Everything above is ORDINARY OCP structure: the
    // runtime, ipmc and the whole eq/ineq machinery see a plain problem.
    //
    // What is left for the interface at solve time is a pure vector scatter
    // between the caller's z-space and the rewritten one, because ubs, lbg
    // and ubg are runtime inputs.  The maps below drive it (and its codegen
    // mirror in the runtime header); they are all empty when n_lift_ == 0,
    // and then nothing at all is different from the un-lifted path.
    casadi_int n_lift_;                  // helper components, 2 per lifted column
    casadi_int nxt_, nat_;               // REWRITTEN nx / ng totals
    std::vector<casadi_int> lift_x_;     // [nx_]  caller x index -> rewritten
    std::vector<char> lift_xfree_;       // [nx_]  its bound moved to a path row
    std::vector<casadi_int> lift_row_;   // [nat_] rewritten row -> the caller's
                                         //        Z-space index (i for a bound,
                                         //        nx_+r for a row), -1 if none
    std::vector<char> lift_kind_;        // [nat_] 0 plain, 1 lower half,
                                         //        2 upper half, 3 helper dyn
    std::vector<casadi_int> lift_ent_;   // [nat_] rewritten row -> entity, -1
    std::vector<casadi_int> lift_mflat_; // [nat_] rewritten row -> flat helper
                                         //        index (stage*ne + entity), -1
    std::vector<casadi_int> lift_m_;     // [(N_+1)*n_lift_] entity e of stage k
    std::vector<casadi_int> lift_slot_;  // [n_lift_] entity -> slack slot
    // The rewritten jacobian/hessian patterns and their per-nonzero source
    // maps (caller nonzero, or a constant; see the runtime's rewrite block).
    // Built by lift_patterns() from the caller patterns and the maps above;
    // empty when n_lift_ == 0.
    Sparsity jacg_lift_sp_, hess_lift_sp_;
    std::vector<casadi_int> lift_a_src_, lift_h_src_;
    void lift_patterns();
    // The REWRITTEN nlp problem sizes, as the ipmc runtime reads them.
    casadi_nlpsol_prob<double> p_nlp_lift_;

    /* The CALLER's partition, kept only so that stats() reports the horizon
       the caller wrote rather than the rewritten one.  Equal to nxs_/ngs_
       whenever nothing was rewritten. */
    std::vector<casadi_int> nxs_u_, ngs_u_;

    // Create nlp_f / nlp_g / nlp_grad_f / nlp_jac_g / nlp_hess_l from a
    // given oracle, and record jacg_sp_ / hesslag_sp_.
    void create_ipmc_functions(const Function& orc);

    // An enum field for the structure detection
    enum StructureDetection {
      STRUCTURE_NONE,
      STRUCTURE_AUTO,
      STRUCTURE_MANUAL
    };
    StructureDetection structure_detection_;

    std::vector<casadi_int> AB_offsets_, CD_offsets_, RSQ_offsets_, I_offsets_;
    bool debug_;
  };

} // namespace casadi
/// \endcond

#endif // CASADI_IPMC_INTERFACE_HPP
