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


#ifndef CASADI_NLPSOL_IMPL_HPP
#define CASADI_NLPSOL_IMPL_HPP

#include "nlpsol.hpp"
#include "oracle_function.hpp"
#include "plugin_interface.hpp"


/// \cond INTERNAL
namespace casadi {

  /** \brief Integrator memory

      \identifier{1ni} */
  struct CASADI_EXPORT NlpsolMemory : public OracleMemory {
    // Problem data structure
    casadi_nlpsol_data<double> d_nlp;
    // number of iterations
    casadi_int n_iter;
    // Success?
    bool success;
    // Return status
    UnifiedReturnStatus unified_return_status;
    // Slack data; not part of casadi_nlpsol_data, which is unaware of slacks.
    // Pointers into arg/res, followed by scratch.
    const double *s0, *ubs, *lam_s0;
    double *s, *lam_s;
    // Scratch for the S_g*s_l product (size ngs_). Expansion mode only.
    double* slack_w;
    // Scratch holding the user-facing quantities for the iteration callback.
    // Expansion mode only.
    double* slack_cb;
    // Native-slack mode only: writable, never-null buffers of size 2*ns_ each.
    // The slacks are NOT part of casadi_nlpsol_data::z in this mode, so this is
    // where they live for the duration of solve(). slack_s/slack_lam_s are
    // pre-filled with the s0/lam_s0 guesses and copied out to res[NLPSOL_S] /
    // res[NLPSOL_LAM_S] afterwards (which may be null); slack_lbs is all zeros
    // and slack_ubs is ubs with +inf substituted for a missing input.
    double *slack_s, *slack_lam_s, *slack_lbs, *slack_ubs;
  };

  /** \brief NLP solver storage class

      @copydoc Nlpsol_doc
      \author Joel Andersson
      \date 2010-2013

      \identifier{1nj} */
  class CASADI_EXPORT
  Nlpsol : public OracleFunction, public PluginInterface<Nlpsol> {
  public:

    // Memory structure
    casadi_nlpsol_prob<double> p_nlp_;

    /// Number of variables
    casadi_int nx_;

    /// Number of constraints
    casadi_int ng_;

    /// Number of parameters
    casadi_int np_;

    ///@{
    /** \brief Slack ("soft constraint") layer

        Two regimes, selected by slack_native_.

        EXPANSION (slack_native_ == false, the default and what every plugin
        that does not opt in gets): the oracle handed to the plugin is a plain
        NLP -- the slacks have been folded into the decision variables and the
        relaxed bounds into extra constraint rows. nx_ and ng_ are therefore
        the *augmented* sizes, the ones the plugin sees, whereas nxu_/ngu_ are
        the user-facing ones.

        NATIVE (slack_native_ == true): the plugin declared
        Exposed::handles_slacks, so create_oracle left the structure intact and
        the oracle is (x, p, s) -> (f, g, f_s) instead. 'x'/'g' are exactly the
        user's, so nx_ == nxu_ and ng_ == ngu_, and 'f' is the user's objective
        WITHOUT the slack penalty. slack_S_ and the index sets below keep their
        meaning; the plugin is responsible for the slack variables themselves
        (see NlpsolMemory::slack_s and friends) and for writing the TOTAL
        objective f + f_s into casadi_nlpsol_data::objective.

        Layout of the augmented decision vector (expansion mode):
          z = [ x (nxu_) ; s (2*ns_) ; g_aug (ng_) ]
        Layout of g_aug:
          [ g[Gh]                    (ngh_ = ngu_ - ngs_)
            g[Gs] + (S_g*s_l)[Gs]    (ngs_)   lower relaxation
            g[Gs] - (S_g*s_u)[Gs]    (ngs_)   upper relaxation
            x[Xs] + (S_x*s_l)[Xs]    (nxs_)   lower relaxation
            x[Xs] - (S_x*s_u)[Xs] ]  (nxs_)   upper relaxation

        \identifier{2k1} */
    /// Slack incidence matrix [S_g;S_x], (ngu_+nxu_) x ns_. The single source
    /// of truth for the slack layer; everything below is derived from it by
    /// set_slack_prob(). Structural only: its entries are implicitly 1.
    Sparsity slack_S_;
    /// True when the plugin handles the slack layer natively, i.e. when the
    /// oracle is (x,p,s)->(f,g,f_s) rather than the de-sugared plain NLP.
    /// Derived from the oracle's arity, never stored, so it cannot desync
    /// across serialization.
    bool slack_native_;
    /// Number of slacks per side; there are 2*ns_ slack variables
    casadi_int ns_;
    /// Number of user-facing decision variables:
    /// nx_ == nxu_ + 2*ns_ (expansion) / nx_ == nxu_ (native)
    casadi_int nxu_;
    /// Number of user-facing constraints
    casadi_int ngu_;
    /// Indices into user g of the softened constraints (sorted, size ngs_)
    std::vector<casadi_int> slack_target_g_;
    /// Indices into user x of the softened simple bounds (sorted, size nxs_)
    std::vector<casadi_int> slack_target_x_;
    /// Indices into user g of the hard constraints (complement, size ngh_)
    std::vector<casadi_int> slack_hard_g_;
    /// S_g restricted to its softened rows, ngs_ x ns_, and its all-ones values
    Sparsity slack_sp_sg_;
    std::vector<double> slack_sg_ones_;
    /// Shorthands
    casadi_int ngs_, nxs_, ngh_;
    ///@}

    /// callback function, executed at each iteration
    Function fcallback_;

    /// Execute the callback function only after this amount of iterations
    casadi_int callback_step_;

    /// Linear solver and options
    std::string sens_linsol_;
    Dict sens_linsol_options_;

    std::vector<char> detect_simple_bounds_is_simple_;
    Function detect_simple_bounds_parts_;
    std::vector<casadi_int> detect_simple_bounds_target_x_;
    std::vector<casadi_int> detect_simple_bounds_target_g_;

    ///@{
    /** \brief Options

        \identifier{1nk} */
    bool eval_errors_fatal_;
    bool warn_initial_bounds_;
    bool iteration_callback_ignore_errors_;
    bool calc_multipliers_;
    bool calc_lam_x_, calc_lam_p_, calc_f_, calc_g_;
    bool bound_consistency_;
    double min_lam_;
    bool no_nlp_grad_;
    std::vector<bool> discrete_;
    std::vector<bool> equality_;
    ///@}

    // Mixed integer problem?
    bool mi_;

    /// Cache for KKT function
    mutable WeakRef kkt_;

#ifdef CASADI_WITH_THREADSAFE_SYMBOLICS
    /// Mutex for thread safety
    mutable std::mutex kkt_mtx_;
#endif // CASADI_WITH_THREADSAFE_SYMBOLICS

    /** \brief Serialize an object without type information

        \identifier{1nl} */
    void serialize_body(SerializingStream &s) const override;
    /** \brief Serialize type information

        \identifier{1nm} */
    void serialize_type(SerializingStream &s) const override;

    /** \brief Deserialize into MX

        \identifier{1nn} */
    static ProtoFunction* deserialize(DeserializingStream& s);

    /** \brief String used to identify the immediate FunctionInternal subclass

        \identifier{1no} */
    std::string serialize_base_function() const override { return "Nlpsol"; }

    /// Constructor
    Nlpsol(const std::string& name, const Function& oracle);

    /// Destructor
    ~Nlpsol() override = 0;

    ///@{
    /** \brief Number of function inputs and outputs

        \identifier{1np} */
    size_t get_n_in() override { return NLPSOL_NUM_IN;}
    size_t get_n_out() override { return NLPSOL_NUM_OUT;}
    ///@}

    /// @{
    /** \brief Sparsities of function inputs and outputs

        \identifier{1nq} */
    Sparsity get_sparsity_in(casadi_int i) override;
    Sparsity get_sparsity_out(casadi_int i) override;
    /// @}

    ///@{
    /** \brief Names of function input and outputs

        \identifier{1nr} */
    std::string get_name_in(casadi_int i) override { return nlpsol_in(i);}
    std::string get_name_out(casadi_int i) override { return nlpsol_out(i);}
    /// @}

    ///@{
    /** \brief Options

        \identifier{1ns} */
    static const Options options_;
    const Options& get_options() const override { return options_;}
    ///@}

    /** \brief  Print description

        \identifier{1nt} */
    void disp_more(std::ostream& stream) const override;

    /// Initialize
    void init(const Dict& opts) override;

    /** \brief Create memory block

        \identifier{1nu} */
    void* alloc_mem() const override { return new NlpsolMemory();}

    /** \brief Initalize memory block

        \identifier{1nv} */
    int init_mem(void* mem) const override;

    /** \brief Free memory block

        \identifier{1nw} */
    void free_mem(void *mem) const override { delete static_cast<NlpsolMemory*>(mem);}

    /** \brief Check if the inputs correspond to a well-posed problem

        \identifier{1nx} */
    virtual void check_inputs(void* mem) const;

    ///@{
    /** \brief Reference implementation of the slack layer

        slack_expand scatters the user inputs (x0/lbx/ubx/s0/ubs/lbg/ubg and
        the multiplier guesses) into the augmented z/lbz/ubz/lam;
        slack_collect gathers the augmented solution back into
        x/s/g/lam_x/lam_s/lam_g. Both are no-ops when ns_==0.

        \identifier{2k2} */
    void slack_expand(NlpsolMemory* m) const;
    void slack_collect(NlpsolMemory* m, double* x, double* s, double* g,
                       double* lam_x, double* lam_s, double* lam_g) const;
    void codegen_slack_expand(CodeGenerator& g, const std::string& d_nlp) const;
    void codegen_slack_collect(CodeGenerator& g, const std::string& d_nlp) const;
    ///@}

    ///@{
    /** \brief Native-slack mode: hand the slacks to the plugin as-is

        No scatter/gather is needed (x/g are the user's own), only the slack
        block itself: slack_native_enter materialises s0/lam_s0/ubs into the
        NlpsolMemory scratch, slack_native_exit copies s/lam_s back out.
        Both are no-ops unless slack_native_ && ns_>0.

        \identifier{2k3} */
    void slack_native_enter(NlpsolMemory* m) const;
    void slack_native_exit(NlpsolMemory* m) const;
    void codegen_slack_native_enter(CodeGenerator& g) const;
    void codegen_slack_native_exit(CodeGenerator& g) const;
    ///@}

    /** \brief Get default input value

        \identifier{1ny} */
    double get_default_in(casadi_int ind) const override { return nlpsol_default_in(ind);}

    /// Can discrete variables be treated
    virtual bool integer_support() const { return false;}

    /** \brief Set the (persistent) work vectors

        \identifier{1nz} */
    void set_work(void* mem, const double**& arg, double**& res,
                          casadi_int*& iw, double*& w) const override;

    // Evaluate numerically
    int eval(const double** arg, double** res, casadi_int* iw, double* w, void* mem) const final;

    // Solve the NLP
    virtual int solve(void* mem) const = 0;

    /** \brief Generate code for the declarations of the C function

        \identifier{27j} */
    void codegen_declarations(CodeGenerator& g) const override;

    /** \brief Generate code for the function body

        \identifier{1o0} */
    void codegen_body_enter(CodeGenerator& g) const override;

    /** \brief Generate code for the function body

        \identifier{27k} */
    void codegen_body_exit(CodeGenerator& g) const override;

    // Plugins that prefer to keep d_nlp/p_nlp/d_oracle storage on their own
    // per-mem-block data struct (instead of letting codegen_body_enter emit
    // them as per-call function-scope locals) skip codegen_body_enter and
    // call these directly with their own lvalues. The default-prefix call
    // chain through codegen_body_enter / codegen_body_exit is preserved
    // byte-for-byte for plugins that don't opt in.
    //   codegen_setup_constants: emit p_nlp constants + d_nlp <-> p_nlp /
    //     d_oracle wiring. Runs once -- call from codegen_init_mem.
    //   codegen_setup_per_call: emit arg/res-dependent d_nlp wiring +
    //     casadi_nlpsol_init + copy_default's. Runs every call -- call
    //     from codegen_body.
    //   codegen_post_solve: emit nlp_grad post-solve + bound consistency +
    //     final NLPSOL_X/F/G/LAM_* fan-out. Per-call -- call from codegen_body.
    void codegen_setup_constants(CodeGenerator& g,
        const std::string& d_nlp, const std::string& p_nlp,
        const std::string& d_oracle) const;
    void codegen_setup_per_call(CodeGenerator& g, const std::string& d_nlp) const;
    void codegen_post_solve(CodeGenerator& g, const std::string& d_nlp) const;

    /** \brief Do the derivative functions need nondifferentiated outputs?

        \identifier{1o1} */
    bool uses_output() const override {return true;}

    /// Get all statistics
    Dict get_stats(void* mem) const override;

    ///@{
    /** \brief Generate a function that calculates forward mode derivatives

        \identifier{1o2} */
    bool has_forward(casadi_int nfwd) const override { return true;}
    Function get_forward(casadi_int nfwd, const std::string& name,
                         const std::vector<std::string>& inames,
                         const std::vector<std::string>& onames,
                         const Dict& opts) const override;
    ///@}

    ///@{
    /** \brief Generate a function that calculates reverse mode derivatives

        \identifier{1o3} */
    bool has_reverse(casadi_int nadj) const override { return true;}
    Function get_reverse(casadi_int nadj, const std::string& name,
                         const std::vector<std::string>& inames,
                         const std::vector<std::string>& onames,
                         const Dict& opts) const override;
    ///@}

    // Call the callback function
    int callback(NlpsolMemory* m) const;

    // Get KKT function
    Function kkt() const;

    // Make sure primal-dual solution is consistent with bounds
    static void bound_consistency(casadi_int n, double* z, double* lam,
                                  const double* lbz, const double* ubz);

    // Creator function for internal class
    typedef Nlpsol* (*Creator)(const std::string& name, const Function& oracle);

    // Per-plugin capabilities, consulted before the plugin object exists
    struct Exposed{
      /// The plugin implements the slack ("soft constraint") formulation
      /// itself and wants the oracle handed over unexpanded, i.e.
      /// (x,p,s)->(f,g,f_s). Flips the default of the 'expand_slacks' option.
      bool handles_slacks = false;
    };

    /// Collection of solvers
    static std::map<std::string, Plugin> solvers_;

#ifdef CASADI_WITH_THREADSAFE_SYMBOLICS
    static std::mutex mutex_solvers_;
#endif // CASADI_WITH_THREADSAFE_SYMBOLICS

    /// Infix
    static const std::string infix_;

    /// Short name
    static std::string shortname() { return "nlpsol";}

    /** \brief Get type name

        \identifier{1o4} */
    std::string class_name() const override {return "Nlpsol";}

    /** \brief Check if the function is of a particular type

        \identifier{1o5} */
    bool is_a(const std::string& type, bool recursive) const override;

    // Get reduced Hessian
    virtual DM getReducedHessian();

    /// Read options from parameter xml
    virtual void setOptionsFromFile(const std::string & file);

    /// WORKAROUND: Add an element to an std::vector stored in a GenericType:
    template<typename Type> static void append_to_vec(GenericType& t, Type el) {
      std::vector<Type> v = t;
      v.push_back(el);
      t = v;
    }

    /** \brief Convert dictionary to Problem

        When the problem carries a slack part ('s'/'f_s' in \a d and the
        'S' option in \a opts) and \a expand_slacks is true, it is de-sugared
        here into a plain (x, p) -> (f, g) oracle by augmenting the decision
        variables, the objective and the constraints. Nlpsol::set_slack_prob
        rederives the bookkeeping needed to map back from that same 'S'.

        With \a expand_slacks false the slack structure is left intact and
        handed to the plugin as extra oracle IO instead:
        (x, p, s) -> (f, g, f_s), i.e. NL_INPUTS_S / NL_OUTPUTS_S. Only a
        plugin declaring Exposed::handles_slacks may ask for this.

        \identifier{2k4} */
    template<typename XType>
      static Function create_oracle(const std::map<std::string, XType>& d,
                                    const Dict& opts,
                                    bool expand_slacks = true);

  protected:
    /** \brief Deserializing constructor

        \identifier{1o6} */
    explicit Nlpsol(DeserializingStream& s);
  private:
    void set_nlpsol_prob();
    // Derive the whole slack layer (sizes and index sets) from slack_S_,
    // nx_ and ng_. Called from init() and from the deserializing constructor.
    void set_slack_prob();
  };

} // namespace casadi
/// \endcond
#endif // CASADI_NLPSOL_IMPL_HPP
