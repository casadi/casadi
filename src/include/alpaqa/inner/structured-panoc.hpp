#pragma once

#include <alpaqa/accelerators/lbfgs.hpp>
#include <alpaqa/export.hpp>
#include <alpaqa/inner/internal/lipschitz.hpp>
#include <alpaqa/inner/internal/panoc-helpers.hpp>
#include <alpaqa/inner/internal/panoc-stop-crit.hpp>
#include <alpaqa/inner/internal/solverstatus.hpp>
#include <alpaqa/problem/problem.hpp>
#include <alpaqa/util/atomic-stop-signal.hpp>

#include <chrono>
#include <vector>

namespace alpaqa {

/// Tuning parameters for the structured PANOC algorithm.
template <Config Conf = DefaultConfig>
struct StructuredPANOCLBFGSParams {
    USING_ALPAQA_CONFIG(Conf);

    /// Parameters related to the Lipschitz constant estimate and step size.
    LipschitzEstimateParams<config_t> Lipschitz;
    /// Maximum number of inner PANOC iterations.
    unsigned max_iter = 100;
    /// Maximum duration.
    std::chrono::microseconds max_time = std::chrono::minutes(5);
    /// Minimum weight factor between Newton step and projected gradient step.
    real_t τ_min = 1. / 256;
    /// Minimum Lipschitz constant estimate.
    real_t L_min = 1e-5;
    /// Maximum Lipschitz constant estimate.
    real_t L_max = 1e20;
    /// Factor used in update for exponentially weighted nonmonotone line search.
    /// Zero means monotone line search.
    real_t nonmonotone_linesearch = 0;
    /// Check the FPR norm (with γ = 1) to quickly accept steps without line
    /// search.
    real_t fpr_shortcut_accept_factor = 0.999;
    /// If greater than one, allows nonmonotone FPR decrease when accepting the
    /// FPR shortcut steps. Cannot be zero.
    unsigned fpr_shortcut_history = 1;
    /// What stopping criterion to use.
    PANOCStopCrit stop_crit = PANOCStopCrit::ApproxKKT;
    /// Maximum number of iterations without any progress before giving up.
    unsigned max_no_progress = 10;

    /// When to print progress. If set to zero, nothing will be printed.
    /// If set to N != 0, progress is printed every N iterations.
    unsigned print_interval = 0;
    /// The precision of the floating point values printed by the solver.
    int print_precision = std::numeric_limits<real_t>::max_digits10 / 2;

    real_t quadratic_upperbound_tolerance_factor =
        10 * std::numeric_limits<real_t>::epsilon();

    bool update_lipschitz_in_linesearch = true;
    bool alternative_linesearch_cond    = false;

    bool hessian_vec                    = true;
    bool hessian_vec_finite_differences = true;
    bool full_augmented_hessian         = true;

    unsigned hessian_step_size_heuristic = 0;
};

template <Config Conf = DefaultConfig>
struct StructuredPANOCLBFGSStats {
    USING_ALPAQA_CONFIG(Conf);

    SolverStatus status = SolverStatus::Busy;
    real_t ε            = inf<config_t>;
    std::chrono::microseconds elapsed_time;
    unsigned iterations          = 0;
    unsigned linesearch_failures = 0;
    unsigned lbfgs_failures      = 0;
    unsigned lbfgs_rejected      = 0;
    unsigned τ_1_accepted        = 0;
    unsigned count_τ             = 0;
    real_t sum_τ                 = 0;
    unsigned fpr_shortcuts       = 0;
};

template <Config Conf = DefaultConfig>
struct StructuredPANOCLBFGSProgressInfo {
    USING_ALPAQA_CONFIG(Conf);

    unsigned k;
    crvec x;
    crvec p;
    real_t norm_sq_p;
    crvec x̂;
    crvec q;
    crindexvec J;
    real_t φγ;
    real_t ψ;
    crvec grad_ψ;
    real_t ψ_hat;
    crvec grad_ψ_hat;
    real_t L;
    real_t γ;
    real_t τ;
    real_t ε;
    crvec Σ;
    crvec y;
    const ProblemBase<config_t> &problem;
    const StructuredPANOCLBFGSParams<config_t> &params;
};

/// Second order PANOC solver for ALM.
/// @ingroup    grp_InnerSolvers
template <Config Conf = DefaultConfig>
class StructuredPANOCLBFGSSolver {
  public:
    USING_ALPAQA_CONFIG(Conf);

    using Problem      = alpaqa::ProblemBase<config_t>;
    using Params       = StructuredPANOCLBFGSParams<config_t>;
    using Stats        = StructuredPANOCLBFGSStats<config_t>;
    using ProgressInfo = StructuredPANOCLBFGSProgressInfo<config_t>;
    using LBFGS        = alpaqa::LBFGS<config_t>;
    using LBFGSParams  = alpaqa::LBFGSParams<config_t>;

    StructuredPANOCLBFGSSolver(Params params, LBFGSParams lbfgsparams)
        : params(params), lbfgs(lbfgsparams) {}

    Stats operator()(const Problem &problem,        // in
                     crvec Σ,                       // in
                     real_t ε,                      // in
                     bool always_overwrite_results, // in
                     rvec x,                        // inout
                     rvec y,                        // inout
                     rvec err_z);                   // out

    StructuredPANOCLBFGSSolver &
    set_progress_callback(std::function<void(const ProgressInfo &)> cb) {
        this->progress_cb = cb;
        return *this;
    }

    std::string get_name() const { return "StructuredPANOCLBFGSSolver"; }

    void stop() { stop_signal.stop(); }

    const Params &get_params() const { return params; }

  private:
    Params params;
    AtomicStopSignal stop_signal;
    std::function<void(const ProgressInfo &)> progress_cb;
    using Helpers     = detail::PANOCHelpers<config_t>;
    using indexstdvec = std::vector<index_t>;

  private:
    void compute_quasi_newton_step(const Params &params, const Problem &problem,
                                   real_t γₖ, crvec xₖ, crvec y, crvec Σ,
                                   crvec grad_ψₖ, crvec pₖ, rvec qₖ,
                                   indexstdvec &J, rvec HqK, LBFGS &lbfgs,
                                   rvec work_n, rvec work_n2, rvec work_m);

  public:
    LBFGS lbfgs;
};

template <class InnerSolverStats>
struct InnerStatsAccumulator;

template <Config Conf>
struct InnerStatsAccumulator<StructuredPANOCLBFGSStats<Conf>> {
    USING_ALPAQA_CONFIG(Conf);

    /// Total elapsed time in the inner solver.
    std::chrono::microseconds elapsed_time;
    /// Total number of inner PANOC iterations.
    unsigned iterations = 0;
    /// Total number of PANOC line search failures.
    unsigned linesearch_failures = 0;
    /// Total number of times that the L-BFGS direction was not finite.
    unsigned lbfgs_failures = 0;
    /// Total number of times that the L-BFGS update was rejected (i.e. it
    /// could have resulted in a non-positive definite Hessian estimate).
    unsigned lbfgs_rejected = 0;
    /// Total number of times that a line search parameter of @f$ \tau = 1 @f$
    /// was accepted (i.e. no backtracking necessary).
    unsigned τ_1_accepted = 0;
    /// The total number of line searches performed (used for computing the
    /// average value of @f$ \tau @f$).
    unsigned count_τ = 0;
    /// The sum of the line search parameter @f$ \tau @f$ in all iterations
    /// (used for computing the average value of @f$ \tau @f$).
    real_t sum_τ = 0;
    /// The number of shortcuts taken because of sufficient FPR decrease.
    unsigned fpr_shortcuts = 0;
};

template <Config Conf>
InnerStatsAccumulator<StructuredPANOCLBFGSStats<Conf>> &
operator+=(InnerStatsAccumulator<StructuredPANOCLBFGSStats<Conf>> &acc,
           const StructuredPANOCLBFGSStats<Conf> &s) {
    acc.iterations += s.iterations;
    acc.elapsed_time += s.elapsed_time;
    acc.linesearch_failures += s.linesearch_failures;
    acc.lbfgs_failures += s.lbfgs_failures;
    acc.lbfgs_rejected += s.lbfgs_rejected;
    acc.τ_1_accepted += s.τ_1_accepted;
    acc.count_τ += s.count_τ;
    acc.sum_τ += s.sum_τ;
    return acc;
}

// clang-format off
ALPAQA_EXPORT_EXTERN_TEMPLATE(struct, StructuredPANOCLBFGSParams, DefaultConfig);
ALPAQA_EXPORT_EXTERN_TEMPLATE(struct, StructuredPANOCLBFGSParams, EigenConfigf);
ALPAQA_EXPORT_EXTERN_TEMPLATE(struct, StructuredPANOCLBFGSParams, EigenConfigd);
ALPAQA_EXPORT_EXTERN_TEMPLATE(struct, StructuredPANOCLBFGSParams, EigenConfigl);
#ifdef ALPAQA_WITH_QUAD_PRECISION
ALPAQA_EXPORT_EXTERN_TEMPLATE(struct, StructuredPANOCLBFGSParams, EigenConfigq);
#endif

ALPAQA_EXPORT_EXTERN_TEMPLATE(struct, StructuredPANOCLBFGSStats, DefaultConfig);
ALPAQA_EXPORT_EXTERN_TEMPLATE(struct, StructuredPANOCLBFGSStats, EigenConfigf);
ALPAQA_EXPORT_EXTERN_TEMPLATE(struct, StructuredPANOCLBFGSStats, EigenConfigd);
ALPAQA_EXPORT_EXTERN_TEMPLATE(struct, StructuredPANOCLBFGSStats, EigenConfigl);
#ifdef ALPAQA_WITH_QUAD_PRECISION
ALPAQA_EXPORT_EXTERN_TEMPLATE(struct, StructuredPANOCLBFGSStats, EigenConfigq);
#endif

ALPAQA_EXPORT_EXTERN_TEMPLATE(struct, StructuredPANOCLBFGSProgressInfo, DefaultConfig);
ALPAQA_EXPORT_EXTERN_TEMPLATE(struct, StructuredPANOCLBFGSProgressInfo, EigenConfigf);
ALPAQA_EXPORT_EXTERN_TEMPLATE(struct, StructuredPANOCLBFGSProgressInfo, EigenConfigd);
ALPAQA_EXPORT_EXTERN_TEMPLATE(struct, StructuredPANOCLBFGSProgressInfo, EigenConfigl);
#ifdef ALPAQA_WITH_QUAD_PRECISION
ALPAQA_EXPORT_EXTERN_TEMPLATE(struct, StructuredPANOCLBFGSProgressInfo, EigenConfigq);
#endif

ALPAQA_EXPORT_EXTERN_TEMPLATE(class, StructuredPANOCLBFGSSolver, DefaultConfig);
ALPAQA_EXPORT_EXTERN_TEMPLATE(class, StructuredPANOCLBFGSSolver, EigenConfigf);
ALPAQA_EXPORT_EXTERN_TEMPLATE(class, StructuredPANOCLBFGSSolver, EigenConfigd);
ALPAQA_EXPORT_EXTERN_TEMPLATE(class, StructuredPANOCLBFGSSolver, EigenConfigl);
#ifdef ALPAQA_WITH_QUAD_PRECISION
ALPAQA_EXPORT_EXTERN_TEMPLATE(class, StructuredPANOCLBFGSSolver, EigenConfigq);
#endif
// clang-format on

} // namespace alpaqa
