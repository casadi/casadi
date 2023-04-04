#pragma once

#include <alpaqa/export.hpp>
#include <alpaqa/inner/directions/panoc-direction-update.hpp>
#include <alpaqa/inner/inner-solve-options.hpp>
#include <alpaqa/inner/internal/lipschitz.hpp>
#include <alpaqa/inner/internal/panoc-helpers.hpp>
#include <alpaqa/inner/internal/panoc-stop-crit.hpp>
#include <alpaqa/inner/internal/solverstatus.hpp>
#include <alpaqa/problem/type-erased-problem.hpp>
#include <alpaqa/util/atomic-stop-signal.hpp>

#include <chrono>
#include <iostream>
#include <limits>
#include <string>
#include <type_traits>

namespace alpaqa {

/// Tuning parameters for the ZeroFPR algorithm.
template <Config Conf = DefaultConfig>
struct ZeroFPRParams {
    USING_ALPAQA_CONFIG(Conf);

    /// Parameters related to the Lipschitz constant estimate and step size.
    LipschitzEstimateParams<config_t> Lipschitz;
    /// Maximum number of inner ZeroFPR iterations.
    unsigned max_iter = 100;
    /// Maximum duration.
    std::chrono::nanoseconds max_time = std::chrono::minutes(5);
    /// Minimum weight factor between Newton step and projected gradient step.
    real_t min_linesearch_coefficient = real_t(1. / 256);
    /// Ignore the line search condition and always accept the accelerated step.
    /// (For testing purposes only).
    bool force_linesearch = false;
    /// Parameter β used in the line search (see Algorithm 2 in
    /// @cite de_marchi_proximal_2022). @f$ 0 < \beta < 1 @f$
    real_t linesearch_strictness_factor = real_t(0.95);
    /// Minimum Lipschitz constant estimate.
    real_t L_min = real_t(1e-5);
    /// Maximum Lipschitz constant estimate.
    real_t L_max = real_t(1e20);
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
    real_t linesearch_tolerance_factor =
        10 * std::numeric_limits<real_t>::epsilon();

    bool update_direction_in_candidate              = false;
    bool recompute_last_prox_step_after_lbfgs_flush = false;
    bool update_direction_from_prox_step            = false;
};

template <Config Conf = DefaultConfig>
struct ZeroFPRStats {
    USING_ALPAQA_CONFIG(Conf);

    SolverStatus status = SolverStatus::Busy;
    real_t ε            = inf<config_t>;
    std::chrono::nanoseconds elapsed_time{};
    std::chrono::nanoseconds time_progress_callback{};
    unsigned iterations            = 0;
    unsigned linesearch_failures   = 0;
    unsigned linesearch_backtracks = 0;
    unsigned stepsize_backtracks   = 0;
    unsigned lbfgs_failures        = 0;
    unsigned lbfgs_rejected        = 0;
    unsigned τ_1_accepted          = 0;
    unsigned count_τ               = 0;
    real_t sum_τ                   = 0;
    real_t final_γ                 = 0;
    real_t final_ψ                 = 0;
    real_t final_h                 = 0;
    real_t final_φγ                = 0;
};

template <Config Conf = DefaultConfig>
struct ZeroFPRProgressInfo {
    USING_ALPAQA_CONFIG(Conf);

    unsigned k;
    SolverStatus status;
    crvec x;
    crvec p;
    real_t norm_sq_p;
    crvec x̂;
    real_t φγ;
    real_t ψ;
    crvec grad_ψ;
    real_t ψ_hat;
    crvec grad_ψ_hat;
    crvec q;
    real_t L;
    real_t γ;
    real_t τ;
    real_t ε;
    crvec Σ;
    crvec y;
    unsigned outer_iter;
    const TypeErasedProblem<config_t> *problem;
    const ZeroFPRParams<config_t> *params;
};

/// ZeroFPR solver for ALM.
/// @ingroup    grp_InnerSolvers
template <class DirectionT>
class ZeroFPRSolver {
  public:
    USING_ALPAQA_CONFIG_TEMPLATE(DirectionT::config_t);

    using Problem      = TypeErasedProblem<config_t>;
    using Params       = ZeroFPRParams<config_t>;
    using Direction    = DirectionT;
    using Stats        = ZeroFPRStats<config_t>;
    using ProgressInfo = ZeroFPRProgressInfo<config_t>;
    using SolveOptions = InnerSolveOptions<config_t>;

    ZeroFPRSolver(const Params &params)
        requires std::default_initializable<Direction>
        : params(params) {}
    ZeroFPRSolver(const Params &params, Direction &&direction)
        : params(params), direction(std::move(direction)) {}
    ZeroFPRSolver(const Params &params, const Direction &direction)
        : params(params), direction(direction) {}

    Stats operator()(const Problem &problem,   // in
                     const SolveOptions &opts, // in
                     rvec x,                   // inout
                     rvec y,                   // inout
                     crvec Σ,                  // in
                     rvec err_z);              // out

    template <class P>
    Stats operator()(const P &problem, const SolveOptions &opts, rvec x, rvec y,
                     crvec Σ, rvec e) {
        return operator()(Problem::template make<P>(problem), opts, x, y, Σ, e);
    }

    template <class P>
    Stats operator()(const P &problem, const SolveOptions &opts, rvec x) {
        if (problem.get_m() != 0)
            throw std::invalid_argument("Missing arguments y, Σ, e");
        mvec y{nullptr, 0}, Σ{nullptr, 0}, e{nullptr, 0};
        return operator()(problem, opts, x, y, Σ, e);
    }

    /// Specify a callable that is invoked with some intermediate results on
    /// each iteration of the algorithm.
    /// @see @ref ProgressInfo
    ZeroFPRSolver &
    set_progress_callback(std::function<void(const ProgressInfo &)> cb) {
        this->progress_cb = cb;
        return *this;
    }

    std::string get_name() const;

    void stop() { stop_signal.stop(); }

    const Params &get_params() const { return params; }

  private:
    Params params;
    AtomicStopSignal stop_signal;
    std::function<void(const ProgressInfo &)> progress_cb;
    using Helpers = detail::PANOCHelpers<config_t>;

  public:
    Direction direction;
    std::ostream *os = &std::cout;
};

template <class InnerSolverStats>
struct InnerStatsAccumulator;

template <Config Conf>
struct InnerStatsAccumulator<ZeroFPRStats<Conf>> {
    USING_ALPAQA_CONFIG(Conf);

    /// Total elapsed time in the inner solver.
    std::chrono::nanoseconds elapsed_time{};
    /// Total time spent in the user-provided progress callback.
    std::chrono::nanoseconds time_progress_callback{};
    /// Total number of inner ZeroFPR iterations.
    unsigned iterations = 0;
    /// Total number of ZeroFPR line search failures.
    unsigned linesearch_failures = 0;
    /// Total number of ZeroFPR line search backtracking steps.
    unsigned linesearch_backtracks = 0;
    /// Total number of ZeroFPR step size reductions.
    unsigned stepsize_backtracks = 0;
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
    /// The final ZeroFPR step size γ.
    real_t final_γ = 0;
    /// Final value of the smooth cost @f$ \psi(\hat x) @f$.
    real_t final_ψ = 0;
    /// Final value of the nonsmooth cost @f$ h(\hat x) @f$.
    real_t final_h = 0;
    /// Final value of the forward-backward envelope, @f$ \varphi_\gamma(x) @f$
    /// (note that this is in the point @f$ x @f$, not @f$ \hat x @f$).
    real_t final_φγ = 0;
};

template <Config Conf>
InnerStatsAccumulator<ZeroFPRStats<Conf>> &
operator+=(InnerStatsAccumulator<ZeroFPRStats<Conf>> &acc,
           const ZeroFPRStats<Conf> &s) {
    acc.iterations += s.iterations;
    acc.elapsed_time += s.elapsed_time;
    acc.time_progress_callback += s.time_progress_callback;
    acc.linesearch_failures += s.linesearch_failures;
    acc.linesearch_backtracks += s.linesearch_backtracks;
    acc.stepsize_backtracks += s.stepsize_backtracks;
    acc.lbfgs_failures += s.lbfgs_failures;
    acc.lbfgs_rejected += s.lbfgs_rejected;
    acc.τ_1_accepted += s.τ_1_accepted;
    acc.count_τ += s.count_τ;
    acc.sum_τ += s.sum_τ;
    acc.final_γ  = s.final_γ;
    acc.final_ψ  = s.final_ψ;
    acc.final_h  = s.final_h;
    acc.final_φγ = s.final_φγ;
    return acc;
}

ALPAQA_EXPORT_EXTERN_TEMPLATE(struct, ZeroFPRParams, DefaultConfig);
ALPAQA_EXPORT_EXTERN_TEMPLATE(struct, ZeroFPRParams, EigenConfigf);
ALPAQA_EXPORT_EXTERN_TEMPLATE(struct, ZeroFPRParams, EigenConfigd);
ALPAQA_EXPORT_EXTERN_TEMPLATE(struct, ZeroFPRParams, EigenConfigl);
#ifdef ALPAQA_WITH_QUAD_PRECISION
ALPAQA_EXPORT_EXTERN_TEMPLATE(struct, ZeroFPRParams, EigenConfigq);
#endif

ALPAQA_EXPORT_EXTERN_TEMPLATE(struct, ZeroFPRStats, DefaultConfig);
ALPAQA_EXPORT_EXTERN_TEMPLATE(struct, ZeroFPRStats, EigenConfigf);
ALPAQA_EXPORT_EXTERN_TEMPLATE(struct, ZeroFPRStats, EigenConfigd);
ALPAQA_EXPORT_EXTERN_TEMPLATE(struct, ZeroFPRStats, EigenConfigl);
#ifdef ALPAQA_WITH_QUAD_PRECISION
ALPAQA_EXPORT_EXTERN_TEMPLATE(struct, ZeroFPRStats, EigenConfigq);
#endif

ALPAQA_EXPORT_EXTERN_TEMPLATE(struct, ZeroFPRProgressInfo, DefaultConfig);
ALPAQA_EXPORT_EXTERN_TEMPLATE(struct, ZeroFPRProgressInfo, EigenConfigf);
ALPAQA_EXPORT_EXTERN_TEMPLATE(struct, ZeroFPRProgressInfo, EigenConfigd);
ALPAQA_EXPORT_EXTERN_TEMPLATE(struct, ZeroFPRProgressInfo, EigenConfigl);
#ifdef ALPAQA_WITH_QUAD_PRECISION
ALPAQA_EXPORT_EXTERN_TEMPLATE(struct, ZeroFPRProgressInfo, EigenConfigq);
#endif

} // namespace alpaqa
