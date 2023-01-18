#pragma once

#include <alpaqa/export.hpp>
#include <alpaqa/inner/directions/panoc-direction-update.hpp>
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

/// Tuning parameters for the PANOC algorithm.
template <Config Conf = DefaultConfig>
struct PANOCParams {
    USING_ALPAQA_CONFIG(Conf);

    /// Parameters related to the Lipschitz constant estimate and step size.
    LipschitzEstimateParams<config_t> Lipschitz;
    /// Maximum number of inner PANOC iterations.
    unsigned max_iter = 100;
    /// Maximum duration.
    std::chrono::microseconds max_time = std::chrono::minutes(5);
    /// Minimum weight factor between Newton step and projected gradient step.
    real_t τ_min = 1. / 256;
    /// Parameter β used in the line search (see Algorithm 2 in
    /// @cite de_marchi_proximal_2022). @f$ 0 < \beta < 1 @f$
    real_t β = 0.95;
    /// Minimum Lipschitz constant estimate.
    real_t L_min = 1e-5;
    /// Maximum Lipschitz constant estimate.
    real_t L_max = 1e20;
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
};

template <Config Conf = DefaultConfig>
struct PANOCStats {
    USING_ALPAQA_CONFIG(Conf);

    SolverStatus status = SolverStatus::Busy;
    real_t ε            = inf<config_t>;
    std::chrono::nanoseconds elapsed_time{};
    unsigned iterations          = 0;
    unsigned linesearch_failures = 0;
    unsigned lbfgs_failures      = 0;
    unsigned lbfgs_rejected      = 0;
    unsigned τ_1_accepted        = 0;
    unsigned count_τ             = 0;
    real_t sum_τ                 = 0;
    real_t final_γ               = 0;
};

template <Config Conf = DefaultConfig>
struct PANOCProgressInfo {
    USING_ALPAQA_CONFIG(Conf);

    unsigned k;
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
    const TypeErasedProblem<config_t> &problem;
    const PANOCParams<config_t> &params;
};

/// PANOC solver for ALM.
/// @ingroup    grp_InnerSolvers
template <class DirectionT>
class PANOCSolver {
  public:
    USING_ALPAQA_CONFIG_TEMPLATE(DirectionT::config_t);

    using Problem      = TypeErasedProblem<config_t>;
    using Params       = PANOCParams<config_t>;
    using Direction    = DirectionT;
    using Stats        = PANOCStats<config_t>;
    using ProgressInfo = PANOCProgressInfo<config_t>;

    PANOCSolver(const Params &params)
        requires std::default_initializable<Direction>
        : params(params) {}
    PANOCSolver(const Params &params, Direction &&direction)
        : params(params), direction(std::move(direction)) {}
    PANOCSolver(const Params &params, const Direction &direction)
        : params(params), direction(direction) {}

    Stats operator()(const Problem &problem,        // in
                     crvec Σ,                       // in
                     real_t ε,                      // in
                     bool always_overwrite_results, // in
                     rvec x,                        // inout
                     rvec y,                        // inout
                     rvec err_z);                   // out

    /// Specify a callable that is invoked with some intermediate results on
    /// each iteration of the algorithm.
    /// @see @ref ProgressInfo
    PANOCSolver &
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
struct InnerStatsAccumulator<PANOCStats<Conf>> {
    USING_ALPAQA_CONFIG(Conf);

    /// Total elapsed time in the inner solver.
    std::chrono::nanoseconds elapsed_time{};
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
    /// The final PANOC step size γ.
    real_t final_γ = 0;
};

template <Config Conf>
InnerStatsAccumulator<PANOCStats<Conf>> &
operator+=(InnerStatsAccumulator<PANOCStats<Conf>> &acc,
           const PANOCStats<Conf> &s) {
    acc.iterations += s.iterations;
    acc.elapsed_time += s.elapsed_time;
    acc.linesearch_failures += s.linesearch_failures;
    acc.lbfgs_failures += s.lbfgs_failures;
    acc.lbfgs_rejected += s.lbfgs_rejected;
    acc.τ_1_accepted += s.τ_1_accepted;
    acc.count_τ += s.count_τ;
    acc.sum_τ += s.sum_τ;
    acc.final_γ = s.final_γ;
    return acc;
}

ALPAQA_EXPORT_EXTERN_TEMPLATE(struct, PANOCParams, DefaultConfig);
ALPAQA_EXPORT_EXTERN_TEMPLATE(struct, PANOCParams, EigenConfigf);
ALPAQA_EXPORT_EXTERN_TEMPLATE(struct, PANOCParams, EigenConfigd);
ALPAQA_EXPORT_EXTERN_TEMPLATE(struct, PANOCParams, EigenConfigl);
#ifdef ALPAQA_WITH_QUAD_PRECISION
ALPAQA_EXPORT_EXTERN_TEMPLATE(struct, PANOCParams, EigenConfigq);
#endif

ALPAQA_EXPORT_EXTERN_TEMPLATE(struct, PANOCStats, DefaultConfig);
ALPAQA_EXPORT_EXTERN_TEMPLATE(struct, PANOCStats, EigenConfigf);
ALPAQA_EXPORT_EXTERN_TEMPLATE(struct, PANOCStats, EigenConfigd);
ALPAQA_EXPORT_EXTERN_TEMPLATE(struct, PANOCStats, EigenConfigl);
#ifdef ALPAQA_WITH_QUAD_PRECISION
ALPAQA_EXPORT_EXTERN_TEMPLATE(struct, PANOCStats, EigenConfigq);
#endif

ALPAQA_EXPORT_EXTERN_TEMPLATE(struct, PANOCProgressInfo, DefaultConfig);
ALPAQA_EXPORT_EXTERN_TEMPLATE(struct, PANOCProgressInfo, EigenConfigf);
ALPAQA_EXPORT_EXTERN_TEMPLATE(struct, PANOCProgressInfo, EigenConfigd);
ALPAQA_EXPORT_EXTERN_TEMPLATE(struct, PANOCProgressInfo, EigenConfigl);
#ifdef ALPAQA_WITH_QUAD_PRECISION
ALPAQA_EXPORT_EXTERN_TEMPLATE(struct, PANOCProgressInfo, EigenConfigq);
#endif

} // namespace alpaqa
