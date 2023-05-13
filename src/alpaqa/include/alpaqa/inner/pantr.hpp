#pragma once

#include <alpaqa/export.hpp>
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

/// Tuning parameters for the PANTR algorithm.
template <Config Conf = DefaultConfig>
struct PANTRParams {
    USING_ALPAQA_CONFIG(Conf);

    /// Parameters related to the Lipschitz constant estimate and step size.
    LipschitzEstimateParams<config_t> Lipschitz;
    /// Maximum number of inner PANTR iterations.
    unsigned max_iter = 100;
    /// Maximum duration.
    std::chrono::nanoseconds max_time = std::chrono::minutes(5);
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
    real_t TR_tolerance_factor = 10 * std::numeric_limits<real_t>::epsilon();

    /// Minimal TR ratio to be accepted (successful).
    real_t ratio_threshold_acceptable = real_t(0.2);
    /// Minimal TR ratio to increase radius (very successful).
    real_t ratio_threshold_good = real_t(0.8);

    /// TR radius decrease coefficient when unsuccessful.
    real_t radius_factor_rejected = real_t(0.35);
    /// TR radius decrease coefficient when successful.
    real_t radius_factor_acceptable = real_t(0.999);
    /// TR radius increase coefficient when very successful.
    real_t radius_factor_good = real_t(2.5);

    /// Initial trust radius.
    real_t initial_radius = NaN<config_t>;
    /// Minimum trust radius.
    real_t min_radius = 100 * std::numeric_limits<real_t>::epsilon();

    /// Check the quadratic upperbound and update γ before computing the
    /// reduction of the TR step.
    bool compute_ratio_using_new_stepsize = false;

    bool update_direction_on_prox_step                  = true;
    bool recompute_last_prox_step_after_direction_reset = false;

    /// Don't compute accelerated steps, fall back to forward-backward splitting.
    /// For testing purposes.
    bool disable_acceleration = false;
    /// Compute the trust-region ratio using an approximation of the quadratic
    /// model of the FBE, rather than the quadratic model of the subproblem.
    /// Specifically, when set to false, the quadratic model used is
    /// @f[ q(d) = \tfrac12 \inprod{\mathcal R_\gamma(\hat x) d}{d} +
    /// \inprod{R_\gamma(\hat x)}{d}. @f]
    /// When set to true, the quadratic model used is
    /// @f[ q_\mathrm{approx}(d) = \inv{(1-\alpha)} q(d), @f]
    /// where @f$ \alpha = @f$ @ref LipschitzParams::Lγ_factor.
    /// This is an approximation of the quadratic model of the FBE,
    /// @f[ q_{\varphi_\gamma}(d) = \tfrac12 \inprod{\mathcal Q_\gamma(\hat x)
    /// \mathcal R_\gamma(\hat x) d}{d} + \inprod{\mathcal Q_\gamma(\hat x)
    /// R_\gamma(\hat x)}{d}, @f]
    /// with @f$ \mathcal Q_\gamma(x) = \id - \gamma \nabla^2 \psi(x) @f$.
    bool ratio_approx_fbe_quadratic_model = true;
};

template <Config Conf = DefaultConfig>
struct PANTRStats {
    USING_ALPAQA_CONFIG(Conf);

    SolverStatus status = SolverStatus::Busy;
    real_t ε            = inf<config_t>;
    std::chrono::nanoseconds elapsed_time{};
    std::chrono::nanoseconds time_progress_callback{};
    unsigned iterations                = 0;
    unsigned accelerated_step_rejected = 0;
    unsigned stepsize_backtracks       = 0;
    unsigned direction_failures        = 0;
    unsigned direction_update_rejected = 0;
    real_t final_γ                     = 0;
    real_t final_ψ                     = 0;
    real_t final_h                     = 0;
    real_t final_φγ                    = 0;
};

template <Config Conf = DefaultConfig>
struct PANTRProgressInfo {
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
    real_t Δ;
    real_t ρ;
    real_t ε;
    crvec Σ;
    crvec y;
    unsigned outer_iter;
    const TypeErasedProblem<config_t> *problem;
    const PANTRParams<config_t> *params;
};

/// PANTR solver for ALM.
/// @ingroup    grp_InnerSolvers
template <class DirectionT>
class PANTRSolver {
  public:
    USING_ALPAQA_CONFIG_TEMPLATE(DirectionT::config_t);

    using Problem      = TypeErasedProblem<config_t>;
    using Params       = PANTRParams<config_t>;
    using Direction    = DirectionT;
    using Stats        = PANTRStats<config_t>;
    using ProgressInfo = PANTRProgressInfo<config_t>;
    using SolveOptions = InnerSolveOptions<config_t>;

    PANTRSolver(const Params &params)
        requires std::default_initializable<Direction>
        : params(params) {}
    PANTRSolver(const Params &params, Direction &&direction)
        : params(params), direction(std::move(direction)) {}
    PANTRSolver(const Params &params, const Direction &direction)
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
    PANTRSolver &
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
struct InnerStatsAccumulator<PANTRStats<Conf>> {
    USING_ALPAQA_CONFIG(Conf);

    /// Total elapsed time in the inner solver.
    std::chrono::nanoseconds elapsed_time{};
    /// Total time spent in the user-provided progress callback.
    std::chrono::nanoseconds time_progress_callback{};
    /// Total number of inner PANTR iterations.
    unsigned iterations = 0;
    /// Total number of PANTR rejected accelerated steps.
    unsigned accelerated_step_rejected = 0;
    /// Total number of FB step size reductions.
    unsigned stepsize_backtracks = 0;
    /// Total number of times that the accelerated direction was not finite.
    unsigned direction_failures = 0;
    /// Total number of times that the direction update was rejected (e.g. it
    /// could have resulted in a non-positive definite Hessian estimate).
    unsigned direction_update_rejected = 0;
    /// The final FB step size γ.
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
InnerStatsAccumulator<PANTRStats<Conf>> &
operator+=(InnerStatsAccumulator<PANTRStats<Conf>> &acc,
           const PANTRStats<Conf> &s) {
    acc.elapsed_time += s.elapsed_time;
    acc.time_progress_callback += s.time_progress_callback;
    acc.iterations += s.iterations;
    acc.accelerated_step_rejected += s.accelerated_step_rejected;
    acc.stepsize_backtracks += s.stepsize_backtracks;
    acc.direction_failures += s.direction_failures;
    acc.direction_update_rejected += s.direction_update_rejected;
    acc.final_γ  = s.final_γ;
    acc.final_ψ  = s.final_ψ;
    acc.final_h  = s.final_h;
    acc.final_φγ = s.final_φγ;
    return acc;
}

ALPAQA_EXPORT_EXTERN_TEMPLATE(struct, PANTRParams, DefaultConfig);
ALPAQA_EXPORT_EXTERN_TEMPLATE(struct, PANTRParams, EigenConfigf);
ALPAQA_EXPORT_EXTERN_TEMPLATE(struct, PANTRParams, EigenConfigd);
ALPAQA_EXPORT_EXTERN_TEMPLATE(struct, PANTRParams, EigenConfigl);
#ifdef ALPAQA_WITH_QUAD_PRECISION
ALPAQA_EXPORT_EXTERN_TEMPLATE(struct, PANTRParams, EigenConfigq);
#endif

ALPAQA_EXPORT_EXTERN_TEMPLATE(struct, PANTRStats, DefaultConfig);
ALPAQA_EXPORT_EXTERN_TEMPLATE(struct, PANTRStats, EigenConfigf);
ALPAQA_EXPORT_EXTERN_TEMPLATE(struct, PANTRStats, EigenConfigd);
ALPAQA_EXPORT_EXTERN_TEMPLATE(struct, PANTRStats, EigenConfigl);
#ifdef ALPAQA_WITH_QUAD_PRECISION
ALPAQA_EXPORT_EXTERN_TEMPLATE(struct, PANTRStats, EigenConfigq);
#endif

ALPAQA_EXPORT_EXTERN_TEMPLATE(struct, PANTRProgressInfo, DefaultConfig);
ALPAQA_EXPORT_EXTERN_TEMPLATE(struct, PANTRProgressInfo, EigenConfigf);
ALPAQA_EXPORT_EXTERN_TEMPLATE(struct, PANTRProgressInfo, EigenConfigd);
ALPAQA_EXPORT_EXTERN_TEMPLATE(struct, PANTRProgressInfo, EigenConfigl);
#ifdef ALPAQA_WITH_QUAD_PRECISION
ALPAQA_EXPORT_EXTERN_TEMPLATE(struct, PANTRProgressInfo, EigenConfigq);
#endif

} // namespace alpaqa
