#pragma once

#include <alpaqa/inner/decl/lbfgs-stepsize.hpp>
#include <alpaqa/inner/decl/panoc-fwd.hpp>
#include <alpaqa/inner/decl/panoc-stop-crit.hpp>
#include <alpaqa/inner/directions/decl/panoc-direction-update.hpp>
#include <alpaqa/util/atomic_stop_signal.hpp>
#include <alpaqa/util/lipschitz.hpp>
#include <alpaqa/util/problem.hpp>
#include <alpaqa/util/solverstatus.hpp>

#include <atomic>
#include <chrono>
#include <limits>
#include <string>

namespace alpaqa {

/// Tuning parameters for the PANOC algorithm.
struct PANOCParams {
    /// Parameters related to the Lipschitz constant estimate and step size.
    LipschitzEstimateParams Lipschitz;
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
    /// What stopping criterion to use.
    PANOCStopCrit stop_crit = PANOCStopCrit::ApproxKKT;
    /// Maximum number of iterations without any progress before giving up.
    unsigned max_no_progress = 10;

    /// When to print progress. If set to zero, nothing will be printed.
    /// If set to N != 0, progress is printed every N iterations.
    unsigned print_interval = 0;

    real_t quadratic_upperbound_tolerance_factor =
        10 * std::numeric_limits<real_t>::epsilon();

    bool update_lipschitz_in_linesearch = true;
    bool alternative_linesearch_cond    = false;

    LBFGSStepSize lbfgs_stepsize = LBFGSStepSize::BasedOnCurvature;
};

struct PANOCStats {
    SolverStatus status = SolverStatus::Unknown;
    real_t ε            = inf;
    std::chrono::microseconds elapsed_time;
    unsigned iterations          = 0;
    unsigned linesearch_failures = 0;
    unsigned lbfgs_failures      = 0;
    unsigned lbfgs_rejected      = 0;
    unsigned τ_1_accepted        = 0;
    unsigned count_τ             = 0;
    real_t sum_τ                 = 0;
};

struct PANOCProgressInfo {
    unsigned k;
    crvec x;
    crvec p;
    real_t norm_sq_p;
    crvec x_hat;
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
    const Problem &problem;
    const PANOCParams &params;
};

/// PANOC solver for ALM.
/// @ingroup    grp_InnerSolvers
template <class DirectionProviderT>
class PANOCSolver {
  public:
    using Params            = PANOCParams;
    using DirectionProvider = DirectionProviderT;
    using Stats             = PANOCStats;
    using ProgressInfo      = PANOCProgressInfo;

    PANOCSolver(Params params,
                PANOCDirection<DirectionProvider> &&direction_provider)
        : params(params), direction_provider(std::move(direction_provider)) {}
    PANOCSolver(Params params,
                const PANOCDirection<DirectionProvider> &direction_provider)
        : params(params), direction_provider(direction_provider) {}

    Stats operator()(const Problem &problem,        // in
                     crvec Σ,                       // in
                     real_t ε,                      // in
                     bool always_overwrite_results, // in
                     rvec x,                        // inout
                     rvec y,                        // inout
                     rvec err_z);                   // out

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

  public:
    PANOCDirection<DirectionProvider> direction_provider;
};

template <class InnerSolverStats>
struct InnerStatsAccumulator;

template <>
struct InnerStatsAccumulator<PANOCStats> {
    std::chrono::microseconds elapsed_time;
    unsigned iterations          = 0;
    unsigned linesearch_failures = 0;
    unsigned lbfgs_failures      = 0;
    unsigned lbfgs_rejected      = 0;
    unsigned τ_1_accepted        = 0;
    unsigned count_τ             = 0;
    real_t sum_τ                 = 0;
};

inline InnerStatsAccumulator<PANOCStats> &
operator+=(InnerStatsAccumulator<PANOCStats> &acc, const PANOCStats &s) {
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

} // namespace alpaqa
