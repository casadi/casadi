#pragma once

#include <panoc-alm/inner/decl/lbfgs-stepsize.hpp>
#include <panoc-alm/inner/decl/panoc-fwd.hpp>
#include <panoc-alm/inner/decl/panoc-stop-crit.hpp>
#include <panoc-alm/inner/directions/decl/lbfgs.hpp>
#include <panoc-alm/util/atomic_stop_signal.hpp>
#include <panoc-alm/util/problem.hpp>
#include <panoc-alm/util/solverstatus.hpp>

#include <atomic>
#include <chrono>
#include <limits>
#include <string>

namespace pa {

/// Tuning parameters for the second order PANOC algorithm.
struct SecondOrderPANOCLBFGSParams {
    struct {
        /// Initial estimate of the Lipschitz constant of ∇ψ(x)
        real_t L₀ = 0;
        /// Relative step size for initial finite difference Lipschitz estimate.
        real_t ε = 1e-6;
        /// Minimum step size for initial finite difference Lipschitz estimate.
        real_t δ = 1e-12;
        /// Factor that relates step size γ and Lipschitz constant.
        real_t Lγ_factor = 0.95;
    } Lipschitz; ///< Parameters related to the Lipschitz constant estimate
                 ///  and step size.

    /// Maximum number of inner PANOC iterations.
    unsigned max_iter = 100;
    /// Maximum duration.
    std::chrono::microseconds max_time = std::chrono::minutes(5);
    /// Minimum weight factor between Newton step and projected gradient step.
    real_t τ_min = 1. / 256;
    /// Minimum step size.
    real_t γ_min = 1e-30;
    /// Factor used in update for exponentially weighted nonmonotone line search.
    /// Zero means monotone line search.
    real_t nonmonotone_linesearch = 0;
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

    bool hessian_vec_finited_differences = true;
    bool full_augmented_hessian          = true;

    LBFGSStepSize lbfgs_stepsize = LBFGSStepSize::BasedOnCurvature;
};

struct SecondOrderPANOCLBFGSProgressInfo {
    unsigned k;
    crvec x;
    crvec p;
    real_t norm_sq_p;
    crvec x_hat;
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
    const SecondOrderPANOCLBFGSParams &params;
};

/// Second order PANOC solver for ALM.
/// @ingroup    grp_InnerSolvers
class SecondOrderPANOCLBFGSSolver {
  public:
    using Params = SecondOrderPANOCLBFGSParams;

    struct Stats {
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

    using ProgressInfo = SecondOrderPANOCLBFGSProgressInfo;

    SecondOrderPANOCLBFGSSolver(Params params, LBFGSParams lbfgsparams)
        : params(params), lbfgs(lbfgsparams) {}

    Stats operator()(const Problem &problem,        // in
                     crvec Σ,                       // in
                     real_t ε,                      // in
                     bool always_overwrite_results, // in
                     rvec x,                        // inout
                     rvec y,                        // inout
                     rvec err_z);                   // out

    SecondOrderPANOCLBFGSSolver &
    set_progress_callback(std::function<void(const ProgressInfo &)> cb) {
        this->progress_cb = cb;
        return *this;
    }

    std::string get_name() const { return "SecondOrderPANOCSolverLBFGS"; }

    void stop() { stop_signal.stop(); }

    const Params &get_params() const { return params; }

  private:
    Params params;
    AtomicStopSignal stop_signal;
    std::function<void(const ProgressInfo &)> progress_cb;

  public:
    LBFGS lbfgs;
};

template <class InnerSolver>
struct InnerStatsAccumulator;

template <>
struct InnerStatsAccumulator<SecondOrderPANOCLBFGSSolver> {
    std::chrono::microseconds elapsed_time;
    unsigned iterations          = 0;
    unsigned linesearch_failures = 0;
    unsigned lbfgs_failures      = 0;
    unsigned lbfgs_rejected      = 0;
    unsigned τ_1_accepted        = 0;
    unsigned count_τ             = 0;
    real_t sum_τ                 = 0;
};

inline InnerStatsAccumulator<SecondOrderPANOCLBFGSSolver> &
operator+=(InnerStatsAccumulator<SecondOrderPANOCLBFGSSolver> &acc,
           const SecondOrderPANOCLBFGSSolver::Stats s) {
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

} // namespace pa
