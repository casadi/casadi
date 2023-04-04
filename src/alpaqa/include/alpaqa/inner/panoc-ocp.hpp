#pragma once

#include <alpaqa/accelerators/lbfgs.hpp>
#include <alpaqa/inner/panoc.hpp>
#include <alpaqa/problem/ocproblem.hpp>

#include <chrono>
#include <iostream>
#include <limits>
#include <string>

namespace alpaqa {

/// Tuning parameters for the PANOC algorithm.
template <Config Conf = DefaultConfig>
struct PANOCOCPParams {
    USING_ALPAQA_CONFIG(Conf);

    /// Parameters related to the Lipschitz constant estimate and step size.
    LipschitzEstimateParams<config_t> Lipschitz;
    /// Maximum number of inner PANOC iterations.
    unsigned max_iter = 100;
    /// Maximum duration.
    std::chrono::nanoseconds max_time = std::chrono::minutes(5);
    /// Minimum weight factor between Newton step and projected gradient step,
    /// line search parameter.
    real_t min_linesearch_coefficient = real_t(1. / 256);
    /// Parameter β used in the line search (see Algorithm 2 in
    /// @cite de_marchi_proximal_2022). @f$ 0 < \beta < 1 @f$
    real_t linesearch_strictness_factor = real_t(0.95);
    /// Minimum Lipschitz constant estimate.
    real_t L_min = real_t(1e-5);
    /// Maximum Lipschitz constant estimate.
    real_t L_max = real_t(1e20);
    /// Maximum number of times to double the Lipschitz constant estimate per
    /// iteration.
    unsigned L_max_inc = 16;
    /// What stopping criterion to use.
    PANOCStopCrit stop_crit = PANOCStopCrit::ApproxKKT;
    /// Maximum number of iterations without any progress before giving up.
    unsigned max_no_progress = 10;
    /// How often to use a Gauss-Newton step. Zero to disable GN entirely.
    unsigned gn_interval        = 1;
    bool gn_sticky              = true;
    bool reset_lbfgs_on_gn_step = false;
    bool lqr_factor_cholesky    = true;

    /// L-BFGS parameters (e.g. memory).
    LBFGSParams<config_t> lbfgs_params;

    /// When to print progress. If set to zero, nothing will be printed.
    /// If set to N != 0, progress is printed every N iterations.
    unsigned print_interval = 0;
    /// The precision of the floating point values printed by the solver.
    int print_precision = std::numeric_limits<real_t>::max_digits10 / 2;

    real_t quadratic_upperbound_tolerance_factor =
        real_t(1e2) * std::numeric_limits<real_t>::epsilon();
    real_t linesearch_tolerance_factor =
        real_t(1e2) * std::numeric_limits<real_t>::epsilon();

    bool disable_acceleration = false;
};

template <Config Conf = DefaultConfig>
struct PANOCOCPProgressInfo {
    USING_ALPAQA_CONFIG(Conf);

    unsigned k;
    SolverStatus status;
    crvec xu;
    crvec p;
    real_t norm_sq_p;
    crvec x̂u;
    real_t φγ;
    real_t ψ;
    crvec grad_ψ;
    real_t ψ_hat;
    crvec q;
    bool gn;
    length_t nJ;
    real_t lqr_min_rcond;
    real_t L;
    real_t γ;
    real_t τ;
    real_t ε;
    unsigned outer_iter;
    const TypeErasedControlProblem<config_t> *problem;
    const PANOCOCPParams<config_t> *params;

    [[nodiscard]] vec u() const;
    [[nodiscard]] vec û() const;
    [[nodiscard]] vec x() const;
    [[nodiscard]] vec x̂() const;
};

template <Config Conf = DefaultConfig>
struct PANOCOCPStats {
    USING_ALPAQA_CONFIG(Conf);

    SolverStatus status = SolverStatus::Busy;
    real_t ε            = inf<config_t>;
    std::chrono::nanoseconds elapsed_time{};
    std::chrono::nanoseconds time_prox{};
    std::chrono::nanoseconds time_forward{};
    std::chrono::nanoseconds time_backward{};
    std::chrono::nanoseconds time_jacobians{};
    std::chrono::nanoseconds time_hessians{};
    std::chrono::nanoseconds time_indices{};
    std::chrono::nanoseconds time_lqr_factor{};
    std::chrono::nanoseconds time_lqr_solve{};
    std::chrono::nanoseconds time_lbfgs_indices{};
    std::chrono::nanoseconds time_lbfgs_apply{};
    std::chrono::nanoseconds time_lbfgs_update{};
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

template <Config Conf>
class PANOCOCPSolver {
  public:
    USING_ALPAQA_CONFIG(Conf);

    using Problem      = alpaqa::TypeErasedControlProblem<config_t>;
    using Params       = PANOCOCPParams<config_t>;
    using Stats        = PANOCOCPStats<config_t>;
    using ProgressInfo = PANOCOCPProgressInfo<config_t>;
    using SolveOptions = InnerSolveOptions<config_t>;

    PANOCOCPSolver(const Params &params) : params(params) {}

    Stats operator()(const Problem &problem,   // in
                     const SolveOptions &opts, // in
                     rvec u,                   // inout
                     rvec y,                   // in
                     crvec μ,                  // in
                     rvec err_z);              // out

    template <class P>
    Stats operator()(const P &problem, const SolveOptions &opts, rvec u, rvec y,
                     crvec μ, rvec e) {
        return operator()(Problem::template make<P>(problem), opts, u, y, μ, e);
    }

    template <class P>
    Stats operator()(const P &problem, const SolveOptions &opts, rvec u) {
        if (problem.get_m() != 0)
            throw std::invalid_argument("Missing arguments y, Σ, e");
        mvec y{nullptr, 0}, Σ{nullptr, 0}, e{nullptr, 0};
        return operator()(problem, opts, u, y, Σ, e);
    }

    /// Specify a callable that is invoked with some intermediate results on
    /// each iteration of the algorithm.
    /// @see @ref ProgressInfo
    PANOCOCPSolver &
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
    std::ostream *os = &std::cout;
};

template <Config Conf>
struct InnerStatsAccumulator<PANOCOCPStats<Conf>> {
    USING_ALPAQA_CONFIG(Conf);

    /// Total elapsed time in the inner solver.
    std::chrono::nanoseconds elapsed_time{};
    /// Total time spent computing proximal mappings.
    std::chrono::nanoseconds time_prox{};
    /// Total time spent doing forward simulations.
    std::chrono::nanoseconds time_forward{};
    /// Total time spent doing backward gradient evaluations.
    std::chrono::nanoseconds time_backward{};
    /// Total time spent computing dynamics Jacobians.
    std::chrono::nanoseconds time_jacobians{};
    /// Total time spent computing cost Hessians and Hessian-vector products.
    std::chrono::nanoseconds time_hessians{};
    /// Total time spent determining active indices.
    std::chrono::nanoseconds time_indices{};
    /// Total time spent performing LQR factorizations.
    std::chrono::nanoseconds time_lqr_factor{};
    /// Total time spent solving the (factorized) LQR problem.
    std::chrono::nanoseconds time_lqr_solve{};
    /// Total time spent determining active indices for L-BFGS applications.
    std::chrono::nanoseconds time_lbfgs_indices{};
    /// Total time spent applying L-BFGS estimates.
    std::chrono::nanoseconds time_lbfgs_apply{};
    /// Total time spent updating the L-BFGS estimate.
    std::chrono::nanoseconds time_lbfgs_update{};
    /// Total time spent in the user-provided progress callback.
    std::chrono::nanoseconds time_progress_callback{};
    /// Total number of inner PANOC iterations.
    unsigned iterations = 0;
    /// Total number of PANOC line search failures.
    unsigned linesearch_failures = 0;
    /// Total number of PANOC line search backtracking steps.
    unsigned linesearch_backtracks = 0;
    /// Total number of PANOC step size reductions.
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
    /// The final PANOC step size γ.
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
InnerStatsAccumulator<PANOCOCPStats<Conf>> &
operator+=(InnerStatsAccumulator<PANOCOCPStats<Conf>> &acc,
           const PANOCOCPStats<Conf> &s) {
    acc.iterations += s.iterations;
    acc.elapsed_time += s.elapsed_time;
    acc.time_prox += s.time_prox;
    acc.time_forward += s.time_forward;
    acc.time_backward += s.time_backward;
    acc.time_jacobians += s.time_jacobians;
    acc.time_hessians += s.time_hessians;
    acc.time_indices += s.time_indices;
    acc.time_lqr_factor += s.time_lqr_factor;
    acc.time_lqr_solve += s.time_lqr_solve;
    acc.time_lbfgs_indices += s.time_lbfgs_indices;
    acc.time_lbfgs_apply += s.time_lbfgs_apply;
    acc.time_lbfgs_update += s.time_lbfgs_update;
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

ALPAQA_EXPORT_EXTERN_TEMPLATE(struct, PANOCOCPProgressInfo, DefaultConfig);
ALPAQA_EXPORT_EXTERN_TEMPLATE(struct, PANOCOCPProgressInfo, EigenConfigf);
ALPAQA_EXPORT_EXTERN_TEMPLATE(struct, PANOCOCPProgressInfo, EigenConfigd);
ALPAQA_EXPORT_EXTERN_TEMPLATE(struct, PANOCOCPProgressInfo, EigenConfigl);
#ifdef ALPAQA_WITH_QUAD_PRECISION
ALPAQA_EXPORT_EXTERN_TEMPLATE(struct, PANOCOCPProgressInfo, EigenConfigq);
#endif

ALPAQA_EXPORT_EXTERN_TEMPLATE(class, PANOCOCPSolver, DefaultConfig);
ALPAQA_EXPORT_EXTERN_TEMPLATE(class, PANOCOCPSolver, EigenConfigf);
ALPAQA_EXPORT_EXTERN_TEMPLATE(class, PANOCOCPSolver, EigenConfigd);
ALPAQA_EXPORT_EXTERN_TEMPLATE(class, PANOCOCPSolver, EigenConfigl);
#ifdef ALPAQA_WITH_QUAD_PRECISION
ALPAQA_EXPORT_EXTERN_TEMPLATE(class, PANOCOCPSolver, EigenConfigq);
#endif

} // namespace alpaqa
