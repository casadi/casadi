#pragma once

#include <alpaqa/inner/panoc.hpp>
#include "alpaqa/problem/dynamics.hpp"

#include <chrono>
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
    /// How often to use a Gauss-Newton step. Zero to disable GN entirely.
    unsigned gn_interval        = 1;
    bool gn_sticky              = true;
    bool reset_lbfgs_on_gn_step = false;
    bool lqr_factor_cholesky    = true;

    /// When to print progress. If set to zero, nothing will be printed.
    /// If set to N != 0, progress is printed every N iterations.
    unsigned print_interval = 0;
    /// The precision of the floating point values printed by the solver.
    int print_precision = std::numeric_limits<real_t>::max_digits10 / 2;

    real_t quadratic_upperbound_tolerance_factor =
        10 * std::numeric_limits<real_t>::epsilon();

    bool update_lipschitz_in_linesearch = true;
};

template <Config Conf = DefaultConfig>
struct PANOCOCPProgressInfo {
    USING_ALPAQA_CONFIG(Conf);

    unsigned k;
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
    const TypeErasedControlProblem<config_t> &problem;
    const PANOCOCPParams<config_t> &params;
};

template <Config Conf = DefaultConfig>
struct PANOCOCPStats {
    USING_ALPAQA_CONFIG(Conf);

    SolverStatus status = SolverStatus::Busy;
    real_t ε            = inf<config_t>;
    std::chrono::nanoseconds elapsed_time{};
    std::chrono::nanoseconds time_forward{};
    std::chrono::nanoseconds time_backward{};
    std::chrono::nanoseconds time_backward_jacobians{};
    std::chrono::nanoseconds time_hessians{};
    std::chrono::nanoseconds time_indices{};
    std::chrono::nanoseconds time_lqr_factor{};
    std::chrono::nanoseconds time_lqr_solve{};
    std::chrono::nanoseconds time_lbfgs_indices{};
    std::chrono::nanoseconds time_lbfgs_apply{};
    std::chrono::nanoseconds time_lbfgs_update{};
    unsigned iterations          = 0;
    unsigned linesearch_failures = 0;
    unsigned lbfgs_failures      = 0;
    unsigned lbfgs_rejected      = 0;
    unsigned τ_1_accepted        = 0;
    unsigned count_τ             = 0;
    real_t sum_τ                 = 0;
    real_t final_γ               = 0;
};

template <Config Conf>
class PANOCOCPSolver {
  public:
    USING_ALPAQA_CONFIG(Conf);

    using Problem      = alpaqa::TypeErasedControlProblem<config_t>;
    using Params       = PANOCOCPParams<config_t>;
    using Stats        = PANOCOCPStats<config_t>;
    using ProgressInfo = PANOCOCPProgressInfo<config_t>;

    PANOCOCPSolver(const Params &params) : params(params) {}

    Stats operator()(const Problem &problem, // in
                     real_t ε,               // in
                     rvec u);                // inout

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
};

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
