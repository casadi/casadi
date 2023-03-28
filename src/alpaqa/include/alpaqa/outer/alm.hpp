#pragma once

#include <alpaqa/config/config.hpp>
#include <alpaqa/export.hpp>
#include <alpaqa/inner/internal/solverstatus.hpp>
#include <alpaqa/outer/internal/alm-helpers.hpp>
#include <alpaqa/problem/type-erased-problem.hpp>

#include <chrono>
#include <iostream>
#include <string>

namespace alpaqa {

template <class InnerSolverStats>
struct InnerStatsAccumulator;

/// Parameters for the Augmented Lagrangian solver.
template <Config Conf = DefaultConfig>
struct ALMParams {
    USING_ALPAQA_CONFIG(Conf);

    /// Primal tolerance.
    real_t ε = 1e-5;
    /// Dual tolerance.
    real_t δ = 1e-5;
    /// Factor used in updating the penalty parameters.
    real_t Δ = 10;
    /// Factor to reduce @ref Δ when inner convergence fails.
    real_t Δ_lower = 0.8;
    /// Minimum value for @ref Δ after reduction because of inner convergence
    /// failure.
    real_t Δ_min = 1.1;
    /// Initial penalty parameter. When set to zero (which is the default),
    /// it is computed automatically, based on the constraint violation in the
    /// starting point and the parameter @ref σ_0.
    /// @todo Change default to 0.
    real_t Σ_0 = 1;
    /// Initial penalty parameter factor. Active if @ref Σ_0 is set
    /// to zero.
    real_t σ_0 = 20;
    /// Factor to reduce the initial penalty factor by if convergence fails in
    /// in the first iteration.
    real_t Σ_0_lower = 0.6;
    /// Initial primal tolerance.
    real_t ε_0 = 1;
    /// Factor to increase the initial primal tolerance if convergence fails in
    /// the first iteration.
    real_t ε_0_increase = 1.1;
    /// Update factor for primal tolerance.
    real_t ρ = 1e-1;
    /// Factor to increase the primal tolerance update factor by if convergence
    /// fails.
    real_t ρ_increase = 2;
    /// Maximum value of @ref ρ after increase because of inner convergence
    /// failure.
    real_t ρ_max = 0.5;
    /// Error tolerance for penalty increase
    real_t θ = 0.1;
    /// Lagrange multiplier bound.
    real_t M = 1e9;
    /// Maximum penalty factor.
    real_t Σ_max = 1e9;
    /// Minimum penalty factor (used during initialization).
    real_t Σ_min = 1e-9;
    /// Maximum number of outer ALM iterations.
    unsigned int max_iter = 100;
    /// Maximum duration.
    std::chrono::nanoseconds max_time = std::chrono::minutes(5);

    /// How many times can the initial penalty @ref Σ_0 or
    /// @ref σ_0 and the initial primal tolerance @ref ε_0
    /// be reduced.
    unsigned max_num_initial_retries = 0;
    /// How many times can the penalty update factor @ref Δ and the
    /// primal tolerance factor @ref ρ be reduced.
    unsigned max_num_retries = 0;
    /// Combined limit for @ref max_num_initial_retries and
    /// @ref max_num_retries.
    unsigned max_total_num_retries = 0;

    /// When to print progress. If set to zero, nothing will be printed.
    /// If set to N != 0, progress is printed every N iterations.
    unsigned print_interval = 0;
    /// The precision of the floating point values printed by the solver.
    int print_precision = std::numeric_limits<real_t>::max_digits10 / 2;

    /// Use one penalty factor for all m constraints.
    bool single_penalty_factor = false;
};

/// Augmented Lagrangian Method solver
///
/// @ingroup    grp_ALMSolver
template <class InnerSolverT>
class ALMSolver {
  public:
    USING_ALPAQA_CONFIG_TEMPLATE(InnerSolverT::config_t);

    using Params      = ALMParams<config_t>;
    using InnerSolver = InnerSolverT;
    using Problem     = typename InnerSolver::Problem;

    struct Stats {
        /// Total number of outer ALM iterations (i.e. the number of times that
        /// the inner solver was invoked).
        unsigned outer_iterations = 0;
        /// Total elapsed time.
        std::chrono::nanoseconds elapsed_time{};
        /// The number of times that the initial penalty factor was reduced by
        /// @ref ALMParams::Σ_0_lower and that the initial tolerance was
        /// increased by @ref ALMParams::ε_0_increase because the inner solver
        /// failed to converge in the first ALM iteration. If this number is
        /// greater than zero, consider lowering the initial penalty factor
        /// @ref ALMParams::Σ_0 or @ref ALMParams::σ_0 or increasing the initial
        /// tolerance @ref ALMParams::ε_0 (or both).
        unsigned initial_penalty_reduced = 0;
        /// The number of times that the penalty update factor @ref ALMParams::Δ
        /// was reduced, that the tolerance update factor @ref ALMParams::ρ was
        /// increased, and that an ALM iteration had to be restarted with a
        /// lower penalty factor and a higher tolerance because the inner solver
        /// failed to converge (not applicable in the first ALM iteration).
        /// If this number is greater than zero, consider lowerering the
        /// penalty update factor @ref ALMParams::Δ or increasing the tolerance
        /// update factor (or both). Lowering the initial penalty factor could
        /// help as well.
        unsigned penalty_reduced = 0;
        /// The total number of times that the inner solver failed to converge.
        unsigned inner_convergence_failures = 0;
        /// Final primal tolerance that was reached, depends on the stopping
        /// criterion used by the inner solver, see for example
        /// @ref PANOCStopCrit.
        real_t ε = inf<config_t>;
        /// Final dual tolerance or constraint violation that was reached:
        /// @f[
        /// \delta = \| \Pi_D\left(g(x^k) + \Sigma^{-1}y\right) \|_\infty
        /// @f]
        real_t δ = inf<config_t>;
        /// 2-norm of the final penalty factors @f$ \| \Sigma \|_2 @f$.
        real_t norm_penalty = 0;

        /// Whether the solver converged or not.
        /// @see @ref SolverStatus
        SolverStatus status = SolverStatus::Busy;

        /// The statistics of the inner solver invocations, accumulated over all
        /// ALM iterations.
        InnerStatsAccumulator<typename InnerSolver::Stats> inner;
    };

    ALMSolver(Params params, InnerSolver &&inner_solver)
        : params(params),
          inner_solver(std::forward<InnerSolver>(inner_solver)) {}
    ALMSolver(Params params, const InnerSolver &inner_solver)
        : params(params), inner_solver(inner_solver) {}

    Stats operator()(const Problem &problem, rvec x, rvec y);
    template <class P>
    Stats operator()(const P &problem, rvec x, rvec y) {
        return operator()(Problem::template make<P>(problem), x, y);
    }

    std::string get_name() const {
        return "ALMSolver<" + inner_solver.get_name() + ">";
    }

    /// Abort the computation and return the result so far.
    /// Can be called from other threads or signal handlers.
    void stop() { inner_solver.stop(); }

    const Params &get_params() const { return params; }

  private:
    Params params;
    using Helpers = detail::ALMHelpers<config_t>;

  public:
    InnerSolver inner_solver;
    std::ostream *os = &std::cout;
};

ALPAQA_EXPORT_EXTERN_TEMPLATE(struct, ALMParams, DefaultConfig);
ALPAQA_EXPORT_EXTERN_TEMPLATE(struct, ALMParams, EigenConfigf);
ALPAQA_EXPORT_EXTERN_TEMPLATE(struct, ALMParams, EigenConfigd);
ALPAQA_EXPORT_EXTERN_TEMPLATE(struct, ALMParams, EigenConfigl);
#ifdef ALPAQA_WITH_QUAD_PRECISION
ALPAQA_EXPORT_EXTERN_TEMPLATE(struct, ALMParams, EigenConfigq);
#endif

} // namespace alpaqa
