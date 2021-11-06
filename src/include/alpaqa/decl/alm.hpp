#pragma once

#include <alpaqa/inner/decl/panoc-fwd.hpp>
#include <alpaqa/util/problem.hpp>
#include <alpaqa/util/solverstatus.hpp>

#include <chrono>
#include <string>

namespace alpaqa {

/// Parameters for the Augmented Lagrangian solver.
struct ALMParams {
    /// Primal tolerance.
    real_t ε = 1e-5;
    /// Dual tolerance.
    real_t δ = 1e-5;
    /// Factor used in updating the penalty parameters.
    real_t Δ = 10;
    /// Factor to reduce @ref ALMParams::Δ when inner convergence fails.
    real_t Δ_lower = 0.8;
    /// Initial penalty parameter. When set to zero (which is the default),
    /// it is computed automatically, based on the constraint violation in the
    /// starting point and the parameter @ref ALMParams::σ₀.
    /// @todo Change default to 0.
    real_t Σ₀ = 1;
    /// Initial penalty parameter factor. Active if @ref ALMParams::Σ₀ is set
    /// to zero.
    real_t σ₀ = 20;
    /// Factor to reduce the initial penalty factor by if convergence fails in
    /// in the first iteration.
    real_t Σ₀_lower = 0.6;
    /// Initial primal tolerance.
    real_t ε₀ = 1;
    /// Factor to increase the initial primal tolerance if convergence fails in
    /// the first iteration.
    real_t ε₀_increase = 1.1;
    /// Update factor for primal tolerance.
    real_t ρ = 1e-1;
    /// Factor to increase the primal tolerance update factor by if convergence
    /// fails.
    real_t ρ_increase = 2;
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
    std::chrono::microseconds max_time = std::chrono::minutes(5);

    /// How many times can the initial penalty @ref ALMParams::Σ₀ or
    /// @ref ALMParams::σ₀ and the initial primal tolerance @ref ALMParams::ε₀
    /// be reduced.
    unsigned max_num_initial_retries = 20;
    /// How many times can the penalty update factor @ref ALMParams::Δ and the
    /// primal tolerance factor @ref ALMParams::ρ be reduced.
    unsigned max_num_retries = 20;
    /// Combined limit for @ref ALMParams::max_num_initial_retries and
    /// @ref ALMParams::max_num_retries.
    unsigned max_total_num_retries = 40;

    /// When to print progress. If set to zero, nothing will be printed.
    /// If set to N != 0, progress is printed every N iterations.
    unsigned print_interval = 0;

    /// Apply preconditioning to the problem, based on the gradients in the
    /// starting point.
    bool preconditioning = false;
    /// Use one penalty factor for all m constraints.
    bool single_penalty_factor = false;
};

/// Augmented Lagrangian Method solver
///
/// @ingroup    grp_ALMSolver
template <class InnerSolverT = PANOCSolver<>>
class ALMSolver {
  public:
    using Params      = ALMParams;
    using InnerSolver = InnerSolverT;

    struct Stats {
        /// Total number of outer ALM iterations (i.e. the number of times that
        /// the inner solver was invoked).
        unsigned outer_iterations = 0;
        /// Total elapsed time.
        std::chrono::microseconds elapsed_time;
        /// The number of times that the initial penalty factor was reduced by
        /// @ref ALMParams::Σ₀_lower and that the initial tolerance was 
        /// increased by @ref ALMParams::ε₀_increase because the inner solver 
        /// failed to converge in the first ALM iteration. If this number is 
        /// greater than zero, consider lowering the initial penalty factor 
        /// @ref ALMParams::Σ₀ or @ref ALMParams::σ₀ or increasing the initial 
        /// tolerance @ref ALMParams::ε₀ (or both).
        unsigned initial_penalty_reduced    = 0;
        /// The number of times that the penalty update factor @ref ALMParams::Δ
        /// was reduced, that the tolerance update factor @ref ALMParams::ρ was 
        /// increased, and that an ALM iteration had to be restarted with a 
        /// lower penalty factor and a higher tolerance because the inner solver
        /// failed to converge (not applicable in the first ALM iteration).
        /// If this number is greater than zero, consider lowerering the 
        /// penalty update factor @ref ALMParams::Δ or increasing the tolerance
        /// update factor (or both). Lowering the initial penalty factor could
        /// help as well.
        unsigned penalty_reduced            = 0;
        /// The total number of times that the inner solver failed to converge.
        unsigned inner_convergence_failures = 0;
        /// Final primal tolerance that was reached, depends on the stopping 
        /// criterion used by the inner solver, see for example 
        /// @ref PANOCStopCrit.
        real_t ε                            = inf;
        /// Final dual tolerance or constraint violation that was reached:
        /// @f[
        /// \delta = \| \Pi_D\left(g(x^k) + \Sigma^{-1}y\right) \|_\infty
        /// @f]
        real_t δ                            = inf;
        /// 2-norm of the final penalty factors @f$ \| \Sigma \|_2 @f$.
        real_t norm_penalty                 = 0;

        /// Whether the solver converged or not.
        /// @see @ref SolverStatus
        SolverStatus status = SolverStatus::Unknown;

        /// The statistics of the inner solver invocations, accumulated over all
        /// ALM iterations.
        InnerStatsAccumulator<typename InnerSolver::Stats> inner;
    };

    ALMSolver(Params params, InnerSolver &&inner_solver)
        : params(params),
          inner_solver(std::forward<InnerSolver>(inner_solver)) {}
    ALMSolver(Params params, const InnerSolver &inner_solver)
        : params(params), inner_solver(inner_solver) {}

    Stats operator()(const Problem &problem, rvec y, rvec x);

    std::string get_name() const {
        return "ALMSolver<" + inner_solver.get_name() + ">";
    }

    /// Abort the computation and return the result so far.
    /// Can be called from other threads or signal handlers.
    void stop() { inner_solver.stop(); }

    const Params &get_params() const { return params; }

  private:
    Params params;

  public:
    InnerSolver inner_solver;
};

} // namespace alpaqa
