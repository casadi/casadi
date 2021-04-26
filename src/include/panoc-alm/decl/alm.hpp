#pragma once

#include <panoc-alm/inner/decl/panoc-fwd.hpp>
#include <panoc-alm/util/problem.hpp>
#include <panoc-alm/util/solverstatus.hpp>

#include <chrono>
#include <string>

namespace pa {

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
    real_t Σₘₐₓ = 1e9;
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
    bool preconditioning = true;
};

/// Augmented Lagrangian Method solver
template <class InnerSolverT = PANOCSolver<>>
class ALMSolver {
  public:
    using Params      = ALMParams;
    using InnerSolver = InnerSolverT;

    struct Stats {
        unsigned outer_iterations = 0;
        std::chrono::microseconds elapsed_time;
        unsigned initial_penalty_reduced    = 0;
        unsigned penalty_reduced            = 0;
        unsigned inner_convergence_failures = 0;
        real_t ε                            = inf;
        real_t δ                            = inf;
        real_t norm_penalty                 = 0;

        SolverStatus status = SolverStatus::Unknown;

        InnerStatsAccumulator<InnerSolver> inner;
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

} // namespace pa
