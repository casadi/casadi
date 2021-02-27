#pragma once

#include "panoc-alm/problem.hpp"
#include "panoc-alm/solverstatus.hpp"
#include "panoc.hpp"
#include <iostream>

namespace pa {

struct ALMParams {
    real_t ε    = 1e-5; ///< Primal tolerance.
    real_t δ    = 1e-5; ///< Dual tolerance.
    real_t Δ    = 10;   ///< Factor used in updating the penalty parameters.
    real_t Σ₀   = 1;    ///< Initial penalty parameter.
    real_t σ₀   = 20;   ///< Initial penalty parameter factor.
    real_t ε₀   = 1;    ///< Initial primal tolerance.
    real_t θ    = 0.25; ///< Error tolerance for penalty increase
    real_t ρ    = 1e-1; ///< Update factor for primal tolerance.
    real_t M    = 1e9;
    real_t σₘₐₓ = 1e15;
    /// Maximum number of outer ALM iterations.
    unsigned int max_iter = 100;
    /// Maximum duration.
    std::chrono::microseconds max_time = std::chrono::minutes(5);

    bool preconditioning = true;

    void verify();
};

/// Augmented Lagrangian Method solver
class ALMSolver {
  public:
    using Params = ALMParams;

    struct Stats {
        unsigned inner_iterations = 0;
        unsigned outer_iterations = 0;
        std::chrono::microseconds elapsed_time;
        unsigned inner_convergence_failures = 0;
        unsigned inner_linesearch_failures  = 0;
        unsigned inner_lbfgs_failures       = 0;
        unsigned inner_lbfgs_rejected       = 0;
        real_t ε = std::numeric_limits<real_t>::infinity();
        real_t δ = std::numeric_limits<real_t>::infinity();

        SolverStatus status = SolverStatus::Unknown;
    };

    ALMSolver(Params params, PANOCSolver::Params panoc_params)
        : params(params), panoc(panoc_params) {}

    Stats operator()(const Problem &problem, vec &y, vec &x);

  private:
    Params params;
    PANOCSolver panoc;
};

} // namespace pa
