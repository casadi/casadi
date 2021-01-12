#pragma once

#include "panoc.hpp"
#include <iostream>

namespace pa {

struct ALMParams {
    real_t ε;
    real_t δ;
    real_t Δ;  ///< Factor used in updating the penalty parameters
    real_t Σ₀; ///< Initial penalty parameter
    real_t ε₀; ///< Initial tolerance on x
    real_t ρ;  ///< Update factor for tolerance on x
    real_t θ;
    real_t M;
    real_t σₘₐₓ;
    unsigned int max_iter;

    void verify();
};

/// Augmented Lagrangian Method solver
class ALMSolver {
  public:
    using Params = ALMParams;

    struct Stats {
        unsigned inner_iterations = 0;
        unsigned outer_iterations = 0;
    };

    ALMSolver(Params params, PANOCSolver::Params panoc_params)
        : params(params), panoc(panoc_params) {}

    Stats operator()(const Problem &problem, vec &y, vec &x);

  private:
    Params params;
    PANOCSolver panoc;
};

} // namespace pa
