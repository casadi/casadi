#pragma once

#include "panoc.hpp"

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

class ALMSolver {
  public:
    using Params = ALMParams;

    ALMSolver(Params params, PANOCSolver::Params panoc_params)
        : params(params), panoc(panoc_params) {}

    void operator()(const Problem &problem, vec &y, vec &x) & {
        vec Σ(problem.m);
        vec z(problem.m);
        vec error(problem.m);
        vec error_old(problem.m);

        // Initialize the penalty weights
        Σ.fill(params.Σ₀);
        real_t ε = params.ε₀;

        for (unsigned int i = 0; i < params.max_iter; ++i) {
            detail::project_y(y, problem.D.lowerbound, problem.D.upperbound,
                              params.M);
            panoc(problem, x, z, y, error, Σ, ε);
            real_t norm_e = detail::norm_inf(error);

            if (ε <= params.ε && norm_e <= params.δ) {
                return;
            }
            update_penalty_weights(i == 0, error, error_old, norm_e, Σ);
            ε = params.ρ * ε;
        }
    }

    void update_penalty_weights(bool first_it, vec &e, vec &old_e,
                                real_t norm_e, vec &Σ) {
        for (Eigen::Index i = 0; i < e.rows(); ++i) {
            if (first_it || std::abs(e(i)) > params.θ * std::abs(old_e(i))) {
                Σ(i) = std::fmin(
                    params.σₘₐₓ,
                    std::fmax(params.Δ * std::abs(e(i)) / norm_e, 1) * Σ(i));
            }
        }
    }

  private:
    Params params;
    PANOCSolver panoc;
};

} // namespace pa
