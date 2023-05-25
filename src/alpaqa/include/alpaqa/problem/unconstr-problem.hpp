#pragma once

#include <alpaqa/config/config.hpp>

namespace alpaqa {

/// Implements common problem functions for minimization problems without
/// constraints. Meant to be used as a base class for custom problem
/// implementations.
/// @ingroup grp_Problems
template <Config Conf>
class UnconstrProblem {
  public:
    USING_ALPAQA_CONFIG(Conf);
    /// Number of constraints
    length_t get_m() const { return 0; }

    void eval_g(crvec, rvec) const {}
    void eval_grad_g_prod(crvec, crvec, rvec grad) const { grad.setZero(); }
    void eval_jac_g(crvec, rindexvec, rindexvec, rvec) const {}
    void eval_grad_gi(crvec, index_t, rvec grad_gi) const { grad_gi.setZero(); }

    /// @see @ref TypeErasedProblem::eval_prox_grad_step
    real_t eval_prox_grad_step(real_t γ, crvec x, crvec grad_ψ, rvec x̂, rvec p) const {
        p = -γ * grad_ψ;
        x̂ = x + p;
        return 0;
    }

    /// @see @ref TypeErasedProblem::eval_proj_diff_g
    void eval_proj_diff_g(crvec, rvec) const {}

    /// @see @ref TypeErasedProblem::eval_proj_multipliers
    void eval_proj_multipliers(rvec, real_t) const {}
};

} // namespace alpaqa
