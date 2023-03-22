#pragma once

#include <alpaqa/config/config.hpp>
#include <alpaqa/problem/type-erased-problem.hpp>

namespace alpaqa {

/// This class outlines the interface for direction providers used by PANOC-like
/// algorithms.
///
/// @ingroup grp_DirectionProviders
template <Config Conf>
struct PANOCDirection {
    USING_ALPAQA_CONFIG(Conf);
    using Problem = TypeErasedProblem<config_t>;

    /// Initialize the direction provider.
    ///
    /// @param[in]  problem
    ///             Problem description.
    /// @param[in]  y
    ///             Lagrange multipliers.
    /// @param[in]  Σ
    ///             Penalty factors.
    /// @param[in]  γ_0
    ///             Initial step size.
    /// @param[in]  x_0
    ///             Initial iterate.
    /// @param[in]  x̂_0
    ///             Result of proximal gradient step in @p x_0.
    /// @param[in]  p_0
    ///             Proximal gradient step in @p x_0.
    /// @param[in]  grad_ψx_0
    ///             Gradient of the objective in @p x_0.
    ///
    /// The references @p problem, @p y and @p Σ are guaranteed to remain valid
    /// for subsequent calls to @ref update, @ref apply, @ref changed_γ and
    /// @ref reset.
    void initialize(const Problem &problem, crvec y, crvec Σ, real_t γ_0,
                    crvec x_0, crvec x̂_0, crvec p_0, crvec grad_ψx_0) = delete;

    /// Return whether a direction is available on the very first iteration,
    /// before the first call to @ref update.
    bool has_initial_direction() const = delete;

    /// Update the direction provider when accepting the next iterate.
    ///
    /// @param[in]  γₖ
    ///             Current step size.
    /// @param[in]  γₙₑₓₜ
    ///             Step size for the next iterate.
    /// @param[in]  xₖ
    ///             Current iterate.
    /// @param[in]  xₙₑₓₜ
    ///             Next iterate.
    /// @param[in]  pₖ
    ///             Proximal gradient step in the current iterate.
    /// @param[in]  pₙₑₓₜ
    ///             Proximal gradient step in the next iterate.
    /// @param[in]  grad_ψxₖ
    ///             Gradient of the objective in the current iterate.
    /// @param[in]  grad_ψxₙₑₓₜ
    ///             Gradient of the objective in the next iterate.
    bool update(real_t γₖ, real_t γₙₑₓₜ, crvec xₖ, crvec xₙₑₓₜ, crvec pₖ,
                crvec pₙₑₓₜ, crvec grad_ψxₖ, crvec grad_ψxₙₑₓₜ) = delete;

    /// Apply the direction estimation in the current point.
    ///
    /// @param[in]  γₖ
    ///             Current step size.
    /// @param[in]  xₖ
    ///             Current iterate.
    /// @param[in]  x̂ₖ
    ///             Result of proximal gradient step in @p xₖ.
    /// @param[in]  pₖ
    ///             Proximal gradient step in @p xₖ.
    /// @param[in]  grad_ψxₖ
    ///             Gradient of the objective at @p xₖ.
    /// @param[out] qₖ
    ///             Resulting step.
    bool apply(real_t γₖ, crvec xₖ, crvec x̂ₖ, crvec pₖ, crvec grad_ψxₖ,
               rvec qₖ) const = delete;

    /// Called when the PANOC step size changes.
    void changed_γ(real_t γₖ, real_t old_γₖ) = delete;

    /// Called when using the direction failed. A possible behavior could be to
    /// flush the buffers, hopefully yielding a better direction on the next
    /// iteration.
    void reset() = delete;

    /// Get a human-readable name for this direction provider.
    std::string get_name() const = delete;
};

} // namespace alpaqa
