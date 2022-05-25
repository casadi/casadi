#pragma once

#include <alpaqa/config/config.hpp>
#include <alpaqa/problem/box.hpp>

namespace alpaqa {

template <class DirectionProviderT>
struct PANOCDirection {
    USING_ALPAQA_CONFIG_TEMPLATE(DirectionProviderT::config_t);

    static void initialize(crvec x₀, crvec x̂₀, crvec p₀, crvec grad₀) = delete;

    static bool update(crvec xₖ, crvec xₖ₊₁, crvec pₖ, crvec pₖ₊₁,
                       crvec grad_new, const Box<config_t> &C,
                       real_t γ_new) = delete;

    /// Apply the direction estimation in the current point.
    /// @param[in]  xₖ
    ///             Current iterate.
    /// @param[in]  x̂ₖ
    ///             Result of proximal gradient step in current iterate.
    /// @param[in]  pₖ
    ///             Proximal gradient step between x̂ₖ and xₖ.
    /// @param[in]  γ
    ///             H₀ = γI for L-BFGS
    /// @param[out] qₖ
    ///             Resulting step.
    static bool apply(crvec xₖ, crvec x̂ₖ, crvec pₖ, real_t γ, rvec qₖ) = delete;

    static void changed_γ(real_t γₖ, real_t old_γₖ) = delete;

    static void reset() = delete;
};

} // namespace alpaqa
