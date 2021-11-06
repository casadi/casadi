#pragma once

#include <alpaqa/util/box.hpp>

namespace alpaqa {

template <class DirectionProviderT>
struct PANOCDirection {

    static void initialize(DirectionProviderT &dp, crvec x₀, crvec x̂₀,
                           crvec p₀, crvec grad₀) = delete;

    static bool update(DirectionProviderT &dp, crvec xₖ, crvec xₖ₊₁,
                       crvec pₖ, crvec pₖ₊₁, crvec gradₖ₊₁,
                       const Box &C, real_t γₖ₊₁) = delete;

    /// Apply the direction estimation in the current point.
    /// @param[in]  dp
    ///             Direction provider (e.g. LBFGS, Anderson Acceleration).
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
    static bool apply(DirectionProviderT &dp, crvec xₖ, crvec x̂ₖ,
                      crvec pₖ, real_t γ, rvec qₖ) = delete;

    static void changed_γ(DirectionProviderT &dp, real_t γₖ,
                          real_t old_γₖ) = delete;
};

} // namespace alpaqa
