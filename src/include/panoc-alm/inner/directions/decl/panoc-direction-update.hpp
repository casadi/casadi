#pragma once

#include <panoc-alm/util/box.hpp>

namespace pa {

template <class DirectionProviderT>
struct PANOCDirection {

    static void initialize(DirectionProviderT &dp, const vec &x₀, const vec &x̂₀,
                           const vec &p₀, const vec &grad₀) = delete;

    static bool update(DirectionProviderT &dp, const vec &xₖ, const vec &xₖ₊₁,
                       const vec &pₖ, const vec &pₖ₊₁, const vec &gradₖ₊₁,
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
    static bool apply(DirectionProviderT &dp, const vec &xₖ, const vec &x̂ₖ,
                      const vec &pₖ, real_t γ, vec &qₖ) = delete;

    static void changed_γ(DirectionProviderT &dp, real_t γₖ,
                          real_t old_γₖ) = delete;
};

} // namespace pa
