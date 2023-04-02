#pragma once

#include <alpaqa/config/config.hpp>

namespace alpaqa {

template <Config Conf = DefaultConfig>
struct LipschitzEstimateParams {
    USING_ALPAQA_CONFIG(Conf);

    /// Initial estimate of the Lipschitz constant of ∇ψ(x)
    real_t L_0 = 0;
    /// Relative step size for initial finite difference Lipschitz estimate.
    real_t ε = real_t(1e-6);
    /// Minimum step size for initial finite difference Lipschitz estimate.
    real_t δ = real_t(1e-12);
    /// Factor that relates step size γ and Lipschitz constant.
    /// Parameter α in Algorithm 2 of @cite de_marchi_proximal_2022.
    /// @f$ 0 < \alpha < 1 @f$
    real_t Lγ_factor = real_t(0.95);

    void verify() const {
        // TODO
    }
};

} // namespace alpaqa