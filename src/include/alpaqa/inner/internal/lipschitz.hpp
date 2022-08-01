#pragma once

#include <alpaqa/config/config.hpp>

namespace alpaqa {

template <Config Conf = DefaultConfig>
struct LipschitzEstimateParams {
    USING_ALPAQA_CONFIG(Conf);

    /// Initial estimate of the Lipschitz constant of ∇ψ(x)
    real_t L_0 = 0;
    /// Relative step size for initial finite difference Lipschitz estimate.
    real_t ε = 1e-6;
    /// Minimum step size for initial finite difference Lipschitz estimate.
    real_t δ = 1e-12;
    /// Factor that relates step size γ and Lipschitz constant.
    real_t Lγ_factor = 0.95;

    void verify() const {
        // TODO
    }
};

} // namespace alpaqa