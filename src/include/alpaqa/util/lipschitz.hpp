#pragma once

#include <alpaqa/util/vec.hpp>

namespace alpaqa {

struct LipschitzEstimateParams {
    /// Initial estimate of the Lipschitz constant of ∇ψ(x)
    real_t L₀ = 0;
    /// Relative step size for initial finite difference Lipschitz estimate.
    real_t ε = 1e-6;
    /// Minimum step size for initial finite difference Lipschitz estimate.
    real_t δ = 1e-12;
    /// Factor that relates step size γ and Lipschitz constant.
    real_t Lγ_factor = 0.95;
};

} // namespace alpaqa