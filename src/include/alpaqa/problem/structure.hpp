#pragma once

namespace alpaqa {

enum class CostStructure {
    General = 1,
    DiagonalHessian,
    Quadratic,
    DiagonalQuadratic,
};

} // namespace alpaqa