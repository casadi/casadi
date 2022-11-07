#pragma once

namespace alpaqa {

enum class CostStructure {
    General = 1,
    Separable,
    Quadratic,
    SeparableQuadratic,
};

} // namespace alpaqa