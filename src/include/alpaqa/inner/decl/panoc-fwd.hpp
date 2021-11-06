#pragma once

#include <alpaqa/inner/directions/decl/lbfgs-fwd.hpp>
namespace alpaqa {
template <class DirectionProviderT = LBFGS>
class PANOCSolver;
template <class InnerSolverStats>
struct InnerStatsAccumulator;
}