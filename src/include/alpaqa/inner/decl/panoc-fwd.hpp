#pragma once

#include <alpaqa/inner/directions/decl/lbfgs-fwd.hpp>
namespace pa {
template <class DirectionProviderT = LBFGS>
class PANOCSolver;
template <class InnerSolverStats>
struct InnerStatsAccumulator;
}