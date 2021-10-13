#pragma once

#include <panoc-alm/inner/directions/decl/lbfgs-fwd.hpp>
namespace pa {
template <class DirectionProviderT = LBFGS>
class PANOCSolver;
template <class InnerSolverStats>
struct InnerStatsAccumulator;
}