#include <alpaqa/config/config.hpp>
#include <alpaqa/interop/casadi/src/experimental-CasADiControlProblem.tpp>

namespace alpaqa::experimental {
CASADI_LOADER_EXPORT_TEMPLATE(class, CasADiControlProblem, EigenConfigd);
CASADI_LOADER_EXPORT_TEMPLATE(class, CasADiControlProblem, DefaultConfig);
} // namespace alpaqa::experimental