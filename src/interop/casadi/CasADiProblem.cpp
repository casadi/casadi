#include <alpaqa/config/config.hpp>
#include <alpaqa/interop/casadi/src/CasADiProblem.tpp>

namespace alpaqa {
CASADI_LOADER_EXPORT_TEMPLATE(class, CasADiProblem, EigenConfigd);
CASADI_LOADER_EXPORT_TEMPLATE(class, CasADiProblem, DefaultConfig);
} // namespace alpaqa