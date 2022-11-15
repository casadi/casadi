#include <alpaqa/config/config.hpp>
#include <alpaqa/interop/casadi/src/CasADiQuadraticControlProblem.tpp>

namespace alpaqa {
CASADI_LOADER_EXPORT_TEMPLATE(class, CasADiQuadraticControlProblem,
                              EigenConfigd);
CASADI_LOADER_EXPORT_TEMPLATE(class, CasADiQuadraticControlProblem,
                              DefaultConfig);
} // namespace alpaqa