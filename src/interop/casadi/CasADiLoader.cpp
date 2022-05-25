#include <alpaqa/config/config.hpp>
#include <alpaqa/interop/casadi/src/CasADiLoader.tpp>

namespace alpaqa {
CASADI_LOADER_EXPORT_TEMPLATE(class, CasADiProblem, EigenConfigd);
CASADI_LOADER_EXPORT_TEMPLATE(class, CasADiProblem, DefaultConfig);
} // namespace alpaqa