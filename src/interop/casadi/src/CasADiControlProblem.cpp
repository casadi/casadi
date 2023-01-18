#include <alpaqa/config/config.hpp>
#include <alpaqa/implementation/casadi/CasADiControlProblem.tpp>

namespace alpaqa {
CASADI_OCP_LOADER_EXPORT_TEMPLATE(class, CasADiControlProblem, EigenConfigd);
CASADI_OCP_LOADER_EXPORT_TEMPLATE(class, CasADiControlProblem, DefaultConfig);
} // namespace alpaqa