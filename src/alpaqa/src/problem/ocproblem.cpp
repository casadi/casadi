#include <alpaqa/implementation/problem/ocproblem.tpp>

namespace alpaqa {

ALPAQA_EXPORT_TEMPLATE(struct, ControlProblemVTable, DefaultConfig);
ALPAQA_EXPORT_TEMPLATE(struct, ControlProblemVTable, EigenConfigf);
ALPAQA_EXPORT_TEMPLATE(struct, ControlProblemVTable, EigenConfigd);
ALPAQA_EXPORT_TEMPLATE(struct, ControlProblemVTable, EigenConfigl);
#ifdef ALPAQA_WITH_QUAD_PRECISION
ALPAQA_EXPORT_TEMPLATE(struct, ControlProblemVTable, EigenConfigq);
#endif

} // namespace alpaqa