#include <alpaqa/implementation/problem/type-erased-problem.tpp>

namespace alpaqa {

ALPAQA_EXPORT_TEMPLATE(struct, ProblemVTable, DefaultConfig);
ALPAQA_EXPORT_TEMPLATE(struct, ProblemVTable, EigenConfigf);
ALPAQA_EXPORT_TEMPLATE(struct, ProblemVTable, EigenConfigd);
ALPAQA_EXPORT_TEMPLATE(struct, ProblemVTable, EigenConfigl);
#ifdef ALPAQA_WITH_QUAD_PRECISION
ALPAQA_EXPORT_TEMPLATE(struct, ProblemVTable, EigenConfigq);
#endif

} // namespace alpaqa