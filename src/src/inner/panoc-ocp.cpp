#include <alpaqa/inner/src/panoc-ocp.tpp>

namespace alpaqa {

ALPAQA_EXPORT_TEMPLATE(class, PANOCOCPSolver, DefaultConfig);
ALPAQA_EXPORT_TEMPLATE(class, PANOCOCPSolver, EigenConfigf);
ALPAQA_EXPORT_TEMPLATE(class, PANOCOCPSolver, EigenConfigd);
ALPAQA_EXPORT_TEMPLATE(class, PANOCOCPSolver, EigenConfigl);
#ifdef ALPAQA_WITH_QUAD_PRECISION
ALPAQA_EXPORT_TEMPLATE(class, PANOCOCPSolver, EigenConfigq);
#endif

} // namespace alpaqa