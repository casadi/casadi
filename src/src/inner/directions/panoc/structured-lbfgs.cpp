#include <alpaqa/inner/directions/panoc/src/structured-lbfgs.tpp>

namespace alpaqa {

ALPAQA_EXPORT_TEMPLATE(struct, StructuredLBFGS, DefaultConfig);
ALPAQA_EXPORT_TEMPLATE(struct, StructuredLBFGS, EigenConfigf);
ALPAQA_EXPORT_TEMPLATE(struct, StructuredLBFGS, EigenConfigd);
ALPAQA_EXPORT_TEMPLATE(struct, StructuredLBFGS, EigenConfigl);
#ifdef ALPAQA_WITH_QUAD_PRECISION
ALPAQA_EXPORT_TEMPLATE(struct, StructuredLBFGS, EigenConfigq);
#endif

} // namespace alpaqa