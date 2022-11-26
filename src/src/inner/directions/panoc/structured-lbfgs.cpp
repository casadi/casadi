#include <alpaqa/inner/directions/panoc/src/structured-lbfgs.tpp>

namespace alpaqa {

ALPAQA_EXPORT_TEMPLATE(struct, StructuredLBFGSDirection, DefaultConfig);
ALPAQA_EXPORT_TEMPLATE(struct, StructuredLBFGSDirection, EigenConfigf);
ALPAQA_EXPORT_TEMPLATE(struct, StructuredLBFGSDirection, EigenConfigd);
ALPAQA_EXPORT_TEMPLATE(struct, StructuredLBFGSDirection, EigenConfigl);
#ifdef ALPAQA_WITH_QUAD_PRECISION
ALPAQA_EXPORT_TEMPLATE(struct, StructuredLBFGSDirection, EigenConfigq);
#endif

} // namespace alpaqa