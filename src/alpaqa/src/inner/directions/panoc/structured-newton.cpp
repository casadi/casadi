#include <alpaqa/inner/directions/panoc/structured-newton.hpp>

namespace alpaqa {

ALPAQA_EXPORT_TEMPLATE(struct, StructuredNewtonDirection, DefaultConfig);
ALPAQA_EXPORT_TEMPLATE(struct, StructuredNewtonDirection, EigenConfigf);
ALPAQA_EXPORT_TEMPLATE(struct, StructuredNewtonDirection, EigenConfigd);
ALPAQA_EXPORT_TEMPLATE(struct, StructuredNewtonDirection, EigenConfigl);
#ifdef ALPAQA_WITH_QUAD_PRECISION
ALPAQA_EXPORT_TEMPLATE(struct, StructuredNewtonDirection, EigenConfigq);
#endif

} // namespace alpaqa