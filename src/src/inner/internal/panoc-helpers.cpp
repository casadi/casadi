#include <alpaqa/inner/src/panoc-helpers.tpp>

namespace alpaqa::detail {

ALPAQA_EXPORT_TEMPLATE(struct, PANOCHelpers, DefaultConfig);
ALPAQA_EXPORT_TEMPLATE(struct, PANOCHelpers, EigenConfigf);
ALPAQA_EXPORT_TEMPLATE(struct, PANOCHelpers, EigenConfigd);
ALPAQA_EXPORT_TEMPLATE(struct, PANOCHelpers, EigenConfigl);
#ifdef ALPAQA_WITH_QUAD_PRECISION
ALPAQA_EXPORT_TEMPLATE(struct, PANOCHelpers, EigenConfigq);
#endif

} // namespace alpaqa::detail