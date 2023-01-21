#include <alpaqa/implementation/outer/internal/alm-helpers.tpp>

namespace alpaqa::detail {

ALPAQA_EXPORT_TEMPLATE(struct, ALMHelpers, DefaultConfig);
ALPAQA_EXPORT_TEMPLATE(struct, ALMHelpers, EigenConfigf);
ALPAQA_EXPORT_TEMPLATE(struct, ALMHelpers, EigenConfigd);
ALPAQA_EXPORT_TEMPLATE(struct, ALMHelpers, EigenConfigl);
#ifdef ALPAQA_WITH_QUAD_PRECISION
ALPAQA_EXPORT_TEMPLATE(struct, ALMHelpers, EigenConfigq);
#endif

} // namespace alpaqa