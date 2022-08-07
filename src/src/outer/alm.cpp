#include <alpaqa/outer/alm.hpp>

namespace alpaqa {

ALPAQA_EXPORT_TEMPLATE(struct, ALMParams, DefaultConfig);
ALPAQA_EXPORT_TEMPLATE(struct, ALMParams, EigenConfigf);
ALPAQA_EXPORT_TEMPLATE(struct, ALMParams, EigenConfigd);
ALPAQA_EXPORT_TEMPLATE(struct, ALMParams, EigenConfigl);
#ifdef ALPAQA_WITH_QUAD_PRECISION
ALPAQA_EXPORT_TEMPLATE(struct, ALMParams, EigenConfigq);
#endif

} // namespace alpaqa