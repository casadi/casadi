#include <alpaqa/implementation/inner/pantr.tpp>

namespace alpaqa {

ALPAQA_EXPORT_TEMPLATE(struct, PANTRParams, DefaultConfig);
ALPAQA_EXPORT_TEMPLATE(struct, PANTRParams, EigenConfigf);
ALPAQA_EXPORT_TEMPLATE(struct, PANTRParams, EigenConfigd);
ALPAQA_EXPORT_TEMPLATE(struct, PANTRParams, EigenConfigl);
#ifdef ALPAQA_WITH_QUAD_PRECISION
ALPAQA_EXPORT_TEMPLATE(struct, PANTRParams, EigenConfigq);
#endif

ALPAQA_EXPORT_TEMPLATE(struct, PANTRStats, DefaultConfig);
ALPAQA_EXPORT_TEMPLATE(struct, PANTRStats, EigenConfigf);
ALPAQA_EXPORT_TEMPLATE(struct, PANTRStats, EigenConfigd);
ALPAQA_EXPORT_TEMPLATE(struct, PANTRStats, EigenConfigl);
#ifdef ALPAQA_WITH_QUAD_PRECISION
ALPAQA_EXPORT_TEMPLATE(struct, PANTRStats, EigenConfigq);
#endif

ALPAQA_EXPORT_TEMPLATE(struct, PANTRProgressInfo, DefaultConfig);
ALPAQA_EXPORT_TEMPLATE(struct, PANTRProgressInfo, EigenConfigf);
ALPAQA_EXPORT_TEMPLATE(struct, PANTRProgressInfo, EigenConfigd);
ALPAQA_EXPORT_TEMPLATE(struct, PANTRProgressInfo, EigenConfigl);
#ifdef ALPAQA_WITH_QUAD_PRECISION
ALPAQA_EXPORT_TEMPLATE(struct, PANTRProgressInfo, EigenConfigq);
#endif

} // namespace alpaqa