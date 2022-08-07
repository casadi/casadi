#include <alpaqa/inner/src/panoc.tpp>

namespace alpaqa {

ALPAQA_EXPORT_TEMPLATE(struct, PANOCParams, DefaultConfig);
ALPAQA_EXPORT_TEMPLATE(struct, PANOCParams, EigenConfigf);
ALPAQA_EXPORT_TEMPLATE(struct, PANOCParams, EigenConfigd);
ALPAQA_EXPORT_TEMPLATE(struct, PANOCParams, EigenConfigl);
#ifdef ALPAQA_WITH_QUAD_PRECISION
ALPAQA_EXPORT_TEMPLATE(struct, PANOCParams, EigenConfigq);
#endif

ALPAQA_EXPORT_TEMPLATE(struct, PANOCStats, DefaultConfig);
ALPAQA_EXPORT_TEMPLATE(struct, PANOCStats, EigenConfigf);
ALPAQA_EXPORT_TEMPLATE(struct, PANOCStats, EigenConfigd);
ALPAQA_EXPORT_TEMPLATE(struct, PANOCStats, EigenConfigl);
#ifdef ALPAQA_WITH_QUAD_PRECISION
ALPAQA_EXPORT_TEMPLATE(struct, PANOCStats, EigenConfigq);
#endif

ALPAQA_EXPORT_TEMPLATE(struct, PANOCProgressInfo, DefaultConfig);
ALPAQA_EXPORT_TEMPLATE(struct, PANOCProgressInfo, EigenConfigf);
ALPAQA_EXPORT_TEMPLATE(struct, PANOCProgressInfo, EigenConfigd);
ALPAQA_EXPORT_TEMPLATE(struct, PANOCProgressInfo, EigenConfigl);
#ifdef ALPAQA_WITH_QUAD_PRECISION
ALPAQA_EXPORT_TEMPLATE(struct, PANOCProgressInfo, EigenConfigq);
#endif

} // namespace alpaqa