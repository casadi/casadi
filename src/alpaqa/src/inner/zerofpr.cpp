#include <alpaqa/implementation/inner/zerofpr.tpp>

namespace alpaqa {

ALPAQA_EXPORT_TEMPLATE(struct, ZeroFPRParams, DefaultConfig);
ALPAQA_EXPORT_TEMPLATE(struct, ZeroFPRParams, EigenConfigf);
ALPAQA_EXPORT_TEMPLATE(struct, ZeroFPRParams, EigenConfigd);
ALPAQA_EXPORT_TEMPLATE(struct, ZeroFPRParams, EigenConfigl);
#ifdef ALPAQA_WITH_QUAD_PRECISION
ALPAQA_EXPORT_TEMPLATE(struct, ZeroFPRParams, EigenConfigq);
#endif

ALPAQA_EXPORT_TEMPLATE(struct, ZeroFPRStats, DefaultConfig);
ALPAQA_EXPORT_TEMPLATE(struct, ZeroFPRStats, EigenConfigf);
ALPAQA_EXPORT_TEMPLATE(struct, ZeroFPRStats, EigenConfigd);
ALPAQA_EXPORT_TEMPLATE(struct, ZeroFPRStats, EigenConfigl);
#ifdef ALPAQA_WITH_QUAD_PRECISION
ALPAQA_EXPORT_TEMPLATE(struct, ZeroFPRStats, EigenConfigq);
#endif

ALPAQA_EXPORT_TEMPLATE(struct, ZeroFPRProgressInfo, DefaultConfig);
ALPAQA_EXPORT_TEMPLATE(struct, ZeroFPRProgressInfo, EigenConfigf);
ALPAQA_EXPORT_TEMPLATE(struct, ZeroFPRProgressInfo, EigenConfigd);
ALPAQA_EXPORT_TEMPLATE(struct, ZeroFPRProgressInfo, EigenConfigl);
#ifdef ALPAQA_WITH_QUAD_PRECISION
ALPAQA_EXPORT_TEMPLATE(struct, ZeroFPRProgressInfo, EigenConfigq);
#endif

} // namespace alpaqa