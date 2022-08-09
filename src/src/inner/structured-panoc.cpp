#include <alpaqa/inner/src/structured-panoc.tpp>

namespace alpaqa {

// clang-format off
ALPAQA_EXPORT_TEMPLATE(struct, StructuredPANOCLBFGSParams, DefaultConfig);
ALPAQA_EXPORT_TEMPLATE(struct, StructuredPANOCLBFGSParams, EigenConfigf);
ALPAQA_EXPORT_TEMPLATE(struct, StructuredPANOCLBFGSParams, EigenConfigd);
ALPAQA_EXPORT_TEMPLATE(struct, StructuredPANOCLBFGSParams, EigenConfigl);
#ifdef ALPAQA_WITH_QUAD_PRECISION
ALPAQA_EXPORT_TEMPLATE(struct, StructuredPANOCLBFGSParams, EigenConfigq);
#endif

ALPAQA_EXPORT_TEMPLATE(struct, StructuredPANOCLBFGSStats, DefaultConfig);
ALPAQA_EXPORT_TEMPLATE(struct, StructuredPANOCLBFGSStats, EigenConfigf);
ALPAQA_EXPORT_TEMPLATE(struct, StructuredPANOCLBFGSStats, EigenConfigd);
ALPAQA_EXPORT_TEMPLATE(struct, StructuredPANOCLBFGSStats, EigenConfigl);
#ifdef ALPAQA_WITH_QUAD_PRECISION
ALPAQA_EXPORT_TEMPLATE(struct, StructuredPANOCLBFGSStats, EigenConfigq);
#endif

ALPAQA_EXPORT_TEMPLATE(struct, StructuredPANOCLBFGSProgressInfo, DefaultConfig);
ALPAQA_EXPORT_TEMPLATE(struct, StructuredPANOCLBFGSProgressInfo, EigenConfigf);
ALPAQA_EXPORT_TEMPLATE(struct, StructuredPANOCLBFGSProgressInfo, EigenConfigd);
ALPAQA_EXPORT_TEMPLATE(struct, StructuredPANOCLBFGSProgressInfo, EigenConfigl);
#ifdef ALPAQA_WITH_QUAD_PRECISION
ALPAQA_EXPORT_TEMPLATE(struct, StructuredPANOCLBFGSProgressInfo, EigenConfigq);
#endif

ALPAQA_EXPORT_TEMPLATE(class, StructuredPANOCLBFGSSolver, DefaultConfig);
ALPAQA_EXPORT_TEMPLATE(class, StructuredPANOCLBFGSSolver, EigenConfigf);
ALPAQA_EXPORT_TEMPLATE(class, StructuredPANOCLBFGSSolver, EigenConfigd);
ALPAQA_EXPORT_TEMPLATE(class, StructuredPANOCLBFGSSolver, EigenConfigl);
#ifdef ALPAQA_WITH_QUAD_PRECISION
ALPAQA_EXPORT_TEMPLATE(class, StructuredPANOCLBFGSSolver, EigenConfigq);
#endif
// clang-format on

} // namespace alpaqa