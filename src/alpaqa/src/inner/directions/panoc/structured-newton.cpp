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

#include <alpaqa/implementation/inner/panoc.tpp>

namespace alpaqa {

// clang-format off
ALPAQA_EXPORT_TEMPLATE(class, PANOCSolver, StructuredNewtonDirection<DefaultConfig>);
ALPAQA_EXPORT_TEMPLATE(class, PANOCSolver, StructuredNewtonDirection<EigenConfigf>);
ALPAQA_EXPORT_TEMPLATE(class, PANOCSolver, StructuredNewtonDirection<EigenConfigd>);
ALPAQA_EXPORT_TEMPLATE(class, PANOCSolver, StructuredNewtonDirection<EigenConfigl>);
#ifdef ALPAQA_WITH_QUAD_PRECISION
ALPAQA_EXPORT_TEMPLATE(class, PANOCSolver, StructuredNewtonDirection<EigenConfigq>);
#endif
// clang-format on

} // namespace alpaqa