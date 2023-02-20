#include <alpaqa/implementation/inner/directions/panoc/structured-lbfgs.tpp>
#include <alpaqa/implementation/inner/zerofpr.tpp>
#include <alpaqa/implementation/outer/alm.tpp>
#include <alpaqa/structured-zerofpr-alm.hpp>

namespace alpaqa {

// clang-format off
ALPAQA_EXPORT_TEMPLATE(class, ZeroFPRSolver, StructuredLBFGSDirection<DefaultConfig>);
ALPAQA_EXPORT_TEMPLATE(class, ZeroFPRSolver, StructuredLBFGSDirection<EigenConfigf>);
ALPAQA_EXPORT_TEMPLATE(class, ZeroFPRSolver, StructuredLBFGSDirection<EigenConfigd>);
ALPAQA_EXPORT_TEMPLATE(class, ZeroFPRSolver, StructuredLBFGSDirection<EigenConfigl>);
#ifdef ALPAQA_WITH_QUAD_PRECISION
ALPAQA_EXPORT_TEMPLATE(class, ZeroFPRSolver, StructuredLBFGSDirection<EigenConfigq>);
#endif

ALPAQA_EXPORT_TEMPLATE(class, ALMSolver, ZeroFPRSolver<StructuredLBFGSDirection<DefaultConfig>>);
ALPAQA_EXPORT_TEMPLATE(class, ALMSolver, ZeroFPRSolver<StructuredLBFGSDirection<EigenConfigf>>);
ALPAQA_EXPORT_TEMPLATE(class, ALMSolver, ZeroFPRSolver<StructuredLBFGSDirection<EigenConfigd>>);
ALPAQA_EXPORT_TEMPLATE(class, ALMSolver, ZeroFPRSolver<StructuredLBFGSDirection<EigenConfigl>>);
#ifdef ALPAQA_WITH_QUAD_PRECISION
ALPAQA_EXPORT_TEMPLATE(class, ALMSolver, ZeroFPRSolver<StructuredLBFGSDirection<EigenConfigq>>);
#endif
// clang-format on

} // namespace alpaqa