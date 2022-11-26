#include <alpaqa/inner/src/panoc.tpp>
#include <alpaqa/outer/src/alm.tpp>
#include <alpaqa/panoc-alm.hpp>

namespace alpaqa {

// clang-format off
ALPAQA_EXPORT_TEMPLATE(class, PANOCSolver, LBFGSDirection<DefaultConfig>);
ALPAQA_EXPORT_TEMPLATE(class, PANOCSolver, LBFGSDirection<EigenConfigf>);
ALPAQA_EXPORT_TEMPLATE(class, PANOCSolver, LBFGSDirection<EigenConfigd>);
ALPAQA_EXPORT_TEMPLATE(class, PANOCSolver, LBFGSDirection<EigenConfigl>);
#ifdef ALPAQA_WITH_QUAD_PRECISION
ALPAQA_EXPORT_TEMPLATE(class, PANOCSolver, LBFGSDirection<EigenConfigq>);
#endif

ALPAQA_EXPORT_TEMPLATE(class, ALMSolver, PANOCSolver<LBFGSDirection<DefaultConfig>>);
ALPAQA_EXPORT_TEMPLATE(class, ALMSolver, PANOCSolver<LBFGSDirection<EigenConfigf>>);
ALPAQA_EXPORT_TEMPLATE(class, ALMSolver, PANOCSolver<LBFGSDirection<EigenConfigd>>);
ALPAQA_EXPORT_TEMPLATE(class, ALMSolver, PANOCSolver<LBFGSDirection<EigenConfigl>>);
#ifdef ALPAQA_WITH_QUAD_PRECISION
ALPAQA_EXPORT_TEMPLATE(class, ALMSolver, PANOCSolver<LBFGSDirection<EigenConfigq>>);
#endif
// clang-format on

} // namespace alpaqa