#include <alpaqa/implementation/inner/zerofpr.tpp>
#include <alpaqa/implementation/outer/alm.tpp>
#include <alpaqa/zerofpr-alm.hpp>

namespace alpaqa {

// clang-format off
ALPAQA_EXPORT_TEMPLATE(class, ZeroFPRSolver, LBFGSDirection<DefaultConfig>);
ALPAQA_EXPORT_TEMPLATE(class, ZeroFPRSolver, LBFGSDirection<EigenConfigf>);
ALPAQA_EXPORT_TEMPLATE(class, ZeroFPRSolver, LBFGSDirection<EigenConfigd>);
ALPAQA_EXPORT_TEMPLATE(class, ZeroFPRSolver, LBFGSDirection<EigenConfigl>);
#ifdef ALPAQA_WITH_QUAD_PRECISION
ALPAQA_EXPORT_TEMPLATE(class, ZeroFPRSolver, LBFGSDirection<EigenConfigq>);
#endif

ALPAQA_EXPORT_TEMPLATE(class, ALMSolver, ZeroFPRSolver<LBFGSDirection<DefaultConfig>>);
ALPAQA_EXPORT_TEMPLATE(class, ALMSolver, ZeroFPRSolver<LBFGSDirection<EigenConfigf>>);
ALPAQA_EXPORT_TEMPLATE(class, ALMSolver, ZeroFPRSolver<LBFGSDirection<EigenConfigd>>);
ALPAQA_EXPORT_TEMPLATE(class, ALMSolver, ZeroFPRSolver<LBFGSDirection<EigenConfigl>>);
#ifdef ALPAQA_WITH_QUAD_PRECISION
ALPAQA_EXPORT_TEMPLATE(class, ALMSolver, ZeroFPRSolver<LBFGSDirection<EigenConfigq>>);
#endif
// clang-format on

} // namespace alpaqa