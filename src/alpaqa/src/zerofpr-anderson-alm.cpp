#include <alpaqa/implementation/inner/zerofpr.tpp>
#include <alpaqa/implementation/outer/alm.tpp>
#include <alpaqa/zerofpr-anderson-alm.hpp>

namespace alpaqa {

// clang-format off
ALPAQA_EXPORT_TEMPLATE(class, ZeroFPRSolver, AndersonDirection<DefaultConfig>);
ALPAQA_EXPORT_TEMPLATE(class, ZeroFPRSolver, AndersonDirection<EigenConfigf>);
ALPAQA_EXPORT_TEMPLATE(class, ZeroFPRSolver, AndersonDirection<EigenConfigd>);
ALPAQA_EXPORT_TEMPLATE(class, ZeroFPRSolver, AndersonDirection<EigenConfigl>);
#ifdef ALPAQA_WITH_QUAD_PRECISION
ALPAQA_EXPORT_TEMPLATE(class, ZeroFPRSolver, AndersonDirection<EigenConfigq>);
#endif

ALPAQA_EXPORT_TEMPLATE(class, ALMSolver, ZeroFPRSolver<AndersonDirection<DefaultConfig>>);
ALPAQA_EXPORT_TEMPLATE(class, ALMSolver, ZeroFPRSolver<AndersonDirection<EigenConfigf>>);
ALPAQA_EXPORT_TEMPLATE(class, ALMSolver, ZeroFPRSolver<AndersonDirection<EigenConfigd>>);
ALPAQA_EXPORT_TEMPLATE(class, ALMSolver, ZeroFPRSolver<AndersonDirection<EigenConfigl>>);
#ifdef ALPAQA_WITH_QUAD_PRECISION
ALPAQA_EXPORT_TEMPLATE(class, ALMSolver, ZeroFPRSolver<AndersonDirection<EigenConfigq>>);
#endif
// clang-format on

} // namespace alpaqa