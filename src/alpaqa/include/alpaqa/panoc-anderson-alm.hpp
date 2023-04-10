#pragma once

#include <alpaqa/inner/directions/panoc/anderson.hpp>
#include <alpaqa/inner/panoc.hpp>
#include <alpaqa/outer/alm.hpp>

namespace alpaqa {

// clang-format off
ALPAQA_EXPORT_EXTERN_TEMPLATE(class, PANOCSolver, AndersonDirection<DefaultConfig>);
ALPAQA_EXPORT_EXTERN_TEMPLATE(class, PANOCSolver, AndersonDirection<EigenConfigf>);
ALPAQA_EXPORT_EXTERN_TEMPLATE(class, PANOCSolver, AndersonDirection<EigenConfigd>);
ALPAQA_EXPORT_EXTERN_TEMPLATE(class, PANOCSolver, AndersonDirection<EigenConfigl>);
#ifdef ALPAQA_WITH_QUAD_PRECISION
ALPAQA_EXPORT_EXTERN_TEMPLATE(class, PANOCSolver, AndersonDirection<EigenConfigq>);
#endif

ALPAQA_EXPORT_EXTERN_TEMPLATE(class, ALMSolver, PANOCSolver<AndersonDirection<DefaultConfig>>);
ALPAQA_EXPORT_EXTERN_TEMPLATE(class, ALMSolver, PANOCSolver<AndersonDirection<EigenConfigf>>);
ALPAQA_EXPORT_EXTERN_TEMPLATE(class, ALMSolver, PANOCSolver<AndersonDirection<EigenConfigd>>);
ALPAQA_EXPORT_EXTERN_TEMPLATE(class, ALMSolver, PANOCSolver<AndersonDirection<EigenConfigl>>);
#ifdef ALPAQA_WITH_QUAD_PRECISION
ALPAQA_EXPORT_EXTERN_TEMPLATE(class, ALMSolver, PANOCSolver<AndersonDirection<EigenConfigq>>);
#endif
// clang-format on

} // namespace alpaqa
