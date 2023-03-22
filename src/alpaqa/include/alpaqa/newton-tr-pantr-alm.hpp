#pragma once

#include <alpaqa/inner/directions/pantr/newton-tr.hpp>
#include <alpaqa/inner/pantr.hpp>
#include <alpaqa/outer/alm.hpp>

namespace alpaqa {

// clang-format off
ALPAQA_EXPORT_EXTERN_TEMPLATE(class, PANTRSolver, NewtonTRDirection<DefaultConfig>);
ALPAQA_EXPORT_EXTERN_TEMPLATE(class, PANTRSolver, NewtonTRDirection<EigenConfigf>);
ALPAQA_EXPORT_EXTERN_TEMPLATE(class, PANTRSolver, NewtonTRDirection<EigenConfigd>);
ALPAQA_EXPORT_EXTERN_TEMPLATE(class, PANTRSolver, NewtonTRDirection<EigenConfigl>);
#ifdef ALPAQA_WITH_QUAD_PRECISION
ALPAQA_EXPORT_EXTERN_TEMPLATE(class, PANTRSolver, NewtonTRDirection<EigenConfigq>);
#endif

ALPAQA_EXPORT_EXTERN_TEMPLATE(class, ALMSolver, PANTRSolver<NewtonTRDirection<DefaultConfig>>);
ALPAQA_EXPORT_EXTERN_TEMPLATE(class, ALMSolver, PANTRSolver<NewtonTRDirection<EigenConfigf>>);
ALPAQA_EXPORT_EXTERN_TEMPLATE(class, ALMSolver, PANTRSolver<NewtonTRDirection<EigenConfigd>>);
ALPAQA_EXPORT_EXTERN_TEMPLATE(class, ALMSolver, PANTRSolver<NewtonTRDirection<EigenConfigl>>);
#ifdef ALPAQA_WITH_QUAD_PRECISION
ALPAQA_EXPORT_EXTERN_TEMPLATE(class, ALMSolver, PANTRSolver<NewtonTRDirection<EigenConfigq>>);
#endif
// clang-format on

} // namespace alpaqa
