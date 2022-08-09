#pragma once

#include <alpaqa/inner/directions/panoc/lbfgs.hpp>
#include <alpaqa/inner/panoc.hpp>
#include <alpaqa/outer/alm.hpp>

namespace alpaqa {

// clang-format off
ALPAQA_EXPORT_EXTERN_TEMPLATE(class, PANOCSolver, LBFGS<DefaultConfig>);
ALPAQA_EXPORT_EXTERN_TEMPLATE(class, PANOCSolver, LBFGS<EigenConfigf>);
ALPAQA_EXPORT_EXTERN_TEMPLATE(class, PANOCSolver, LBFGS<EigenConfigd>);
ALPAQA_EXPORT_EXTERN_TEMPLATE(class, PANOCSolver, LBFGS<EigenConfigl>);
#ifdef ALPAQA_WITH_QUAD_PRECISION
ALPAQA_EXPORT_EXTERN_TEMPLATE(class, PANOCSolver, LBFGS<EigenConfigq>);
#endif

ALPAQA_EXPORT_EXTERN_TEMPLATE(class, ALMSolver, PANOCSolver<LBFGS<DefaultConfig>>);
ALPAQA_EXPORT_EXTERN_TEMPLATE(class, ALMSolver, PANOCSolver<LBFGS<EigenConfigf>>);
ALPAQA_EXPORT_EXTERN_TEMPLATE(class, ALMSolver, PANOCSolver<LBFGS<EigenConfigd>>);
ALPAQA_EXPORT_EXTERN_TEMPLATE(class, ALMSolver, PANOCSolver<LBFGS<EigenConfigl>>);
#ifdef ALPAQA_WITH_QUAD_PRECISION
ALPAQA_EXPORT_EXTERN_TEMPLATE(class, ALMSolver, PANOCSolver<LBFGS<EigenConfigq>>);
#endif

} // namespace alpaqa
