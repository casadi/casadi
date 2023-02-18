#pragma once

#include <alpaqa/inner/directions/panoc/structured-lbfgs.hpp>
#include <alpaqa/inner/panoc.hpp>
#include <alpaqa/outer/alm.hpp>

namespace alpaqa {

// clang-format off
ALPAQA_EXPORT_EXTERN_TEMPLATE(class, PANOCSolver, StructuredLBFGSDirection<DefaultConfig>);
ALPAQA_EXPORT_EXTERN_TEMPLATE(class, PANOCSolver, StructuredLBFGSDirection<EigenConfigf>);
ALPAQA_EXPORT_EXTERN_TEMPLATE(class, PANOCSolver, StructuredLBFGSDirection<EigenConfigd>);
ALPAQA_EXPORT_EXTERN_TEMPLATE(class, PANOCSolver, StructuredLBFGSDirection<EigenConfigl>);
#ifdef ALPAQA_WITH_QUAD_PRECISION
ALPAQA_EXPORT_EXTERN_TEMPLATE(class, PANOCSolver, StructuredLBFGSDirection<EigenConfigq>);
#endif

ALPAQA_EXPORT_EXTERN_TEMPLATE(class, ALMSolver, PANOCSolver<StructuredLBFGSDirection<DefaultConfig>>);
ALPAQA_EXPORT_EXTERN_TEMPLATE(class, ALMSolver, PANOCSolver<StructuredLBFGSDirection<EigenConfigf>>);
ALPAQA_EXPORT_EXTERN_TEMPLATE(class, ALMSolver, PANOCSolver<StructuredLBFGSDirection<EigenConfigd>>);
ALPAQA_EXPORT_EXTERN_TEMPLATE(class, ALMSolver, PANOCSolver<StructuredLBFGSDirection<EigenConfigl>>);
#ifdef ALPAQA_WITH_QUAD_PRECISION
ALPAQA_EXPORT_EXTERN_TEMPLATE(class, ALMSolver, PANOCSolver<StructuredLBFGSDirection<EigenConfigq>>);
#endif
// clang-format on

} // namespace alpaqa
