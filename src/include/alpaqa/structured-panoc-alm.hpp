#pragma once

#include <alpaqa/inner/panoc.hpp>
#include <alpaqa/inner/directions/panoc/structured-lbfgs.hpp>
#include <alpaqa/outer/alm.hpp>

namespace alpaqa {

// clang-format off
ALPAQA_EXPORT_EXTERN_TEMPLATE(class, ALMSolver, PANOCSolver<StructuredLBFGS<DefaultConfig>>);
ALPAQA_EXPORT_EXTERN_TEMPLATE(class, ALMSolver, PANOCSolver<StructuredLBFGS<EigenConfigf>>);
ALPAQA_EXPORT_EXTERN_TEMPLATE(class, ALMSolver, PANOCSolver<StructuredLBFGS<EigenConfigd>>);
ALPAQA_EXPORT_EXTERN_TEMPLATE(class, ALMSolver, PANOCSolver<StructuredLBFGS<EigenConfigl>>);
#ifdef ALPAQA_WITH_QUAD_PRECISION
ALPAQA_EXPORT_EXTERN_TEMPLATE(class, ALMSolver, PANOCSolver<StructuredLBFGS<EigenConfigq>>);
#endif
// clang-format on

} // namespace alpaqa
