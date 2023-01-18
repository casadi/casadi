#include <alpaqa/implementation/inner/directions/panoc/structured-lbfgs.tpp>
#include <alpaqa/implementation/inner/panoc.tpp>
#include <alpaqa/implementation/outer/alm.tpp>
#include <alpaqa/structured-panoc-alm.hpp>

namespace alpaqa {

// clang-format off
ALPAQA_EXPORT_TEMPLATE(class, ALMSolver, PANOCSolver<StructuredLBFGSDirection<DefaultConfig>>);
ALPAQA_EXPORT_TEMPLATE(class, ALMSolver, PANOCSolver<StructuredLBFGSDirection<EigenConfigf>>);
ALPAQA_EXPORT_TEMPLATE(class, ALMSolver, PANOCSolver<StructuredLBFGSDirection<EigenConfigd>>);
ALPAQA_EXPORT_TEMPLATE(class, ALMSolver, PANOCSolver<StructuredLBFGSDirection<EigenConfigl>>);
#ifdef ALPAQA_WITH_QUAD_PRECISION
ALPAQA_EXPORT_TEMPLATE(class, ALMSolver, PANOCSolver<StructuredLBFGSDirection<EigenConfigq>>);
#endif
// clang-format on

} // namespace alpaqa