#include <alpaqa/inner/src/structured-panoc.tpp>
#include <alpaqa/outer/src/alm.tpp>
#include <alpaqa/structured-panoc-alm.hpp>

namespace alpaqa {

// clang-format off
ALPAQA_EXPORT_TEMPLATE(class, ALMSolver, StructuredPANOCLBFGSSolver<DefaultConfig>);
ALPAQA_EXPORT_TEMPLATE(class, ALMSolver, StructuredPANOCLBFGSSolver<EigenConfigf>);
ALPAQA_EXPORT_TEMPLATE(class, ALMSolver, StructuredPANOCLBFGSSolver<EigenConfigd>);
ALPAQA_EXPORT_TEMPLATE(class, ALMSolver, StructuredPANOCLBFGSSolver<EigenConfigl>);
#ifdef ALPAQA_WITH_QUAD_PRECISION
ALPAQA_EXPORT_TEMPLATE(class, ALMSolver, StructuredPANOCLBFGSSolver<EigenConfigq>);
#endif
// clang-format on

} // namespace alpaqa