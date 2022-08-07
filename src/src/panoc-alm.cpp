#include <alpaqa/inner/src/panoc.tpp>
#include <alpaqa/outer/src/alm.tpp>
#include <alpaqa/panoc-alm.hpp>

namespace alpaqa {

ALPAQA_EXPORT_TEMPLATE(class, PANOCSolver, LBFGS<DefaultConfig>);
ALPAQA_EXPORT_TEMPLATE(class, PANOCSolver, LBFGS<EigenConfigf>);
ALPAQA_EXPORT_TEMPLATE(class, PANOCSolver, LBFGS<EigenConfigd>);
ALPAQA_EXPORT_TEMPLATE(class, PANOCSolver, LBFGS<EigenConfigl>);
#ifdef ALPAQA_WITH_QUAD_PRECISION
ALPAQA_EXPORT_TEMPLATE(class, PANOCSolver, LBFGS<EigenConfigq>);
#endif

ALPAQA_EXPORT_TEMPLATE(class, ALMSolver, PANOCSolver<LBFGS<DefaultConfig>>);
ALPAQA_EXPORT_TEMPLATE(class, ALMSolver, PANOCSolver<LBFGS<EigenConfigf>>);
ALPAQA_EXPORT_TEMPLATE(class, ALMSolver, PANOCSolver<LBFGS<EigenConfigd>>);
ALPAQA_EXPORT_TEMPLATE(class, ALMSolver, PANOCSolver<LBFGS<EigenConfigl>>);
#ifdef ALPAQA_WITH_QUAD_PRECISION
ALPAQA_EXPORT_TEMPLATE(class, ALMSolver, PANOCSolver<LBFGS<EigenConfigq>>);
#endif

} // namespace alpaqa