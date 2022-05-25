#include <alpaqa/accelerators/lbfgs.hpp>
#include <alpaqa/inner/directions/panoc/lbfgs.hpp>
#include <alpaqa/inner/src/panoc.tpp>
#include <alpaqa/outer/src/alm.tpp>

#include <alpaqa/config/config.hpp>

namespace alpaqa {

ALPAQA_EXPORT_TEMPLATE(struct, ALMParams, DefaultConfig);
ALPAQA_EXPORT_TEMPLATE(struct, ALMParams, EigenConfigf);
ALPAQA_EXPORT_TEMPLATE(struct, ALMParams, EigenConfigd);
ALPAQA_EXPORT_TEMPLATE(struct, ALMParams, EigenConfigl);
#ifdef ALPAQA_WITH_QUAD_PRECISION
ALPAQA_EXPORT_TEMPLATE(struct, ALMParams, EigenConfigq);
#endif

ALPAQA_EXPORT_TEMPLATE(class, ALMSolver, PANOCSolver<LBFGS<DefaultConfig>>);
ALPAQA_EXPORT_TEMPLATE(class, ALMSolver, PANOCSolver<LBFGS<EigenConfigf>>);
ALPAQA_EXPORT_TEMPLATE(class, ALMSolver, PANOCSolver<LBFGS<EigenConfigd>>);
ALPAQA_EXPORT_TEMPLATE(class, ALMSolver, PANOCSolver<LBFGS<EigenConfigl>>);
#ifdef ALPAQA_WITH_QUAD_PRECISION
ALPAQA_EXPORT_TEMPLATE(class, ALMSolver, PANOCSolver<LBFGS<EigenConfigq>>);
#endif

} // namespace alpaqa