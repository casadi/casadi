#include <alpaqa/implementation/lbfgsb-adapter.tpp>

namespace alpaqa::lbfgspp {

ALPAQA_LBFGSPP_EXPORT_TEMPLATE(struct, LBFGSBStats, DefaultConfig);
ALPAQA_LBFGSPP_EXPORT_TEMPLATE(struct, LBFGSBStats, EigenConfigf);
ALPAQA_LBFGSPP_EXPORT_TEMPLATE(struct, LBFGSBStats, EigenConfigd);
ALPAQA_LBFGSPP_EXPORT_TEMPLATE(struct, LBFGSBStats, EigenConfigl);
#ifdef ALPAQA_WITH_QUAD_PRECISION
ALPAQA_LBFGSPP_EXPORT_TEMPLATE(struct, LBFGSBStats, EigenConfigq);
#endif

ALPAQA_LBFGSPP_EXPORT_TEMPLATE(class, LBFGSBSolver, DefaultConfig);
ALPAQA_LBFGSPP_EXPORT_TEMPLATE(class, LBFGSBSolver, EigenConfigf);
ALPAQA_LBFGSPP_EXPORT_TEMPLATE(class, LBFGSBSolver, EigenConfigd);
ALPAQA_LBFGSPP_EXPORT_TEMPLATE(class, LBFGSBSolver, EigenConfigl);
#ifdef ALPAQA_WITH_QUAD_PRECISION
ALPAQA_LBFGSPP_EXPORT_TEMPLATE(class, LBFGSBSolver, EigenConfigq);
#endif

} // namespace alpaqa::lbfgspp