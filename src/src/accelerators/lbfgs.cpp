#include <alpaqa/config/config.hpp>

#include <alpaqa/accelerators/src/lbfgs.tpp>

namespace alpaqa {

ALPAQA_EXPORT_TEMPLATE(struct, CBFGSParams, DefaultConfig);
ALPAQA_EXPORT_TEMPLATE(struct, CBFGSParams, EigenConfigf);
ALPAQA_EXPORT_TEMPLATE(struct, CBFGSParams, EigenConfigd);
ALPAQA_EXPORT_TEMPLATE(struct, CBFGSParams, EigenConfigl);
#ifdef ALPAQA_WITH_QUAD_PRECISION
ALPAQA_EXPORT_TEMPLATE(struct, CBFGSParams, EigenConfigq);
#endif

ALPAQA_EXPORT_TEMPLATE(struct, LBFGSParams, DefaultConfig);
ALPAQA_EXPORT_TEMPLATE(struct, LBFGSParams, EigenConfigf);
ALPAQA_EXPORT_TEMPLATE(struct, LBFGSParams, EigenConfigd);
ALPAQA_EXPORT_TEMPLATE(struct, LBFGSParams, EigenConfigl);
#ifdef ALPAQA_WITH_QUAD_PRECISION
ALPAQA_EXPORT_TEMPLATE(struct, LBFGSParams, EigenConfigq);
#endif

ALPAQA_EXPORT_TEMPLATE(class, LBFGS, DefaultConfig);
ALPAQA_EXPORT_TEMPLATE(class, LBFGS, EigenConfigf);
ALPAQA_EXPORT_TEMPLATE(class, LBFGS, EigenConfigd);
ALPAQA_EXPORT_TEMPLATE(class, LBFGS, EigenConfigl);
#ifdef ALPAQA_WITH_QUAD_PRECISION
ALPAQA_EXPORT_TEMPLATE(class, LBFGS, EigenConfigq);
#endif

} // namespace alpaqa