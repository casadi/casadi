#include "params.hpp"

PARAMS_TABLE_DEF(alpaqa::LBFGSDirectionParams<Conf>,    //
                 PARAMS_MEMBER(rescale_on_step_size_changes), //
);

PARAMS_TABLE_INST(alpaqa::LBFGSDirectionParams<alpaqa::EigenConfigf>);
PARAMS_TABLE_INST(alpaqa::LBFGSDirectionParams<alpaqa::EigenConfigd>);
PARAMS_TABLE_INST(alpaqa::LBFGSDirectionParams<alpaqa::EigenConfigl>);
#ifdef ALPAQA_WITH_QUAD_PRECISION
PARAMS_TABLE_INST(alpaqa::LBFGSDirectionParams<alpaqa::EigenConfigq>);
#endif