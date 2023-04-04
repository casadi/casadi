#include "params.hpp"

PARAMS_TABLE_DEF(alpaqa::AndersonDirectionParams<Conf>, //
                 PARAMS_MEMBER(rescale_on_step_size_changes), //
);

PARAMS_TABLE_INST(alpaqa::AndersonDirectionParams<alpaqa::EigenConfigf>);
PARAMS_TABLE_INST(alpaqa::AndersonDirectionParams<alpaqa::EigenConfigd>);
PARAMS_TABLE_INST(alpaqa::AndersonDirectionParams<alpaqa::EigenConfigl>);
#ifdef ALPAQA_WITH_QUAD_PRECISION
PARAMS_TABLE_INST(alpaqa::AndersonDirectionParams<alpaqa::EigenConfigq>);
#endif
