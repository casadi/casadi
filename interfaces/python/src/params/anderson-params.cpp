#include <params/anderson-params.hpp>

PARAMS_TABLE_DEF(alpaqa::AndersonAccelParams<Conf>, //
                 PARAMS_MEMBER(memory),             //
                 PARAMS_MEMBER(min_div_fac),        //
);

PARAMS_TABLE_INST(alpaqa::AndersonAccelParams<alpaqa::EigenConfigf>);
PARAMS_TABLE_INST(alpaqa::AndersonAccelParams<alpaqa::EigenConfigd>);
PARAMS_TABLE_INST(alpaqa::AndersonAccelParams<alpaqa::EigenConfigl>);
#ifdef ALPAQA_WITH_QUAD_PRECISION
PARAMS_TABLE_INST(alpaqa::AndersonAccelParams<alpaqa::EigenConfigq>);
#endif
