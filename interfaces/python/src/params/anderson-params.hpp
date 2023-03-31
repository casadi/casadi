#pragma once

#include <alpaqa/accelerators/anderson.hpp>
#include <kwargs-to-struct.hpp>

PARAMS_TABLE_DECL(alpaqa::AndersonAccelParams<Conf>);

extern PARAMS_TABLE_INST(alpaqa::AndersonAccelParams<alpaqa::EigenConfigf>);
extern PARAMS_TABLE_INST(alpaqa::AndersonAccelParams<alpaqa::EigenConfigd>);
extern PARAMS_TABLE_INST(alpaqa::AndersonAccelParams<alpaqa::EigenConfigl>);
#ifdef ALPAQA_WITH_QUAD_PRECISION
extern PARAMS_TABLE_INST(alpaqa::AndersonAccelParams<alpaqa::EigenConfigq>);
#endif
