#pragma once

#include <alpaqa/inner/pantr.hpp>
#include <params/params.hpp>

PARAMS_TABLE_DECL(alpaqa::PANTRParams<Conf>);

extern PARAMS_TABLE_INST(alpaqa::PANTRParams<alpaqa::EigenConfigf>);
extern PARAMS_TABLE_INST(alpaqa::PANTRParams<alpaqa::EigenConfigd>);
extern PARAMS_TABLE_INST(alpaqa::PANTRParams<alpaqa::EigenConfigl>);
#ifdef ALPAQA_WITH_QUAD_PRECISION
extern PARAMS_TABLE_INST(alpaqa::PANTRParams<alpaqa::EigenConfigq>);
#endif
