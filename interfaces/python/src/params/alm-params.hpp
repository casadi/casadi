#pragma once

#include <alpaqa/outer/alm.hpp>

#include <dict/kwargs-to-struct.hpp>

PARAMS_TABLE_DECL(alpaqa::ALMParams<Conf>);

extern PARAMS_TABLE_INST(alpaqa::ALMParams<alpaqa::EigenConfigf>);
extern PARAMS_TABLE_INST(alpaqa::ALMParams<alpaqa::EigenConfigd>);
extern PARAMS_TABLE_INST(alpaqa::ALMParams<alpaqa::EigenConfigl>);
#ifdef ALPAQA_WITH_QUAD_PRECISION
extern PARAMS_TABLE_INST(alpaqa::ALMParams<alpaqa::EigenConfigq>);
#endif
