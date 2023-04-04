#pragma once

#include <alpaqa/inner/zerofpr.hpp>
#include <params/params.hpp>

PARAMS_TABLE_DECL(alpaqa::ZeroFPRParams<Conf>);

extern PARAMS_TABLE_INST(alpaqa::ZeroFPRParams<alpaqa::EigenConfigf>);
extern PARAMS_TABLE_INST(alpaqa::ZeroFPRParams<alpaqa::EigenConfigd>);
extern PARAMS_TABLE_INST(alpaqa::ZeroFPRParams<alpaqa::EigenConfigl>);
#ifdef ALPAQA_WITH_QUAD_PRECISION
extern PARAMS_TABLE_INST(alpaqa::ZeroFPRParams<alpaqa::EigenConfigq>);
#endif
