#pragma once

#include <alpaqa/accelerators/lbfgs.hpp>
#include <kwargs-to-struct.hpp>

PARAMS_TABLE_DECL(alpaqa::LBFGSParams<Conf>);
PARAMS_TABLE_DECL(alpaqa::CBFGSParams<Conf>);

extern PARAMS_TABLE_INST(alpaqa::LBFGSParams<alpaqa::EigenConfigf>);
extern PARAMS_TABLE_INST(alpaqa::LBFGSParams<alpaqa::EigenConfigd>);
extern PARAMS_TABLE_INST(alpaqa::LBFGSParams<alpaqa::EigenConfigl>);
#ifdef ALPAQA_WITH_QUAD_PRECISION
extern PARAMS_TABLE_INST(alpaqa::LBFGSParams<alpaqa::EigenConfigq>);
#endif

extern PARAMS_TABLE_INST(alpaqa::CBFGSParams<alpaqa::EigenConfigf>);
extern PARAMS_TABLE_INST(alpaqa::CBFGSParams<alpaqa::EigenConfigd>);
extern PARAMS_TABLE_INST(alpaqa::CBFGSParams<alpaqa::EigenConfigl>);
#ifdef ALPAQA_WITH_QUAD_PRECISION
extern PARAMS_TABLE_INST(alpaqa::CBFGSParams<alpaqa::EigenConfigq>);
#endif
