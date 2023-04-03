#pragma once

#include <alpaqa/inner/directions/panoc/lbfgs.hpp>
#include <dict/kwargs-to-struct.hpp>

PARAMS_TABLE_DECL(alpaqa::LBFGSDirectionParams<Conf>);

extern PARAMS_TABLE_INST(alpaqa::LBFGSDirectionParams<alpaqa::EigenConfigf>);
extern PARAMS_TABLE_INST(alpaqa::LBFGSDirectionParams<alpaqa::EigenConfigd>);
extern PARAMS_TABLE_INST(alpaqa::LBFGSDirectionParams<alpaqa::EigenConfigl>);
#ifdef ALPAQA_WITH_QUAD_PRECISION
extern PARAMS_TABLE_INST(alpaqa::LBFGSDirectionParams<alpaqa::EigenConfigq>);
#endif
