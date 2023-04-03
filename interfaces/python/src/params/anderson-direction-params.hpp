#pragma once

#include <alpaqa/inner/directions/panoc/anderson.hpp>
#include <dict/kwargs-to-struct.hpp>

PARAMS_TABLE_DECL(alpaqa::AndersonDirectionParams<Conf>);

extern PARAMS_TABLE_INST(alpaqa::AndersonDirectionParams<alpaqa::EigenConfigf>);
extern PARAMS_TABLE_INST(alpaqa::AndersonDirectionParams<alpaqa::EigenConfigd>);
extern PARAMS_TABLE_INST(alpaqa::AndersonDirectionParams<alpaqa::EigenConfigl>);
#ifdef ALPAQA_WITH_QUAD_PRECISION
extern PARAMS_TABLE_INST(alpaqa::AndersonDirectionParams<alpaqa::EigenConfigq>);
#endif
