#pragma once

#include <alpaqa/inner/directions/pantr/newton-tr.hpp>
#include <params/params.hpp>

PARAMS_TABLE_DECL(alpaqa::NewtonTRDirectionParams<Conf>);

extern PARAMS_TABLE_INST(alpaqa::NewtonTRDirectionParams<alpaqa::EigenConfigf>);
extern PARAMS_TABLE_INST(alpaqa::NewtonTRDirectionParams<alpaqa::EigenConfigd>);
extern PARAMS_TABLE_INST(alpaqa::NewtonTRDirectionParams<alpaqa::EigenConfigl>);
#ifdef ALPAQA_WITH_QUAD_PRECISION
extern PARAMS_TABLE_INST(alpaqa::NewtonTRDirectionParams<alpaqa::EigenConfigq>);
#endif
