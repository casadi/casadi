#pragma once

#include <alpaqa/inner/panoc-ocp.hpp>
#include <kwargs-to-struct.hpp>

PARAMS_TABLE_DECL(alpaqa::PANOCOCPParams<Conf>);

extern PARAMS_TABLE_INST(alpaqa::PANOCOCPParams<alpaqa::EigenConfigf>);
extern PARAMS_TABLE_INST(alpaqa::PANOCOCPParams<alpaqa::EigenConfigd>);
extern PARAMS_TABLE_INST(alpaqa::PANOCOCPParams<alpaqa::EigenConfigl>);
#ifdef ALPAQA_WITH_QUAD_PRECISION
extern PARAMS_TABLE_INST(alpaqa::PANOCOCPParams<alpaqa::EigenConfigq>);
#endif
