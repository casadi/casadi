#pragma once

#include <alpaqa/inner/panoc.hpp>
#include <kwargs-to-struct.hpp>

PARAMS_TABLE_DECL(alpaqa::PANOCParams<Conf>);
PARAMS_TABLE_DECL(alpaqa::LipschitzEstimateParams<Conf>);

extern PARAMS_TABLE_INST(alpaqa::PANOCParams<alpaqa::EigenConfigf>);
extern PARAMS_TABLE_INST(alpaqa::PANOCParams<alpaqa::EigenConfigd>);
extern PARAMS_TABLE_INST(alpaqa::PANOCParams<alpaqa::EigenConfigl>);
#ifdef ALPAQA_WITH_QUAD_PRECISION
extern PARAMS_TABLE_INST(alpaqa::PANOCParams<alpaqa::EigenConfigq>);
#endif

extern PARAMS_TABLE_INST(alpaqa::LipschitzEstimateParams<alpaqa::EigenConfigf>);
extern PARAMS_TABLE_INST(alpaqa::LipschitzEstimateParams<alpaqa::EigenConfigd>);
extern PARAMS_TABLE_INST(alpaqa::LipschitzEstimateParams<alpaqa::EigenConfigl>);
#ifdef ALPAQA_WITH_QUAD_PRECISION
extern PARAMS_TABLE_INST(alpaqa::LipschitzEstimateParams<alpaqa::EigenConfigq>);
#endif
