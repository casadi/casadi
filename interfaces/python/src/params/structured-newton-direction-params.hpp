#pragma once

#include <alpaqa/inner/directions/panoc/structured-newton.hpp>
#include <dict/kwargs-to-struct.hpp>

PARAMS_TABLE_DECL(alpaqa::StructuredNewtonRegularizationParams<Conf>);
PARAMS_TABLE_DECL(alpaqa::StructuredNewtonDirectionParams<Conf>);

extern PARAMS_TABLE_INST(alpaqa::StructuredNewtonRegularizationParams<alpaqa::EigenConfigf>);
extern PARAMS_TABLE_INST(alpaqa::StructuredNewtonRegularizationParams<alpaqa::EigenConfigd>);
extern PARAMS_TABLE_INST(alpaqa::StructuredNewtonRegularizationParams<alpaqa::EigenConfigl>);
#ifdef ALPAQA_WITH_QUAD_PRECISION
extern PARAMS_TABLE_INST(alpaqa::StructuredNewtonRegularizationParams<alpaqa::EigenConfigq>);
#endif

extern PARAMS_TABLE_INST(alpaqa::StructuredNewtonDirectionParams<alpaqa::EigenConfigf>);
extern PARAMS_TABLE_INST(alpaqa::StructuredNewtonDirectionParams<alpaqa::EigenConfigd>);
extern PARAMS_TABLE_INST(alpaqa::StructuredNewtonDirectionParams<alpaqa::EigenConfigl>);
#ifdef ALPAQA_WITH_QUAD_PRECISION
extern PARAMS_TABLE_INST(alpaqa::StructuredNewtonDirectionParams<alpaqa::EigenConfigq>);
#endif
