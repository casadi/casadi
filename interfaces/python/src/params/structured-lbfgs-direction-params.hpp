#pragma once

#include <alpaqa/inner/directions/panoc/structured-lbfgs.hpp>
#include <kwargs-to-struct.hpp>

PARAMS_TABLE_DECL(alpaqa::StructuredLBFGSDirectionParams<Conf>);

extern PARAMS_TABLE_INST(alpaqa::StructuredLBFGSDirectionParams<alpaqa::EigenConfigf>);
extern PARAMS_TABLE_INST(alpaqa::StructuredLBFGSDirectionParams<alpaqa::EigenConfigd>);
extern PARAMS_TABLE_INST(alpaqa::StructuredLBFGSDirectionParams<alpaqa::EigenConfigl>);
#ifdef ALPAQA_WITH_QUAD_PRECISION
extern PARAMS_TABLE_INST(alpaqa::StructuredLBFGSDirectionParams<alpaqa::EigenConfigq>);
#endif
