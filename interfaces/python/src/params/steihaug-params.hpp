#pragma once

#include <alpaqa/accelerators/steihaugcg.hpp>
#include <params/params.hpp>

PARAMS_TABLE_DECL(alpaqa::SteihaugCGParams<Conf>);

extern PARAMS_TABLE_INST(alpaqa::SteihaugCGParams<alpaqa::EigenConfigf>);
extern PARAMS_TABLE_INST(alpaqa::SteihaugCGParams<alpaqa::EigenConfigd>);
extern PARAMS_TABLE_INST(alpaqa::SteihaugCGParams<alpaqa::EigenConfigl>);
#ifdef ALPAQA_WITH_QUAD_PRECISION
extern PARAMS_TABLE_INST(alpaqa::SteihaugCGParams<alpaqa::EigenConfigq>);
#endif
