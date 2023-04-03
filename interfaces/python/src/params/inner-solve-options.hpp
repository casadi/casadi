#pragma once

#include <alpaqa/inner/inner-solve-options.hpp>
#include <dict/kwargs-to-struct.hpp>

PARAMS_TABLE_DECL(alpaqa::InnerSolveOptions<Conf>);

extern PARAMS_TABLE_INST(alpaqa::InnerSolveOptions<alpaqa::EigenConfigf>);
extern PARAMS_TABLE_INST(alpaqa::InnerSolveOptions<alpaqa::EigenConfigd>);
extern PARAMS_TABLE_INST(alpaqa::InnerSolveOptions<alpaqa::EigenConfigl>);
#ifdef ALPAQA_WITH_QUAD_PRECISION
extern PARAMS_TABLE_INST(alpaqa::InnerSolveOptions<alpaqa::EigenConfigq>);
#endif
