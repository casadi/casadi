#include "params.hpp"

PARAMS_TABLE_DEF(alpaqa::InnerSolveOptions<Conf>,
                 PARAMS_MEMBER(always_overwrite_results), //
                 PARAMS_MEMBER(max_time),                 //
                 PARAMS_MEMBER(tolerance),                //
);

PARAMS_TABLE_INST(alpaqa::InnerSolveOptions<alpaqa::EigenConfigf>);
PARAMS_TABLE_INST(alpaqa::InnerSolveOptions<alpaqa::EigenConfigd>);
PARAMS_TABLE_INST(alpaqa::InnerSolveOptions<alpaqa::EigenConfigl>);
#ifdef ALPAQA_WITH_QUAD_PRECISION
PARAMS_TABLE_INST(alpaqa::InnerSolveOptions<alpaqa::EigenConfigq>);
#endif
