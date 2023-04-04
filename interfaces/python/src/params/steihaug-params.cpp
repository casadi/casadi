#include "steihaug-params.hpp"

PARAMS_TABLE_DEF(alpaqa::SteihaugCGParams<Conf>, //
                 PARAMS_MEMBER(tol_scale),       //
                 PARAMS_MEMBER(tol_scale_root),  //
                 PARAMS_MEMBER(tol_max),         //
                 PARAMS_MEMBER(max_iter_factor), //
);

PARAMS_TABLE_INST(alpaqa::SteihaugCGParams<alpaqa::EigenConfigf>);
PARAMS_TABLE_INST(alpaqa::SteihaugCGParams<alpaqa::EigenConfigd>);
PARAMS_TABLE_INST(alpaqa::SteihaugCGParams<alpaqa::EigenConfigl>);
#ifdef ALPAQA_WITH_QUAD_PRECISION
PARAMS_TABLE_INST(alpaqa::SteihaugCGParams<alpaqa::EigenConfigq>);
#endif
