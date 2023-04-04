#include "newton-tr-direction-params.hpp"

PARAMS_TABLE_DEF(alpaqa::NewtonTRDirectionParams<Conf>, //
                 PARAMS_MEMBER(rescale_when_γ_changes), //
                 PARAMS_MEMBER(hessian_vec_factor),     //
                 PARAMS_MEMBER(finite_diff),            //
                 PARAMS_MEMBER(finite_diff_stepsize),   //
);

PARAMS_TABLE_INST(alpaqa::NewtonTRDirectionParams<alpaqa::EigenConfigf>);
PARAMS_TABLE_INST(alpaqa::NewtonTRDirectionParams<alpaqa::EigenConfigd>);
PARAMS_TABLE_INST(alpaqa::NewtonTRDirectionParams<alpaqa::EigenConfigl>);
#ifdef ALPAQA_WITH_QUAD_PRECISION
PARAMS_TABLE_INST(alpaqa::NewtonTRDirectionParams<alpaqa::EigenConfigq>);
#endif
