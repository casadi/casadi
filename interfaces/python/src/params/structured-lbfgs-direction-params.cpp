#include "structured-lbfgs-direction-params.hpp"

PARAMS_TABLE_DEF(alpaqa::StructuredLBFGSDirectionParams<Conf>,  //
                 PARAMS_MEMBER(hessian_vec),                    //
                 PARAMS_MEMBER(hessian_vec_finite_differences), //
                 PARAMS_MEMBER(full_augmented_hessian),         //
);

PARAMS_TABLE_INST(alpaqa::StructuredLBFGSDirectionParams<alpaqa::EigenConfigf>);
PARAMS_TABLE_INST(alpaqa::StructuredLBFGSDirectionParams<alpaqa::EigenConfigd>);
PARAMS_TABLE_INST(alpaqa::StructuredLBFGSDirectionParams<alpaqa::EigenConfigl>);
#ifdef ALPAQA_WITH_QUAD_PRECISION
PARAMS_TABLE_INST(alpaqa::StructuredLBFGSDirectionParams<alpaqa::EigenConfigq>);
#endif
