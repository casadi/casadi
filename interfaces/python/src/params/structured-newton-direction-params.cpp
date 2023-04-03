#include "structured-newton-direction-params.hpp"

PARAMS_TABLE_DEF(alpaqa::StructuredNewtonRegularizationParams<Conf>, //
                 PARAMS_MEMBER(min_eig),                             //
                 PARAMS_MEMBER(print_eig),                           //
);
PARAMS_TABLE_DEF(alpaqa::StructuredNewtonDirectionParams<Conf>,      //
                 PARAMS_MEMBER(hessian_vec),                         //
);

PARAMS_TABLE_INST(alpaqa::StructuredNewtonRegularizationParams<alpaqa::EigenConfigf>);
PARAMS_TABLE_INST(alpaqa::StructuredNewtonRegularizationParams<alpaqa::EigenConfigd>);
PARAMS_TABLE_INST(alpaqa::StructuredNewtonRegularizationParams<alpaqa::EigenConfigl>);
#ifdef ALPAQA_WITH_QUAD_PRECISION
PARAMS_TABLE_INST(alpaqa::StructuredNewtonRegularizationParams<alpaqa::EigenConfigq>);
#endif

PARAMS_TABLE_INST(alpaqa::StructuredNewtonDirectionParams<alpaqa::EigenConfigf>);
PARAMS_TABLE_INST(alpaqa::StructuredNewtonDirectionParams<alpaqa::EigenConfigd>);
PARAMS_TABLE_INST(alpaqa::StructuredNewtonDirectionParams<alpaqa::EigenConfigl>);
#ifdef ALPAQA_WITH_QUAD_PRECISION
PARAMS_TABLE_INST(alpaqa::StructuredNewtonDirectionParams<alpaqa::EigenConfigq>);
#endif
