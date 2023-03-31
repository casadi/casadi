#include "lbfgs-params.hpp"

PARAMS_TABLE_DEF(alpaqa::LBFGSParams<Conf>,    //
                 PARAMS_MEMBER(memory),        //
                 PARAMS_MEMBER(min_div_fac),   //
                 PARAMS_MEMBER(min_abs_s),     //
                 PARAMS_MEMBER(cbfgs),         //
                 PARAMS_MEMBER(force_pos_def), //
                 PARAMS_MEMBER(stepsize),      //
);

PARAMS_TABLE_DEF(alpaqa::CBFGSParams<Conf>, //
                 PARAMS_MEMBER(α),          //
                 PARAMS_MEMBER(ϵ),          //
);

PARAMS_TABLE_INST(alpaqa::LBFGSParams<alpaqa::EigenConfigf>);
PARAMS_TABLE_INST(alpaqa::LBFGSParams<alpaqa::EigenConfigd>);
PARAMS_TABLE_INST(alpaqa::LBFGSParams<alpaqa::EigenConfigl>);
#ifdef ALPAQA_WITH_QUAD_PRECISION
PARAMS_TABLE_INST(alpaqa::LBFGSParams<alpaqa::EigenConfigq>);
#endif

PARAMS_TABLE_INST(alpaqa::CBFGSParams<alpaqa::EigenConfigf>);
PARAMS_TABLE_INST(alpaqa::CBFGSParams<alpaqa::EigenConfigd>);
PARAMS_TABLE_INST(alpaqa::CBFGSParams<alpaqa::EigenConfigl>);
#ifdef ALPAQA_WITH_QUAD_PRECISION
PARAMS_TABLE_INST(alpaqa::CBFGSParams<alpaqa::EigenConfigq>);
#endif