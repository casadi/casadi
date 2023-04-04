#include "params.hpp"

PARAMS_TABLE_DEF(alpaqa::ALMParams<Conf>,                       //
                 PARAMS_MEMBER(tolerance),                      //
                 PARAMS_MEMBER(dual_tolerance),                 //
                 PARAMS_MEMBER(penalty_update_factor),          //
                 PARAMS_MEMBER(penalty_update_factor_lower),    //
                 PARAMS_MEMBER(min_penalty_update_factor),      //
                 PARAMS_MEMBER(initial_penalty),                //
                 PARAMS_MEMBER(initial_penalty_factor),         //
                 PARAMS_MEMBER(initial_penalty_lower),          //
                 PARAMS_MEMBER(initial_tolerance),              //
                 PARAMS_MEMBER(initial_tolerance_increase),     //
                 PARAMS_MEMBER(tolerance_update_factor),        //
                 PARAMS_MEMBER(ρ_increase),                     //
                 PARAMS_MEMBER(ρ_max),                          //
                 PARAMS_MEMBER(rel_penalty_increase_threshold), //
                 PARAMS_MEMBER(max_multiplier),                 //
                 PARAMS_MEMBER(max_penalty),                    //
                 PARAMS_MEMBER(min_penalty),                    //
                 PARAMS_MEMBER(max_iter),                       //
                 PARAMS_MEMBER(max_time),                       //
                 PARAMS_MEMBER(max_num_initial_retries),        //
                 PARAMS_MEMBER(max_num_retries),                //
                 PARAMS_MEMBER(max_total_num_retries),          //
                 PARAMS_MEMBER(print_interval),                 //
                 PARAMS_MEMBER(single_penalty_factor),          //
);

PARAMS_TABLE_INST(alpaqa::ALMParams<alpaqa::EigenConfigf>);
PARAMS_TABLE_INST(alpaqa::ALMParams<alpaqa::EigenConfigd>);
PARAMS_TABLE_INST(alpaqa::ALMParams<alpaqa::EigenConfigl>);
#ifdef ALPAQA_WITH_QUAD_PRECISION
PARAMS_TABLE_INST(alpaqa::ALMParams<alpaqa::EigenConfigq>);
#endif
