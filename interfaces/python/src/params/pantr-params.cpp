#include "pantr-params.hpp"

PARAMS_TABLE_DEF(alpaqa::PANTRParams<Conf>,                                     //
                 PARAMS_MEMBER(Lipschitz),                                      //
                 PARAMS_MEMBER(max_iter),                                       //
                 PARAMS_MEMBER(max_time),                                       //
                 PARAMS_MEMBER(L_min),                                          //
                 PARAMS_MEMBER(L_max),                                          //
                 PARAMS_MEMBER(stop_crit),                                      //
                 PARAMS_MEMBER(max_no_progress),                                //
                 PARAMS_MEMBER(print_interval),                                 //
                 PARAMS_MEMBER(print_precision),                                //
                 PARAMS_MEMBER(quadratic_upperbound_tolerance_factor),          //
                 PARAMS_MEMBER(TR_tolerance_factor),                            //
                 PARAMS_MEMBER(μ1),                                             //
                 PARAMS_MEMBER(μ2),                                             //
                 PARAMS_MEMBER(c1),                                             //
                 PARAMS_MEMBER(c2),                                             //
                 PARAMS_MEMBER(c3),                                             //
                 PARAMS_MEMBER(Δ_0),                                            //
                 PARAMS_MEMBER(Δ_min),                                          //
                 PARAMS_MEMBER(compute_ratio_using_new_γ),                      //
                 PARAMS_MEMBER(update_direction_on_prox_step),                  //
                 PARAMS_MEMBER(recompute_last_prox_step_after_direction_reset), //
                 PARAMS_MEMBER(disable_acceleration),                           //
);

PARAMS_TABLE_INST(alpaqa::PANTRParams<alpaqa::EigenConfigf>);
PARAMS_TABLE_INST(alpaqa::PANTRParams<alpaqa::EigenConfigd>);
PARAMS_TABLE_INST(alpaqa::PANTRParams<alpaqa::EigenConfigl>);
#ifdef ALPAQA_WITH_QUAD_PRECISION
PARAMS_TABLE_INST(alpaqa::PANTRParams<alpaqa::EigenConfigq>);
#endif
