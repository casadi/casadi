#include "zerofpr-params.hpp"

PARAMS_TABLE_DEF(alpaqa::ZeroFPRParams<Conf>,                               //
                 PARAMS_MEMBER(Lipschitz),                                  //
                 PARAMS_MEMBER(max_iter),                                   //
                 PARAMS_MEMBER(max_time),                                   //
                 PARAMS_MEMBER(min_linesearch_coefficient),                 //
                 PARAMS_MEMBER(force_linesearch),                           //
                 PARAMS_MEMBER(linesearch_strictness_factor),               //
                 PARAMS_MEMBER(L_min),                                      //
                 PARAMS_MEMBER(L_max),                                      //
                 PARAMS_MEMBER(stop_crit),                                  //
                 PARAMS_MEMBER(max_no_progress),                            //
                 PARAMS_MEMBER(print_interval),                             //
                 PARAMS_MEMBER(print_precision),                            //
                 PARAMS_MEMBER(quadratic_upperbound_tolerance_factor),      //
                 PARAMS_MEMBER(linesearch_tolerance_factor),                //
                 PARAMS_MEMBER(update_direction_in_candidate),              //
                 PARAMS_MEMBER(recompute_last_prox_step_after_lbfgs_flush), //
                 PARAMS_MEMBER(update_direction_from_prox_step),            //
);

PARAMS_TABLE_INST(alpaqa::ZeroFPRParams<alpaqa::EigenConfigf>);
PARAMS_TABLE_INST(alpaqa::ZeroFPRParams<alpaqa::EigenConfigd>);
PARAMS_TABLE_INST(alpaqa::ZeroFPRParams<alpaqa::EigenConfigl>);
#ifdef ALPAQA_WITH_QUAD_PRECISION
PARAMS_TABLE_INST(alpaqa::ZeroFPRParams<alpaqa::EigenConfigq>);
#endif
