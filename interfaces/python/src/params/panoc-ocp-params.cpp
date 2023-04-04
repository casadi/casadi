#include "panoc-ocp-params.hpp"
#include "lbfgs-params.hpp"

PARAMS_TABLE_DEF(alpaqa::PANOCOCPParams<Conf>,
                 PARAMS_MEMBER(Lipschitz),                             //
                 PARAMS_MEMBER(max_iter),                              //
                 PARAMS_MEMBER(max_time),                              //
                 PARAMS_MEMBER(min_linesearch_coefficient),            //
                 PARAMS_MEMBER(linesearch_strictness_factor),          //
                 PARAMS_MEMBER(L_min),                                 //
                 PARAMS_MEMBER(L_max),                                 //
                 PARAMS_MEMBER(L_max_inc),                             //
                 PARAMS_MEMBER(stop_crit),                             //
                 PARAMS_MEMBER(max_no_progress),                       //
                 PARAMS_MEMBER(gn_interval),                           //
                 PARAMS_MEMBER(gn_sticky),                             //
                 PARAMS_MEMBER(reset_lbfgs_on_gn_step),                //
                 PARAMS_MEMBER(lqr_factor_cholesky),                   //
                 PARAMS_MEMBER(lbfgs_params),                          //
                 PARAMS_MEMBER(print_interval),                        //
                 PARAMS_MEMBER(print_precision),                       //
                 PARAMS_MEMBER(quadratic_upperbound_tolerance_factor), //
                 PARAMS_MEMBER(linesearch_tolerance_factor),           //
                 PARAMS_MEMBER(disable_acceleration),                  //
);

PARAMS_TABLE_INST(alpaqa::PANOCOCPParams<alpaqa::EigenConfigf>);
PARAMS_TABLE_INST(alpaqa::PANOCOCPParams<alpaqa::EigenConfigd>);
PARAMS_TABLE_INST(alpaqa::PANOCOCPParams<alpaqa::EigenConfigl>);
#ifdef ALPAQA_WITH_QUAD_PRECISION
PARAMS_TABLE_INST(alpaqa::PANOCOCPParams<alpaqa::EigenConfigq>);
#endif
