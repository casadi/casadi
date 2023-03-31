#include "params.hpp"

PARAMS_TABLE_DEF(alpaqa::ALMParams<Conf>,                //
                 PARAMS_MEMBER(ε),                       //
                 PARAMS_MEMBER(δ),                       //
                 PARAMS_MEMBER(Δ),                       //
                 PARAMS_MEMBER(Δ_lower),                 //
                 PARAMS_MEMBER(Δ_min),                   //
                 PARAMS_MEMBER(Σ_0),                     //
                 PARAMS_MEMBER(σ_0),                     //
                 PARAMS_MEMBER(Σ_0_lower),               //
                 PARAMS_MEMBER(ε_0),                     //
                 PARAMS_MEMBER(ε_0_increase),            //
                 PARAMS_MEMBER(ρ),                       //
                 PARAMS_MEMBER(ρ_increase),              //
                 PARAMS_MEMBER(ρ_max),                   //
                 PARAMS_MEMBER(θ),                       //
                 PARAMS_MEMBER(M),                       //
                 PARAMS_MEMBER(Σ_max),                   //
                 PARAMS_MEMBER(Σ_min),                   //
                 PARAMS_MEMBER(max_iter),                //
                 PARAMS_MEMBER(max_time),                //
                 PARAMS_MEMBER(max_num_initial_retries), //
                 PARAMS_MEMBER(max_num_retries),         //
                 PARAMS_MEMBER(max_total_num_retries),   //
                 PARAMS_MEMBER(print_interval),          //
                 PARAMS_MEMBER(single_penalty_factor),   //
);

PARAMS_TABLE_INST(alpaqa::ALMParams<alpaqa::EigenConfigf>);
PARAMS_TABLE_INST(alpaqa::ALMParams<alpaqa::EigenConfigd>);
PARAMS_TABLE_INST(alpaqa::ALMParams<alpaqa::EigenConfigl>);
#ifdef ALPAQA_WITH_QUAD_PRECISION
PARAMS_TABLE_INST(alpaqa::ALMParams<alpaqa::EigenConfigq>);
#endif
