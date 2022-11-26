#include "panoc-ocp-params.hpp"

template <alpaqa::Config Conf>
const dict_to_struct_table_t<alpaqa::PANOCOCPParams<Conf>>
    dict_to_struct_table<alpaqa::PANOCOCPParams<Conf>>::table{
        // clang-format off
        {"Lipschitz", &alpaqa::PANOCOCPParams<Conf>::Lipschitz},
        {"max_iter", &alpaqa::PANOCOCPParams<Conf>::max_iter},
        {"max_time", &alpaqa::PANOCOCPParams<Conf>::max_time},
        {"τ_min", &alpaqa::PANOCOCPParams<Conf>::τ_min},
        {"β", &alpaqa::PANOCOCPParams<Conf>::β},
        {"L_min", &alpaqa::PANOCOCPParams<Conf>::L_min},
        {"L_max", &alpaqa::PANOCOCPParams<Conf>::L_max},
        {"L_max_inc", &alpaqa::PANOCOCPParams<Conf>::L_max_inc},
        {"stop_crit", &alpaqa::PANOCOCPParams<Conf>::stop_crit},
        {"max_no_progress", &alpaqa::PANOCOCPParams<Conf>::max_no_progress},
        {"gn_interval", &alpaqa::PANOCOCPParams<Conf>::gn_interval},
        {"gn_sticky", &alpaqa::PANOCOCPParams<Conf>::gn_sticky},
        {"reset_lbfgs_on_gn_step", &alpaqa::PANOCOCPParams<Conf>::reset_lbfgs_on_gn_step},
        {"lqr_factor_cholesky", &alpaqa::PANOCOCPParams<Conf>::lqr_factor_cholesky},
        {"print_interval", &alpaqa::PANOCOCPParams<Conf>::print_interval},
        {"print_precision", &alpaqa::PANOCOCPParams<Conf>::print_precision},
        {"quadratic_upperbound_tolerance_factor", &alpaqa::PANOCOCPParams<Conf>::quadratic_upperbound_tolerance_factor},
        {"linesearch_tolerance_factor", &alpaqa::PANOCOCPParams<Conf>::linesearch_tolerance_factor},
        {"disable_acceleration", &alpaqa::PANOCOCPParams<Conf>::disable_acceleration},
        // clang-format on
    };

template struct dict_to_struct_table<alpaqa::PANOCOCPParams<alpaqa::EigenConfigf>>;
template struct dict_to_struct_table<alpaqa::PANOCOCPParams<alpaqa::EigenConfigd>>;
template struct dict_to_struct_table<alpaqa::PANOCOCPParams<alpaqa::EigenConfigl>>;
#ifdef ALPAQA_WITH_QUAD_PRECISION
template struct dict_to_struct_table<alpaqa::PANOCOCPParams<alpaqa::EigenConfigf>>;
#endif
