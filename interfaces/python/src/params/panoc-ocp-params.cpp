#include "panoc-ocp-params.hpp"

template <alpaqa::Config Conf>
const dict_to_struct_table_t<alpaqa::experimental::PANOCOCPParams<Conf>>
    dict_to_struct_table<alpaqa::experimental::PANOCOCPParams<Conf>>::table{
        // clang-format off
        {"Lipschitz", &alpaqa::experimental::PANOCOCPParams<Conf>::Lipschitz},
        {"max_iter", &alpaqa::experimental::PANOCOCPParams<Conf>::max_iter},
        {"max_time", &alpaqa::experimental::PANOCOCPParams<Conf>::max_time},
        {"τ_min", &alpaqa::experimental::PANOCOCPParams<Conf>::τ_min},
        {"β", &alpaqa::experimental::PANOCOCPParams<Conf>::β},
        {"L_min", &alpaqa::experimental::PANOCOCPParams<Conf>::L_min},
        {"L_max", &alpaqa::experimental::PANOCOCPParams<Conf>::L_max},
        {"L_max_inc", &alpaqa::experimental::PANOCOCPParams<Conf>::L_max_inc},
        {"stop_crit", &alpaqa::experimental::PANOCOCPParams<Conf>::stop_crit},
        {"max_no_progress", &alpaqa::experimental::PANOCOCPParams<Conf>::max_no_progress},
        {"gn_interval", &alpaqa::experimental::PANOCOCPParams<Conf>::gn_interval},
        {"gn_sticky", &alpaqa::experimental::PANOCOCPParams<Conf>::gn_sticky},
        {"reset_lbfgs_on_gn_step", &alpaqa::experimental::PANOCOCPParams<Conf>::reset_lbfgs_on_gn_step},
        {"lqr_factor_cholesky", &alpaqa::experimental::PANOCOCPParams<Conf>::lqr_factor_cholesky},
        {"print_interval", &alpaqa::experimental::PANOCOCPParams<Conf>::print_interval},
        {"print_precision", &alpaqa::experimental::PANOCOCPParams<Conf>::print_precision},
        {"quadratic_upperbound_tolerance_factor", &alpaqa::experimental::PANOCOCPParams<Conf>::quadratic_upperbound_tolerance_factor},
        {"linesearch_tolerance_factor", &alpaqa::experimental::PANOCOCPParams<Conf>::linesearch_tolerance_factor},
        {"disable_acceleration", &alpaqa::experimental::PANOCOCPParams<Conf>::disable_acceleration},
        // clang-format on
    };

template struct dict_to_struct_table<alpaqa::experimental::PANOCOCPParams<alpaqa::EigenConfigf>>;
template struct dict_to_struct_table<alpaqa::experimental::PANOCOCPParams<alpaqa::EigenConfigd>>;
template struct dict_to_struct_table<alpaqa::experimental::PANOCOCPParams<alpaqa::EigenConfigl>>;
#ifdef ALPAQA_WITH_QUAD_PRECISION
template struct dict_to_struct_table<alpaqa::experimental::PANOCOCPParams<alpaqa::EigenConfigq>>;
#endif
