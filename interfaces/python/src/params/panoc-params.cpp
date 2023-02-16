#include "panoc-params.hpp"

template <alpaqa::Config Conf>
const dict_to_struct_table_t<alpaqa::PANOCParams<Conf>>
    dict_to_struct_table<alpaqa::PANOCParams<Conf>>::table{
        // clang-format off
        {"Lipschitz", &alpaqa::PANOCParams<Conf>::Lipschitz},
        {"max_iter", &alpaqa::PANOCParams<Conf>::max_iter},
        {"max_time", &alpaqa::PANOCParams<Conf>::max_time},
        {"τ_min", &alpaqa::PANOCParams<Conf>::τ_min},
        {"force_linesearch", &alpaqa::PANOCParams<Conf>::force_linesearch},
        {"β", &alpaqa::PANOCParams<Conf>::β},
        {"L_min", &alpaqa::PANOCParams<Conf>::L_min},
        {"L_max", &alpaqa::PANOCParams<Conf>::L_max},
        {"stop_crit", &alpaqa::PANOCParams<Conf>::stop_crit},
        {"max_no_progress", &alpaqa::PANOCParams<Conf>::max_no_progress},
        {"print_interval", &alpaqa::PANOCParams<Conf>::print_interval},
        {"print_precision", &alpaqa::PANOCParams<Conf>::print_precision},
        {"quadratic_upperbound_tolerance_factor", &alpaqa::PANOCParams<Conf>::quadratic_upperbound_tolerance_factor},
        {"linesearch_tolerance_factor", &alpaqa::PANOCParams<Conf>::linesearch_tolerance_factor},
        // clang-format on
    };

template <alpaqa::Config Conf>
const dict_to_struct_table_t<alpaqa::LipschitzEstimateParams<Conf>>
    dict_to_struct_table<alpaqa::LipschitzEstimateParams<Conf>>::table{
        {"L_0", &alpaqa::LipschitzEstimateParams<Conf>::L_0},
        {"δ", &alpaqa::LipschitzEstimateParams<Conf>::δ},
        {"ε", &alpaqa::LipschitzEstimateParams<Conf>::ε},
        {"Lγ_factor", &alpaqa::LipschitzEstimateParams<Conf>::Lγ_factor},
    };

template struct dict_to_struct_table<alpaqa::PANOCParams<alpaqa::EigenConfigf>>;
template struct dict_to_struct_table<alpaqa::PANOCParams<alpaqa::EigenConfigd>>;
template struct dict_to_struct_table<alpaqa::PANOCParams<alpaqa::EigenConfigl>>;
#ifdef ALPAQA_WITH_QUAD_PRECISION
template struct dict_to_struct_table<alpaqa::PANOCParams<alpaqa::EigenConfigq>>;
#endif

template struct dict_to_struct_table<alpaqa::LipschitzEstimateParams<alpaqa::EigenConfigf>>;
template struct dict_to_struct_table<alpaqa::LipschitzEstimateParams<alpaqa::EigenConfigd>>;
template struct dict_to_struct_table<alpaqa::LipschitzEstimateParams<alpaqa::EigenConfigl>>;
#ifdef ALPAQA_WITH_QUAD_PRECISION
template struct dict_to_struct_table<alpaqa::LipschitzEstimateParams<alpaqa::EigenConfigq>>;
#endif
