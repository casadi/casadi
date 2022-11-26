#include <params/lbfgs-params.hpp>

template <alpaqa::Config Conf>
const dict_to_struct_table_t<alpaqa::LBFGSParams<Conf>>
    dict_to_struct_table<alpaqa::LBFGSParams<Conf>>::table{
        {"memory", &alpaqa::LBFGSParams<Conf>::memory},
        {"min_div_fac", &alpaqa::LBFGSParams<Conf>::min_div_fac},
        {"min_abs_s", &alpaqa::LBFGSParams<Conf>::min_abs_s},
        {"cbfgs", &alpaqa::LBFGSParams<Conf>::cbfgs},
        {"force_pos_def", &alpaqa::LBFGSParams<Conf>::force_pos_def},
        {"stepsize", &alpaqa::LBFGSParams<Conf>::stepsize},
    };

template <alpaqa::Config Conf>
const dict_to_struct_table_t<alpaqa::CBFGSParams<Conf>>
    dict_to_struct_table<alpaqa::CBFGSParams<Conf>>::table{
        {"α", &alpaqa::CBFGSParams<Conf>::α},
        {"ϵ", &alpaqa::CBFGSParams<Conf>::ϵ},
    };

template struct dict_to_struct_table<alpaqa::LBFGSParams<alpaqa::EigenConfigf>>;
template struct dict_to_struct_table<alpaqa::LBFGSParams<alpaqa::EigenConfigd>>;
template struct dict_to_struct_table<alpaqa::LBFGSParams<alpaqa::EigenConfigl>>;
#ifdef ALPAQA_WITH_QUAD_PRECISION
template struct dict_to_struct_table<alpaqa::LBFGSParams<alpaqa::EigenConfigf>>;
#endif

template struct dict_to_struct_table<alpaqa::CBFGSParams<alpaqa::EigenConfigf>>;
template struct dict_to_struct_table<alpaqa::CBFGSParams<alpaqa::EigenConfigd>>;
template struct dict_to_struct_table<alpaqa::CBFGSParams<alpaqa::EigenConfigl>>;
#ifdef ALPAQA_WITH_QUAD_PRECISION
template struct dict_to_struct_table<alpaqa::CBFGSParams<alpaqa::EigenConfigf>>;
#endif