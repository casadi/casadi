#include <params/anderson-params.hpp>

template <alpaqa::Config Conf>
const dict_to_struct_table_t<alpaqa::AndersonAccelParams<Conf>>
    dict_to_struct_table<alpaqa::AndersonAccelParams<Conf>>::table{
        {"memory", &alpaqa::AndersonAccelParams<Conf>::memory},
        {"min_div_fac", &alpaqa::AndersonAccelParams<Conf>::min_div_fac},
    };

template struct dict_to_struct_table<alpaqa::AndersonAccelParams<alpaqa::EigenConfigf>>;
template struct dict_to_struct_table<alpaqa::AndersonAccelParams<alpaqa::EigenConfigd>>;
template struct dict_to_struct_table<alpaqa::AndersonAccelParams<alpaqa::EigenConfigl>>;
#ifdef ALPAQA_WITH_QUAD_PRECISION
template struct dict_to_struct_table<alpaqa::AndersonAccelParams<alpaqa::EigenConfigq>>;
#endif
