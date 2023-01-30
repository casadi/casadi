#pragma once

#include <alpaqa/accelerators/anderson.hpp>
#include "kwargs-to-struct.hpp"

template <alpaqa::Config Conf>
struct dict_to_struct_table<alpaqa::AndersonAccelParams<Conf>> {
    static const dict_to_struct_table_t<alpaqa::AndersonAccelParams<Conf>> table;
};

extern template struct dict_to_struct_table<alpaqa::AndersonAccelParams<alpaqa::EigenConfigf>>;
extern template struct dict_to_struct_table<alpaqa::AndersonAccelParams<alpaqa::EigenConfigd>>;
extern template struct dict_to_struct_table<alpaqa::AndersonAccelParams<alpaqa::EigenConfigl>>;
#ifdef ALPAQA_WITH_QUAD_PRECISION
extern template struct dict_to_struct_table<alpaqa::AndersonAccelParams<alpaqa::EigenConfigq>>;
#endif
