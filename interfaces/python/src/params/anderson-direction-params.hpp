#pragma once

#include <alpaqa/inner/directions/panoc/anderson.hpp>
#include <kwargs-to-struct.hpp>

template <alpaqa::Config Conf>
struct dict_to_struct_table<alpaqa::AndersonDirectionParams<Conf>> {
    static const dict_to_struct_table_t<alpaqa::AndersonDirectionParams<Conf>> table;
};

extern template struct dict_to_struct_table<alpaqa::AndersonDirectionParams<alpaqa::EigenConfigf>>;
extern template struct dict_to_struct_table<alpaqa::AndersonDirectionParams<alpaqa::EigenConfigd>>;
extern template struct dict_to_struct_table<alpaqa::AndersonDirectionParams<alpaqa::EigenConfigl>>;
#ifdef ALPAQA_WITH_QUAD_PRECISION
extern template struct dict_to_struct_table<alpaqa::AndersonDirectionParams<alpaqa::EigenConfigq>>;
#endif
