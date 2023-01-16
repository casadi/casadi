#pragma once

#include <alpaqa/inner/experimental-panoc-ocp.hpp>
#include <kwargs-to-struct.hpp>

template <alpaqa::Config Conf>
struct dict_to_struct_table<alpaqa::experimental::PANOCOCPParams<Conf>> {
    static const dict_to_struct_table_t<alpaqa::experimental::PANOCOCPParams<Conf>> table;
};

// clang-format off
extern template struct dict_to_struct_table<alpaqa::experimental::PANOCOCPParams<alpaqa::EigenConfigf>>;
extern template struct dict_to_struct_table<alpaqa::experimental::PANOCOCPParams<alpaqa::EigenConfigd>>;
extern template struct dict_to_struct_table<alpaqa::experimental::PANOCOCPParams<alpaqa::EigenConfigl>>;
#ifdef ALPAQA_WITH_QUAD_PRECISION
extern template struct dict_to_struct_table<alpaqa::experimental::PANOCOCPParams<alpaqa::EigenConfigq>>;
#endif
// clang-format on
