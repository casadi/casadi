#pragma once

#include <alpaqa/inner/directions/panoc/structured-newton.hpp>
#include <kwargs-to-struct.hpp>

template <alpaqa::Config Conf>
struct dict_to_struct_table<alpaqa::StructuredNewtonDirectionParams<Conf>> {
    static const dict_to_struct_table_t<alpaqa::StructuredNewtonDirectionParams<Conf>> table;
};

// clang-format off
extern template struct dict_to_struct_table<alpaqa::StructuredNewtonDirectionParams<alpaqa::EigenConfigf>>;
extern template struct dict_to_struct_table<alpaqa::StructuredNewtonDirectionParams<alpaqa::EigenConfigd>>;
extern template struct dict_to_struct_table<alpaqa::StructuredNewtonDirectionParams<alpaqa::EigenConfigl>>;
#ifdef ALPAQA_WITH_QUAD_PRECISION
extern template struct dict_to_struct_table<alpaqa::StructuredNewtonDirectionParams<alpaqa::EigenConfigq>>;
#endif
// clang-format on
