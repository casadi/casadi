#pragma once

#include <alpaqa/inner/directions/panoc/structured-lbfgs.hpp>
#include <kwargs-to-struct.hpp>

template <alpaqa::Config Conf>
struct dict_to_struct_table<alpaqa::StructuredLBFGSDirectionParams<Conf>> {
    static const dict_to_struct_table_t<alpaqa::StructuredLBFGSDirectionParams<Conf>> table;
};

// clang-format off
extern template struct dict_to_struct_table<alpaqa::StructuredLBFGSDirectionParams<alpaqa::EigenConfigf>>;
extern template struct dict_to_struct_table<alpaqa::StructuredLBFGSDirectionParams<alpaqa::EigenConfigd>>;
extern template struct dict_to_struct_table<alpaqa::StructuredLBFGSDirectionParams<alpaqa::EigenConfigl>>;
#ifdef ALPAQA_WITH_QUAD_PRECISION
extern template struct dict_to_struct_table<alpaqa::StructuredLBFGSDirectionParams<alpaqa::EigenConfigq>>;
#endif
// clang-format on
