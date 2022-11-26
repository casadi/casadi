#pragma once

#include <alpaqa/inner/directions/panoc/lbfgs.hpp>
#include <kwargs-to-struct.hpp>

template <alpaqa::Config Conf>
struct dict_to_struct_table<alpaqa::LBFGSDirectionParams<Conf>> {
    static const dict_to_struct_table_t<alpaqa::LBFGSDirectionParams<Conf>> table;
};

extern template struct dict_to_struct_table<alpaqa::LBFGSDirectionParams<alpaqa::EigenConfigf>>;
extern template struct dict_to_struct_table<alpaqa::LBFGSDirectionParams<alpaqa::EigenConfigd>>;
extern template struct dict_to_struct_table<alpaqa::LBFGSDirectionParams<alpaqa::EigenConfigl>>;
#ifdef ALPAQA_WITH_QUAD_PRECISION
extern template struct dict_to_struct_table<alpaqa::LBFGSDirectionParams<alpaqa::EigenConfigf>>;
#endif
