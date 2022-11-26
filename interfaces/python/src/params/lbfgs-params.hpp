#pragma once

#include <alpaqa/accelerators/lbfgs.hpp>
#include "kwargs-to-struct.hpp"

template <alpaqa::Config Conf>
struct dict_to_struct_table<alpaqa::LBFGSParams<Conf>> {
    static const dict_to_struct_table_t<alpaqa::LBFGSParams<Conf>> table;
};

template <alpaqa::Config Conf>
struct dict_to_struct_table<alpaqa::CBFGSParams<Conf>> {
    static const dict_to_struct_table_t<alpaqa::CBFGSParams<Conf>> table;
};

extern template struct dict_to_struct_table<alpaqa::LBFGSParams<alpaqa::EigenConfigf>>;
extern template struct dict_to_struct_table<alpaqa::LBFGSParams<alpaqa::EigenConfigd>>;
extern template struct dict_to_struct_table<alpaqa::LBFGSParams<alpaqa::EigenConfigl>>;
#ifdef ALPAQA_WITH_QUAD_PRECISION
extern template struct dict_to_struct_table<alpaqa::LBFGSParams<alpaqa::EigenConfigf>>;
#endif

extern template struct dict_to_struct_table<alpaqa::CBFGSParams<alpaqa::EigenConfigf>>;
extern template struct dict_to_struct_table<alpaqa::CBFGSParams<alpaqa::EigenConfigd>>;
extern template struct dict_to_struct_table<alpaqa::CBFGSParams<alpaqa::EigenConfigl>>;
#ifdef ALPAQA_WITH_QUAD_PRECISION
extern template struct dict_to_struct_table<alpaqa::CBFGSParams<alpaqa::EigenConfigf>>;
#endif