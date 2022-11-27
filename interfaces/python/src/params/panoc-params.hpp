#pragma once

#include <alpaqa/inner/panoc.hpp>
#include <kwargs-to-struct.hpp>

template <alpaqa::Config Conf>
struct dict_to_struct_table<alpaqa::PANOCParams<Conf>> {
    static const dict_to_struct_table_t<alpaqa::PANOCParams<Conf>> table;
};

template <alpaqa::Config Conf>
struct dict_to_struct_table<alpaqa::LipschitzEstimateParams<Conf>> {
    static const dict_to_struct_table_t<alpaqa::LipschitzEstimateParams<Conf>> table;
};

extern template struct dict_to_struct_table<alpaqa::PANOCParams<alpaqa::EigenConfigf>>;
extern template struct dict_to_struct_table<alpaqa::PANOCParams<alpaqa::EigenConfigd>>;
extern template struct dict_to_struct_table<alpaqa::PANOCParams<alpaqa::EigenConfigl>>;
#ifdef ALPAQA_WITH_QUAD_PRECISION
extern template struct dict_to_struct_table<alpaqa::PANOCParams<alpaqa::EigenConfigq>>;
#endif

extern template struct dict_to_struct_table<alpaqa::LipschitzEstimateParams<alpaqa::EigenConfigf>>;
extern template struct dict_to_struct_table<alpaqa::LipschitzEstimateParams<alpaqa::EigenConfigd>>;
extern template struct dict_to_struct_table<alpaqa::LipschitzEstimateParams<alpaqa::EigenConfigl>>;
#ifdef ALPAQA_WITH_QUAD_PRECISION
extern template struct dict_to_struct_table<alpaqa::LipschitzEstimateParams<alpaqa::EigenConfigq>>;
#endif
