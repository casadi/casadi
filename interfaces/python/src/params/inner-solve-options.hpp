#pragma once

#include <alpaqa/inner/inner-solve-options.hpp>
#include "kwargs-to-struct.hpp"

template <alpaqa::Config Conf>
struct dict_to_struct_table<alpaqa::InnerSolveOptions<Conf>> {
    static const dict_to_struct_table_t<alpaqa::InnerSolveOptions<Conf>> table;
};

extern template struct dict_to_struct_table<alpaqa::InnerSolveOptions<alpaqa::EigenConfigf>>;
extern template struct dict_to_struct_table<alpaqa::InnerSolveOptions<alpaqa::EigenConfigd>>;
extern template struct dict_to_struct_table<alpaqa::InnerSolveOptions<alpaqa::EigenConfigl>>;
#ifdef ALPAQA_WITH_QUAD_PRECISION
extern template struct dict_to_struct_table<alpaqa::InnerSolveOptions<alpaqa::EigenConfigq>>;
#endif
