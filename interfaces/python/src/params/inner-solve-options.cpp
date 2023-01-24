#include "params.hpp"

template <alpaqa::Config Conf>
const dict_to_struct_table_t<alpaqa::InnerSolveOptions<Conf>>
    dict_to_struct_table<alpaqa::InnerSolveOptions<Conf>>::table{
        // clang-format off
        {"always_overwrite_results", &alpaqa::InnerSolveOptions<Conf>::always_overwrite_results},
        {"max_time", &alpaqa::InnerSolveOptions<Conf>::max_time},
        {"tolerance", &alpaqa::InnerSolveOptions<Conf>::tolerance},
        // clang-format on
    };

template struct dict_to_struct_table<alpaqa::InnerSolveOptions<alpaqa::EigenConfigf>>;
template struct dict_to_struct_table<alpaqa::InnerSolveOptions<alpaqa::EigenConfigd>>;
template struct dict_to_struct_table<alpaqa::InnerSolveOptions<alpaqa::EigenConfigl>>;
#ifdef ALPAQA_WITH_QUAD_PRECISION
template struct dict_to_struct_table<alpaqa::InnerSolveOptions<alpaqa::EigenConfigq>>;
#endif
