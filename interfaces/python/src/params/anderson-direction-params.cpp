#include "params.hpp"

template <alpaqa::Config Conf>
const dict_to_struct_table_t<alpaqa::AndersonDirectionParams<Conf>>
    dict_to_struct_table<alpaqa::AndersonDirectionParams<Conf>>::table{
        // clang-format off
        {"rescale_when_γ_changes", &alpaqa::AndersonDirectionParams<Conf>::rescale_when_γ_changes},
        // clang-format on
    };

template struct dict_to_struct_table<alpaqa::AndersonDirectionParams<alpaqa::EigenConfigf>>;
template struct dict_to_struct_table<alpaqa::AndersonDirectionParams<alpaqa::EigenConfigd>>;
template struct dict_to_struct_table<alpaqa::AndersonDirectionParams<alpaqa::EigenConfigl>>;
#ifdef ALPAQA_WITH_QUAD_PRECISION
template struct dict_to_struct_table<alpaqa::AndersonDirectionParams<alpaqa::EigenConfigq>>;
#endif