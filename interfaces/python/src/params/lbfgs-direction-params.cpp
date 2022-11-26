#include "lbfgs-direction-params.hpp"

template <alpaqa::Config Conf>
const dict_to_struct_table_t<alpaqa::LBFGSDirectionParams<Conf>>
    dict_to_struct_table<alpaqa::LBFGSDirectionParams<Conf>>::table{
        // clang-format off
        {"rescale_when_γ_changes", &alpaqa::LBFGSDirectionParams<Conf>::rescale_when_γ_changes},
        // clang-format on
    };

template struct dict_to_struct_table<alpaqa::LBFGSDirectionParams<alpaqa::EigenConfigf>>;
template struct dict_to_struct_table<alpaqa::LBFGSDirectionParams<alpaqa::EigenConfigd>>;
template struct dict_to_struct_table<alpaqa::LBFGSDirectionParams<alpaqa::EigenConfigl>>;
#ifdef ALPAQA_WITH_QUAD_PRECISION
template struct dict_to_struct_table<alpaqa::LBFGSDirectionParams<alpaqa::EigenConfigf>>;
#endif