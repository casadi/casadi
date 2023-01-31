#include "structured-newton-direction-params.hpp"

template <alpaqa::Config Conf>
const dict_to_struct_table_t<alpaqa::StructuredNewtonDirectionParams<Conf>>
    dict_to_struct_table<alpaqa::StructuredNewtonDirectionParams<Conf>>::table{
        // clang-format off
        {"hessian_vec", &alpaqa::StructuredNewtonDirectionParams<Conf>::hessian_vec},
        {"min_eig", &alpaqa::StructuredNewtonDirectionParams<Conf>::min_eig},
        {"print_eig", &alpaqa::StructuredNewtonDirectionParams<Conf>::print_eig},
        // clang-format on
    };

// clang-format off
template struct dict_to_struct_table<alpaqa::StructuredNewtonDirectionParams<alpaqa::EigenConfigf>>;
template struct dict_to_struct_table<alpaqa::StructuredNewtonDirectionParams<alpaqa::EigenConfigd>>;
template struct dict_to_struct_table<alpaqa::StructuredNewtonDirectionParams<alpaqa::EigenConfigl>>;
#ifdef ALPAQA_WITH_QUAD_PRECISION
template struct dict_to_struct_table<alpaqa::StructuredNewtonDirectionParams<alpaqa::EigenConfigq>>;
#endif
// clang-format on