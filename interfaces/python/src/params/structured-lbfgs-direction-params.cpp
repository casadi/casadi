#include "structured-lbfgs-direction-params.hpp"

template <alpaqa::Config Conf>
const dict_to_struct_table_t<alpaqa::StructuredLBFGSDirectionParams<Conf>>
    dict_to_struct_table<alpaqa::StructuredLBFGSDirectionParams<Conf>>::table{
        // clang-format off
        {"full_augmented_hessian", &alpaqa::StructuredLBFGSDirectionParams<Conf>::full_augmented_hessian},
        {"hessian_vec", &alpaqa::StructuredLBFGSDirectionParams<Conf>::hessian_vec},
        {"hessian_vec_finite_differences", &alpaqa::StructuredLBFGSDirectionParams<Conf>::hessian_vec_finite_differences},
        // clang-format on
    };

// clang-format off
template struct dict_to_struct_table<alpaqa::StructuredLBFGSDirectionParams<alpaqa::EigenConfigf>>;
template struct dict_to_struct_table<alpaqa::StructuredLBFGSDirectionParams<alpaqa::EigenConfigd>>;
template struct dict_to_struct_table<alpaqa::StructuredLBFGSDirectionParams<alpaqa::EigenConfigl>>;
#ifdef ALPAQA_WITH_QUAD_PRECISION
template struct dict_to_struct_table<alpaqa::StructuredLBFGSDirectionParams<alpaqa::EigenConfigf>>;
#endif
// clang-format on