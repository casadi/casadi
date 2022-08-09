#pragma once

#include <alpaqa/accelerators/lbfgs.hpp>
#include "kwargs-to-struct.hpp"

template <alpaqa::Config Conf>
struct kwargs_to_struct_table<alpaqa::LBFGSParams<Conf>> {
    inline static const kwargs_to_struct_table_t<alpaqa::LBFGSParams<Conf>> table {
        {"memory", &alpaqa::LBFGSParams<Conf>::memory},
        {"cbfgs", &alpaqa::LBFGSParams<Conf>::cbfgs},
    };
};