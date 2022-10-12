#pragma once

#include <alpaqa/accelerators/lbfgs.hpp>
#include "kwargs-to-struct.hpp"

template <alpaqa::Config Conf>
struct kwargs_to_struct_table<alpaqa::LBFGSParams<Conf>> {
    inline static const kwargs_to_struct_table_t<alpaqa::LBFGSParams<Conf>> table{
        {"memory", &alpaqa::LBFGSParams<Conf>::memory},
        {"min_div_fac", &alpaqa::LBFGSParams<Conf>::min_div_fac},
        {"min_abs_s", &alpaqa::LBFGSParams<Conf>::min_abs_s},
        {"cbfgs", &alpaqa::LBFGSParams<Conf>::cbfgs},
        {"force_pos_def", &alpaqa::LBFGSParams<Conf>::force_pos_def},
        {"stepsize", &alpaqa::LBFGSParams<Conf>::stepsize},
    };
};