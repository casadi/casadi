#pragma once

#include <alpaqa/util/quadmath/quadmath.hpp>

#include <pybind11/stl.h>
namespace py = pybind11;

#include <alpaqa/config/config.hpp>
#include <alpaqa/inner/inner-solve-options.hpp>
#include <alpaqa/util/check-dim.hpp>

#include "stats-to-dict.hpp"

template <class Solver, class Problem>
auto checked_inner_solve() {
    USING_ALPAQA_CONFIG_TEMPLATE(Solver::config_t);
    return [](Solver &solver, const Problem &problem,
              const alpaqa::InnerSolveOptions<config_t> &opts, std::optional<vec> u,
              std::optional<vec> y, std::optional<vec> Σ, bool async) {
        bool ret_y = y.has_value();
        if (!y && problem.get_m() > 0)
            throw std::invalid_argument("Missing argument y");
        alpaqa::util::check_dim_msg<config_t>(y, problem.get_m(),
                                              "Length of y does not match problem size problem.m");
        if (!Σ && problem.get_m() > 0)
            throw std::invalid_argument("Missing argument Σ");
        alpaqa::util::check_dim_msg<config_t>(Σ, problem.get_m(),
                                              "Length of Σ does not match problem size problem.m");
        vec err_z          = vec::Zero(problem.get_m());
        auto invoke_solver = [&] { return solver(problem, opts, *u, *y, *Σ, err_z); };
        auto &&stats       = async_solve(async, solver, invoke_solver, problem);
        return ret_y ? py::make_tuple(std::move(*u), std::move(*y), std::move(err_z),
                                      alpaqa::conv::stats_to_dict(stats))
                     : py::make_tuple(std::move(*u), alpaqa::conv::stats_to_dict(stats));
    };
}
