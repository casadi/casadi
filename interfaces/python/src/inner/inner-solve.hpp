#pragma once

#include <alpaqa/util/quadmath/quadmath.hpp>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
namespace py = pybind11;
using namespace py::literals;

#include <alpaqa/config/config.hpp>
#include <alpaqa/inner/inner-solve-options.hpp>
#include <alpaqa/util/check-dim.hpp>

#include <dict/stats-to-dict.hpp>
#include <inner/type-erased-inner-solver.hpp>

/// Python interface to the inner solvers, checks the argument sizes and
/// presence, and returns a Python tuple.
template <class Solver, class Problem>
auto checked_inner_solve() {
    USING_ALPAQA_CONFIG_TEMPLATE(Solver::config_t);
    return [](Solver &solver, const Problem &problem,
              const alpaqa::InnerSolveOptions<config_t> &opts, std::optional<vec> x,
              std::optional<vec> y, std::optional<vec> Σ, bool async) {
        alpaqa::util::check_dim_msg<config_t>(x, problem.get_n(),
                                              "Length of x does not match problem size problem.n");
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
        auto invoke_solver = [&] { return solver(problem, opts, *x, *y, *Σ, err_z); };
        auto &&stats       = async_solve(async, solver, invoke_solver, problem);
        return ret_y ? py::make_tuple(std::move(*x), std::move(*y), std::move(err_z),
                                      alpaqa::conv::stats_to_dict(stats))
                     : py::make_tuple(std::move(*x), alpaqa::conv::stats_to_dict(stats));
    };
}

inline const char *checked_inner_solve_doc() {
    return "Solve.\n\n"
           ":param problem: Problem to solve\n"
           ":param opts: Options\n"
           ":param u: Initial guess\n"
           ":param y: Lagrange multipliers\n"
           ":param Σ: Penalty factors\n"
           ":param asynchronous: Release the GIL and run the solver on a separate thread\n"
           ":return: * Solution :math:`u`\n"
           "         * Updated Lagrange multipliers (only if parameter ``y`` was not ``None``)\n"
           "         * Constraint violation (only if parameter ``y`` was not ``None``)\n"
           "         * Statistics\n\n";
}

template <class Solver, class Problem, class InnerSolverType>
void register_inner_solver_methods(py::class_<Solver> &cls) {
    cls.def("__call__", checked_inner_solve<Solver, Problem>(), "problem"_a, "opts"_a = py::dict(),
            "x"_a = py::none(), "y"_a = py::none(), "Σ"_a = py::none(), "asynchronous"_a = true,
            checked_inner_solve_doc())
        .def_property_readonly("name", &Solver::get_name)
        .def("stop", &Solver::stop)
        .def("__str__", &Solver::get_name);
    if constexpr (requires { &Solver::set_progress_callback; })
        cls.def(
            "set_progress_callback", &Solver::set_progress_callback, "callback"_a,
            "Specify a callable that is invoked with some intermediate results on each iteration "
            "of the algorithm.");
    inner_solver_class<InnerSolverType>.template implicitly_convertible_to<Solver>();
}
