#include <alpaqa/inner/directions/pantr/newton-tr.hpp>
#include <alpaqa/inner/pantr.hpp>

#include "alm-driver.hpp"
#include "cancel.hpp"
#include "pantr-driver.hpp"
#include "solver-driver.hpp"
#include "util.hpp"

namespace {

template <class T>
struct tag_t {};

template <template <class Direction> class Solver>
solver_func_t make_pantr_like_solver(std::string_view direction,
                                     [[maybe_unused]] Options &opts) {
    USING_ALPAQA_CONFIG(alpaqa::DefaultConfig);
    auto builder = []<class Direction>(tag_t<Direction>) {
        return [](std::string_view, Options &opts) -> solver_func_t {
            auto inner_solver = make_inner_solver<Solver<Direction>>(opts);
            auto solver       = make_alm_solver(std::move(inner_solver), opts);
            unsigned N_exp    = 0;
            set_params(N_exp, "num_exp", opts);
            return [solver{std::move(solver)},
                    N_exp](LoadedProblem &problem,
                           std::ostream &os) mutable -> SolverResults {
                static std::atomic<decltype(solver) *> solver_to_stop;
                attach_cancellation<solver_to_stop>(solver);
                return run_alm_solver(problem, solver, os, N_exp);
            };
        };
    };
    std::map<std::string_view, solver_builder_func_t> builders{
        {"newtontr", //
         builder(tag_t<alpaqa::NewtonTRDirection<config_t>>())},
    };
    if (direction.empty())
        direction = "newtontr";
    auto builder_it = builders.find(direction);
    if (builder_it != builders.end())
        return builder_it->second(direction, opts);
    else
        throw std::invalid_argument(
            "Unknown direction '" + std::string(direction) + "'\n" +
            "  Available directions: " +
            format_string_list(builders,
                               [](const auto &x) { return x.first; }));
}

} // namespace

solver_func_t make_pantr_driver(std::string_view direction, Options &opts) {
    return make_pantr_like_solver<alpaqa::PANTRSolver>(direction, opts);
}
