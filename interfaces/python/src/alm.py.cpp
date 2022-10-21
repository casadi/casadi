#include <alpaqa/util/quadmath/quadmath.hpp>

#include <pybind11/chrono.h>
#include <pybind11/eigen.h>
#include <pybind11/gil.h>
#include <pybind11/iostream.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
using namespace py::literals;
constexpr auto ret_ref_internal = py::return_value_policy::reference_internal;

#include <chrono>
#include <exception>
#include <future>
#include <stdexcept>
using namespace std::chrono_literals;

#include <alpaqa/inner/panoc.hpp>
#include <alpaqa/inner/src/panoc.tpp>
#include <alpaqa/outer/alm.hpp>
#include <alpaqa/outer/src/alm.tpp>

#include "check-dim.hpp"
#include "kwargs-to-struct.hpp"
#include "type-erased-inner-solver.hpp"
#include "type-erased-panoc-direction.hpp"

template <alpaqa::Config Conf>
struct kwargs_to_struct_table<alpaqa::ALMParams<Conf>> {
    inline static const kwargs_to_struct_table_t<alpaqa::ALMParams<Conf>> table{
        {"ε", &alpaqa::ALMParams<Conf>::ε},
        {"δ", &alpaqa::ALMParams<Conf>::δ},
        {"Δ", &alpaqa::ALMParams<Conf>::Δ},
        {"Δ_lower", &alpaqa::ALMParams<Conf>::Δ_lower},
        {"Δ_min", &alpaqa::ALMParams<Conf>::Δ_min},
        {"Σ_0", &alpaqa::ALMParams<Conf>::Σ_0},
        {"σ_0", &alpaqa::ALMParams<Conf>::σ_0},
        {"Σ_0_lower", &alpaqa::ALMParams<Conf>::Σ_0_lower},
        {"ε_0", &alpaqa::ALMParams<Conf>::ε_0},
        {"ε_0_increase", &alpaqa::ALMParams<Conf>::ε_0_increase},
        {"ρ", &alpaqa::ALMParams<Conf>::ρ},
        {"ρ_increase", &alpaqa::ALMParams<Conf>::ρ_increase},
        {"ρ_max", &alpaqa::ALMParams<Conf>::ρ_max},
        {"θ", &alpaqa::ALMParams<Conf>::θ},
        {"M", &alpaqa::ALMParams<Conf>::M},
        {"Σ_max", &alpaqa::ALMParams<Conf>::Σ_max},
        {"Σ_min", &alpaqa::ALMParams<Conf>::Σ_min},
        {"max_iter", &alpaqa::ALMParams<Conf>::max_iter},
        {"max_time", &alpaqa::ALMParams<Conf>::max_time},
        {"max_num_initial_retries", &alpaqa::ALMParams<Conf>::max_num_initial_retries},
        {"max_num_retries", &alpaqa::ALMParams<Conf>::max_num_retries},
        {"max_total_num_retries", &alpaqa::ALMParams<Conf>::max_total_num_retries},
        {"print_interval", &alpaqa::ALMParams<Conf>::print_interval},
        {"single_penalty_factor", &alpaqa::ALMParams<Conf>::single_penalty_factor},
    };
};

template <alpaqa::Config Conf>
void register_alm(py::module_ &m) {
    USING_ALPAQA_CONFIG(Conf);

    using TypeErasedPANOCDirection   = alpaqa::TypeErasedPANOCDirection<config_t>;
    using PANOCSolver                = alpaqa::PANOCSolver<TypeErasedPANOCDirection>;
    using StructuredPANOCLBFGSSolver = alpaqa::StructuredPANOCLBFGSSolver<config_t>;
    using InnerSolver                = alpaqa::TypeErasedInnerSolver<config_t>;
    py::class_<InnerSolver>(m, "InnerSolver")
        .def(py::init<PANOCSolver>())
        .def("__call__",
             [](InnerSolver &self, const alpaqa::ProblemBase<config_t> &p, crvec Σ, real_t ε,
                bool a, rvec x, rvec y, rvec e) { return self(p, Σ, ε, a, x, y, e).to_dict(); })
        .def_property_readonly("name", &InnerSolver::template get_name<>);

    using ALMSolver = alpaqa::ALMSolver<InnerSolver>;
    using ALMParams = typename ALMSolver::Params;
    py::class_<ALMParams>(m, "ALMParams", "C++ documentation: :cpp:class:`alpaqa::ALMParams`")
        .def(py::init())
        .def(py::init(&kwargs_to_struct<ALMParams>))
        .def("to_dict", &struct_to_dict<ALMParams>)
        .def_readwrite("ε", &ALMParams::ε)
        .def_readwrite("δ", &ALMParams::δ)
        .def_readwrite("Δ", &ALMParams::Δ)
        .def_readwrite("Δ_lower", &ALMParams::Δ_lower)
        .def_readwrite("Δ_min", &ALMParams::Δ_min)
        .def_readwrite("Σ_0", &ALMParams::Σ_0)
        .def_readwrite("σ_0", &ALMParams::σ_0)
        .def_readwrite("Σ_0_lower", &ALMParams::Σ_0_lower)
        .def_readwrite("ε_0", &ALMParams::ε_0)
        .def_readwrite("ε_0_increase", &ALMParams::ε_0_increase)
        .def_readwrite("ρ", &ALMParams::ρ)
        .def_readwrite("ρ_increase", &ALMParams::ρ_increase)
        .def_readwrite("ρ_max", &ALMParams::ρ_max)
        .def_readwrite("θ", &ALMParams::θ)
        .def_readwrite("M", &ALMParams::M)
        .def_readwrite("Σ_max", &ALMParams::Σ_max)
        .def_readwrite("Σ_min", &ALMParams::Σ_min)
        .def_readwrite("max_iter", &ALMParams::max_iter)
        .def_readwrite("max_time", &ALMParams::max_time)
        .def_readwrite("max_num_initial_retries", &ALMParams::max_num_initial_retries)
        .def_readwrite("max_num_retries", &ALMParams::max_num_retries)
        .def_readwrite("max_total_num_retries", &ALMParams::max_total_num_retries)
        .def_readwrite("print_interval", &ALMParams::print_interval)
        .def_readwrite("single_penalty_factor", &ALMParams::single_penalty_factor);

    auto safe_alm_call = [](ALMSolver &solver, const alpaqa::ProblemBase<config_t> &p,
                            std::optional<vec> x, std::optional<vec> y,
                            bool async) -> std::tuple<vec, vec, py::dict> {
        if (!x)
            x = vec::Zero(p.n);
        else
            check_dim_msg(*x, p.n, "Length of x does not match problem size problem.n");
        if (!y)
            y = vec::Zero(p.m);
        else
            check_dim_msg(*y, p.m, "Length of y does not match problem size problem.m");
        check_dim_msg(p.get_C().lowerbound, p.n,
                      "Length of problem.C.lowerbound does not match problem size problem.n");
        check_dim_msg(p.get_C().upperbound, p.n,
                      "Length of problem.C.upperbound does not match problem size problem.n");
        check_dim_msg(p.get_D().lowerbound, p.m,
                      "Length of problem.D.lowerbound does not match problem size problem.m");
        check_dim_msg(p.get_D().upperbound, p.m,
                      "Length of problem.D.upperbound does not match problem size problem.m");

        auto invoke_solver = [&] { return solver(p, *y, *x); };
        if (!async) {
            auto stats = invoke_solver();
            return std::make_tuple(std::move(*x), std::move(*y),
                                   alpaqa::conv::stats_to_dict<InnerSolver>(stats));
        } else {
            auto stats = std::async(std::launch::async, invoke_solver);
            {
                py::gil_scoped_release gil{};
                while (stats.wait_for(50ms) != std::future_status::ready) {
                    py::gil_scoped_acquire gil{};
                    // Check if Python received a signal (e.g. Ctrl+C)
                    if (PyErr_CheckSignals() != 0) {
                        // Nicely ask the solver to stop
                        solver.stop();
                        // It should return a result soon
                        if (py::gil_scoped_release gil{};
                            stats.wait_for(15s) != std::future_status::ready) {
                            // If it doesn't, we terminate the entire program,
                            // because the solver uses variables local to this
                            // function, so we cannot safely return without
                            // waiting for the solver to finish.
                            std::terminate();
                        }
                        if (PyErr_Occurred())
                            throw py::error_already_set();
                        break;
                    }
                }
            }
            return std::make_tuple(std::move(*x), std::move(*y),
                                   alpaqa::conv::stats_to_dict<InnerSolver>(stats.get()));
        }
    };

    py::class_<ALMSolver>(m, "ALMSolver",
                          "Main augmented Lagrangian solver.\n\n"
                          "C++ documentation: :cpp:class:`alpaqa::ALMSolver`")
        // Default constructor
        .def(py::init([] {
                 return std::make_unique<ALMSolver>(
                     ALMParams{}, InnerSolver{StructuredPANOCLBFGSSolver{{}, {}}});
             }),
             "Build an ALM solver using Structured PANOC as inner solver.")
        // Solver only
        .def(py::init([](const PANOCSolver &inner) {
                 return std::make_unique<ALMSolver>(ALMParams{}, InnerSolver{inner});
             }),
             "inner_solver"_a, "Build an ALM solver using PANOC as inner solver.")
        .def(py::init([](const StructuredPANOCLBFGSSolver &inner) {
                 return std::make_unique<ALMSolver>(ALMParams{}, InnerSolver{inner});
             }),
             "inner_solver"_a, "Build an ALM solver using Structured PANOC as inner solver.")
        // Params and solver
        .def(py::init([](params_or_dict<ALMParams> params, const PANOCSolver &inner) {
                 return std::make_unique<ALMSolver>(var_kwargs_to_struct(params),
                                                    InnerSolver{inner});
             }),
             "alm_params"_a, "inner_solver"_a, "Build an ALM solver using PANOC as inner solver.")
        .def(
            py::init([](params_or_dict<ALMParams> params, const StructuredPANOCLBFGSSolver &inner) {
                return std::make_unique<ALMSolver>(var_kwargs_to_struct(params),
                                                   InnerSolver{inner});
            }),
            "alm_params"_a, "inner_solver"_a,
            "Build an ALM solver using Structured PANOC as inner solver.")
        // Other functions and properties
        .def_property_readonly(
            "inner_solver",
            py::cpp_function([](ALMSolver &self) -> InnerSolver & { return self.inner_solver; },
                             ret_ref_internal))
        .def("__call__", safe_alm_call, "problem"_a, "x"_a = std::nullopt, "y"_a = std::nullopt,
             "async_"_a = false,
             py::call_guard<py::scoped_ostream_redirect, py::scoped_estream_redirect>(),
             "Solve.\n\n"
             ":param problem: Problem to solve.\n"
             ":param x: Initial guess for decision variables :math:`x`\n\n"
             ":param y: Initial guess for Lagrange multipliers :math:`y`\n"
             ":param async_: Release the GIL and run the solver on a separate thread\n"
             ":return: * Solution :math:`x`\n"
             "         * Lagrange multipliers :math:`y` at the solution\n"
             "         * Statistics\n\n")
        .def("__str__", &ALMSolver::get_name)
        .def_property_readonly("params", &ALMSolver::get_params);
}

template void register_alm<alpaqa::EigenConfigd>(py::module_ &);
template void register_alm<alpaqa::EigenConfigf>(py::module_ &);
template void register_alm<alpaqa::EigenConfigl>(py::module_ &);
#ifdef ALPAQA_WITH_QUAD_PRECISION
template void register_alm<alpaqa::EigenConfigq>(py::module_ &);
#endif
