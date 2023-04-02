#include <alpaqa/util/quadmath/quadmath.hpp>

#include <pybind11/chrono.h>
#include <pybind11/eigen.h>
#include <pybind11/gil.h>
#include <pybind11/iostream.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
using namespace py::literals;

#include <chrono>
#include <exception>
#include <future>
#include <stdexcept>
using namespace std::chrono_literals;

#include <alpaqa/implementation/inner/panoc.tpp>
#include <alpaqa/implementation/outer/alm.tpp>
#include <alpaqa/inner/directions/panoc/structured-lbfgs.hpp>
#include <alpaqa/inner/panoc.hpp>
#include <alpaqa/outer/alm.hpp>
#include <alpaqa/util/check-dim.hpp>

#include "async.hpp"
#include "copy.hpp"
#include "kwargs-to-struct.hpp"
#include "member.hpp"
#include "params/alm-params.hpp"
#include "type-erased-inner-solver.hpp"
#include "type-erased-panoc-direction.hpp"

template <alpaqa::Config Conf>
void register_alm(py::module_ &m) {
    USING_ALPAQA_CONFIG(Conf);

    using TypeErasedPANOCDirection = alpaqa::TypeErasedPANOCDirection<config_t>;
    using PANOCSolver              = alpaqa::PANOCSolver<TypeErasedPANOCDirection>;
    using TypeErasedProblem        = alpaqa::TypeErasedProblem<config_t>;
    using InnerSolver              = alpaqa::TypeErasedInnerSolver<config_t>;
    using DefaultInnerSolver = alpaqa::PANOCSolver<alpaqa::StructuredLBFGSDirection<config_t>>;
    py::class_<InnerSolver>(m, "InnerSolver")
        .def(py::init<PANOCSolver>())
        .def_property_readonly("name", &InnerSolver::template get_name<>);

    using ALMSolver = alpaqa::ALMSolver<InnerSolver>;
    using ALMParams = typename ALMSolver::Params;
    register_dataclass<ALMParams>(m, "ALMParams",
                                  "C++ documentation: :cpp:class:`alpaqa::ALMParams`");

    auto safe_alm_call = [](ALMSolver &solver, const TypeErasedProblem &p, std::optional<vec> x,
                            std::optional<vec> y, bool async) -> std::tuple<vec, vec, py::dict> {
        alpaqa::util::check_dim_msg<config_t>(x, p.get_n(),
                                              "Length of x does not match problem size problem.n");
        alpaqa::util::check_dim_msg<config_t>(y, p.get_m(),
                                              "Length of y does not match problem size problem.m");
        auto invoke_solver = [&] { return solver(p, *x, *y); };
        auto stats         = async_solve(async, solver, invoke_solver, p);
        return std::make_tuple(std::move(*x), std::move(*y),
                               alpaqa::conv::stats_to_dict<InnerSolver>(std::move(stats)));
    };

    py::class_<ALMSolver> almsolver(m, "ALMSolver",
                                    "Main augmented Lagrangian solver.\n\n"
                                    "C++ documentation: :cpp:class:`alpaqa::ALMSolver`");
    almsolver
        // Default constructor
        .def(py::init([] {
                 return std::make_unique<ALMSolver>(ALMParams{},
                                                    InnerSolver::template make<DefaultInnerSolver>(
                                                        alpaqa::PANOCParams<config_t>{}));
             }),
             "Build an ALM solver using Structured PANOC as inner solver.")
        // Solver only
        .def(py::init([](const PANOCSolver &inner) {
                 return std::make_unique<ALMSolver>(ALMParams{}, InnerSolver{inner});
             }),
             "inner_solver"_a, "Build an ALM solver using PANOC as inner solver.")
        // Params and solver
        .def(py::init([](params_or_dict<ALMParams> params, const PANOCSolver &inner) {
                 return std::make_unique<ALMSolver>(var_kwargs_to_struct(params),
                                                    InnerSolver{inner});
             }),
             "alm_params"_a, "inner_solver"_a, "Build an ALM solver using PANOC as inner solver.")
        // Other functions and properties
        .def_property_readonly("inner_solver", member_ref<&ALMSolver::inner_solver>())
        .def("__call__", safe_alm_call, "problem"_a, "x"_a = std::nullopt, "y"_a = std::nullopt,
             "asynchronous"_a = true,
             "Solve.\n\n"
             ":param problem: Problem to solve.\n"
             ":param x: Initial guess for decision variables :math:`x`\n\n"
             ":param y: Initial guess for Lagrange multipliers :math:`y`\n"
             ":param asynchronous: Release the GIL and run the solver on a separate thread\n"
             ":return: * Solution :math:`x`\n"
             "         * Lagrange multipliers :math:`y` at the solution\n"
             "         * Statistics\n\n")
        .def("__str__", &ALMSolver::get_name)
        .def_property_readonly("params", &ALMSolver::get_params);
    default_copy_methods(almsolver);
}

template void register_alm<alpaqa::EigenConfigd>(py::module_ &);
template void register_alm<alpaqa::EigenConfigf>(py::module_ &);
template void register_alm<alpaqa::EigenConfigl>(py::module_ &);
#ifdef ALPAQA_WITH_QUAD_PRECISION
template void register_alm<alpaqa::EigenConfigq>(py::module_ &);
#endif
