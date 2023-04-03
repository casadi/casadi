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
#include "type-erased-alm-solver.hpp"
#include "type-erased-inner-solver.hpp"

template <alpaqa::Config Conf>
void register_alm(py::module_ &m) {
    USING_ALPAQA_CONFIG(Conf);

    using TEProblem          = alpaqa::TypeErasedProblem<config_t>;
    using TEOCProblem        = alpaqa::TypeErasedControlProblem<config_t>;
    using InnerSolver        = alpaqa::TypeErasedInnerSolver<config_t, TEProblem>;
    using InnerOCPSolver     = alpaqa::TypeErasedInnerSolver<config_t, TEOCProblem>;
    using DefaultInnerSolver = alpaqa::PANOCSolver<alpaqa::StructuredLBFGSDirection<config_t>>;

    using ALMSolver    = alpaqa::ALMSolver<InnerSolver>;
    using ALMOCPSolver = alpaqa::ALMSolver<InnerOCPSolver>;
    using ALMParams    = alpaqa::ALMParams<config_t>;
    using TEALMSolver  = alpaqa::TypeErasedALMSolver<config_t>;
    register_dataclass<ALMParams>(m, "ALMParams",
                                  "C++ documentation: :cpp:class:`alpaqa::ALMParams`");

    py::class_<TEALMSolver> almsolver(m, "ALMSolver",
                                      "Main augmented Lagrangian solver.\n\n"
                                      "C++ documentation: :cpp:class:`alpaqa::ALMSolver`");
    default_copy_methods(almsolver);
    almsolver
        // Default constructor
        .def(py::init([] {
                 return std::make_unique<TEALMSolver>(
                     alpaqa::util::te_in_place<ALMSolver>, ALMParams{},
                     InnerSolver::template make<DefaultInnerSolver>(
                         alpaqa::PANOCParams<config_t>{}));
             }),
             "Build an ALM solver using Structured PANOC as inner solver.")
        // Solver only
        .def(py::init([](const InnerSolver &inner) {
                 return std::make_unique<TEALMSolver>(alpaqa::util::te_in_place<ALMSolver>,
                                                      ALMParams{}, inner);
             }),
             "inner_solver"_a, "Build an ALM solver using the given inner solver.")
        .def(py::init([](const InnerOCPSolver &inner) {
                 return std::make_unique<TEALMSolver>(alpaqa::util::te_in_place<ALMOCPSolver>,
                                                      ALMParams{}, inner);
             }),
             "inner_solver"_a, "Build an ALM solver using the given inner solver.")
        // Params and solver
        .def(py::init([](params_or_dict<ALMParams> params, const InnerSolver &inner) {
                 return std::make_unique<TEALMSolver>(alpaqa::util::te_in_place<ALMSolver>,
                                                      var_kwargs_to_struct(params), inner);
             }),
             "alm_params"_a, "inner_solver"_a, "Build an ALM solver using the given inner solver.")
        .def(py::init([](params_or_dict<ALMParams> params, const InnerOCPSolver &inner) {
                 return std::make_unique<TEALMSolver>(alpaqa::util::te_in_place<ALMOCPSolver>,
                                                      var_kwargs_to_struct(params), inner);
             }),
             "alm_params"_a, "inner_solver"_a, "Build an ALM solver using the given inner solver.")
        // Other functions and properties
        .def_property_readonly("inner_solver", &TEALMSolver::get_inner_solver)
        .def("__call__", &TEALMSolver::operator(), "problem"_a, "x"_a = std::nullopt,
             "y"_a = std::nullopt, "asynchronous"_a = true,
             "Solve.\n\n"
             ":param problem: Problem to solve.\n"
             ":param x: Initial guess for decision variables :math:`x`\n\n"
             ":param y: Initial guess for Lagrange multipliers :math:`y`\n"
             ":param asynchronous: Release the GIL and run the solver on a separate thread\n"
             ":return: * Solution :math:`x`\n"
             "         * Lagrange multipliers :math:`y` at the solution\n"
             "         * Statistics\n\n")
        .def("stop", &TEALMSolver::stop)
        .def_property_readonly("name", &TEALMSolver::get_name)
        .def("__str__", &TEALMSolver::get_name)
        .def_property_readonly("params", &TEALMSolver::get_params);
}

template void register_alm<alpaqa::EigenConfigd>(py::module_ &);
template void register_alm<alpaqa::EigenConfigf>(py::module_ &);
template void register_alm<alpaqa::EigenConfigl>(py::module_ &);
#ifdef ALPAQA_WITH_QUAD_PRECISION
template void register_alm<alpaqa::EigenConfigq>(py::module_ &);
#endif
