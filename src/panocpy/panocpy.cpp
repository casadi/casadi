/**
 * @file
 * This file defines all Python bindings.
 */

#include <panoc-alm/decl/alm.hpp>
#include <panoc-alm/inner/decl/panoc-stop-crit.hpp>
#include <panoc-alm/inner/decl/panoc.hpp>
#include <panoc-alm/inner/directions/lbfgs.hpp>
#include <panoc-alm/inner/guarded-aa-pga.hpp>
#include <panoc-alm/inner/panoc.hpp>
#include <panoc-alm/inner/pga.hpp>
#include <panoc-alm/inner/structured-panoc-lbfgs.hpp>
#include <panoc-alm/standalone/panoc.hpp>
#include <panoc-alm/util/solverstatus.hpp>

#if PANOCPY_HAVE_CASADI
#include <panoc-alm/interop/casadi/CasADiLoader.hpp>
#endif

#include <pybind11/attr.h>
#include <pybind11/cast.h>
#include <pybind11/chrono.h>
#include <pybind11/detail/common.h>
#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include <pybind11/iostream.h>
#include <pybind11/pybind11.h>
#include <pybind11/pytypes.h>
#include <pybind11/stl.h>

#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <tuple>
#include <type_traits>
#include <utility>

#include "kwargs-to-struct.hpp"
#include "polymorphic-inner-solver.hpp"
#include "polymorphic-panoc-direction.hpp"
#include "problem.hpp"

namespace py = pybind11;

template <class DirectionProviderT>
auto PolymorphicPANOCConstructor() {
    return [](const pa::PANOCParams &pp, const DirectionProviderT &dir) {
        using Base = pa::PolymorphicPANOCDirectionBase;
        static_assert(std::is_base_of_v<Base, DirectionProviderT>);
        auto full_python_copy = std::make_shared<py::object>(py::cast(dir));
        auto base_copy        = full_python_copy->template cast<Base *>();
        return std::make_shared<pa::PolymorphicPANOCSolver>(
            pa::PANOCSolver<Base>{
                pp,
                std::shared_ptr<Base>(full_python_copy, base_copy),
            });
    };
}

template <class DirectionProviderT, class... DirectionArgumentsT>
auto PolymorphicPANOCConversion() {
    return [](const pa::PANOCParams &pp, const DirectionArgumentsT &...args) {
        using Base = pa::PolymorphicPANOCDirectionBase;
        static_assert(std::is_base_of_v<Base, DirectionProviderT>);
        static_assert(std::is_constructible_v<DirectionProviderT,
                                              DirectionArgumentsT...>);
        DirectionProviderT dir{args...};
        return PolymorphicPANOCConstructor<DirectionProviderT>()(pp, dir);
    };
}

template <class InnerSolverT>
auto PolymorphicALMConstructor() {
    return [](const pa::ALMParams &pp, const InnerSolverT &inner) {
        using Base = pa::PolymorphicInnerSolverBase;
        static_assert(std::is_base_of_v<Base, InnerSolverT>);
        auto full_python_copy = std::make_shared<py::object>(py::cast(inner));
        auto base_copy        = full_python_copy->template cast<Base *>();
        return pa::PolymorphicALMSolver{
            pp,
            std::shared_ptr<Base>(full_python_copy, base_copy),
        };
    };
}

template <class InnerSolverT>
auto PolymorphicALMConstructorDefaultParams() {
    return [](const InnerSolverT &inner) {
        return PolymorphicALMConstructor<InnerSolverT>()(pa::ALMParams(),
                                                         inner);
    };
}

template <class InnerSolverT, class... InnerSolverArgumentsT>
auto PolymorphicALMConversion() {
    return [](const pa::ALMParams &pp, const InnerSolverArgumentsT &...args) {
        using Base = pa::PolymorphicPANOCDirectionBase;
        static_assert(std::is_base_of_v<Base, InnerSolverT>);
        static_assert(
            std::is_constructible_v<InnerSolverT, InnerSolverArgumentsT...>);
        InnerSolverT inner{args...};
        return PolymorphicALMConstructor<InnerSolverT>()(pp, inner);
    };
}

#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)

PYBIND11_MODULE(PANOCPY_MODULE_NAME, m) {
    using py::operator""_a;

    py::options options;
    options.enable_function_signatures();
    options.enable_user_defined_docstrings();

    m.doc() = "PANOC+ALM solvers"; // TODO

#ifdef VERSION_INFO
    m.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);
#else
    m.attr("__version__") = "dev";
#endif

    py::class_<pa::Box>(m, "Box", "C++ documentation: :cpp:class:`pa::Box`")
        .def(py::init([](unsigned n) {
                 return pa::Box{pa::vec::Constant(n, pa::inf),
                                pa::vec::Constant(n, -pa::inf)};
             }),
             "n"_a,
             "Create an :math:`n`-dimensional box at with bounds at "
             ":math:`\\pm\\infty` (no constraints).")
        .def(py::init([](pa::vec ub, pa::vec lb) {
                 if (ub.size() != lb.size())
                     throw std::invalid_argument(
                         "Upper bound and lower bound dimensions do not "
                         "match");
                 return pa::Box{std::move(ub), std::move(lb)};
             }),
             "ub"_a, "lb"_a, "Create a box with the given bounds.")
        .def_readwrite("upperbound", &pa::Box::upperbound)
        .def_readwrite("lowerbound", &pa::Box::lowerbound);

    py::class_<pa::Problem>(m, "Problem",
                            "C++ documentation: :cpp:class:`pa::Problem`")
        // .def(py::init())
        .def(py::init<unsigned, unsigned>(), "n"_a, "m"_a,
             ":param n: Number of unknowns\n"
             ":param m: Number of constraints")
        .def_readwrite("n", &pa::Problem::n,
                       "Number of unknowns, dimension of :math:`x`")
        .def_readwrite(
            "m", &pa::Problem::m,
            "Number of general constraints, dimension of :math:`g(x)`")
        .def_readwrite("C", &pa::Problem::C, "Box constraints on :math:`x`")
        .def_readwrite("D", &pa::Problem::D, "Box constraints on :math:`g(x)`")
        .def_property("f", prob_getter_f(), prob_setter_f(),
                      "Objective funcion, :math:`f(x)`")
        .def_property(
            "grad_f", prob_getter_grad_f(), prob_setter_grad_f(),
            "Gradient of the objective function, :math:`\\nabla f(x)`")
        .def_property("g", prob_getter_g(), prob_setter_g(),
                      "Constraint function, :math:`g(x)`")
        .def_property("grad_g_prod", prob_getter_grad_g_prod(),
                      prob_setter_grad_g_prod(),
                      "Gradient of constraint function times vector, "
                      ":math:`\\nabla g(x)\\, v`")
        .def_property("grad_gi", prob_getter_grad_gi(), prob_setter_grad_gi(),
                      "Gradient vector of the :math:`i`-th component of the "
                      "constriant function, :math:`\\nabla g_i(x)`")
        .def_property(
            "hess_L", prob_getter_hess_L(), prob_setter_hess_L(),
            "Hessian of the Lagrangian function, :math:`\\nabla^2_{xx} L(x,y)`")
        .def_property("hess_L_prod", prob_getter_hess_L_prod(),
                      prob_setter_hess_L_prod(),
                      "Hessian of the Lagrangian function times vector, "
                      ":math:`\\nabla^2_{xx} L(x,y)\\, v`");

    py::class_<pa::ProblemWithParam, pa::Problem>(
        m, "ProblemWithParam",
        "C++ documentation: :cpp:class:`pa::ProblemWithParam`\n\n"
        "See :py:class:`panocpy._panocpy.Problem` for the full documentation.")
        .def(py::init<unsigned, unsigned, unsigned>(), "n"_a, "m"_a, "p"_a)
        .def_readwrite("n", &pa::Problem::n)
        .def_readwrite("m", &pa::Problem::m)
        .def_readwrite("C", &pa::Problem::C)
        .def_readwrite("D", &pa::Problem::D)
        .def_property("f", prob_getter_f(), prob_setter_f())
        .def_property("grad_f", prob_getter_grad_f(), prob_setter_grad_f())
        .def_property("g", prob_getter_g(), prob_setter_g())
        .def_property("grad_g_prod", prob_getter_grad_g_prod(),
                      prob_setter_grad_g_prod())
        .def_property("grad_gi", prob_getter_grad_gi(), prob_setter_grad_gi())
        .def_property("hess_L", prob_getter_hess_L(), prob_setter_hess_L())
        .def_property("hess_L_prod", prob_getter_hess_L_prod(),
                      prob_setter_hess_L_prod())
        .def_property(
            "param", py::overload_cast<>(&pa::ProblemWithParam::get_param),
            [](pa::ProblemWithParam &p, pa::crvec param) {
                if (param.size() != p.get_param().size())
                    throw std::invalid_argument(
                        "Invalid parameter dimension: got " +
                        std::to_string(param.size()) + ", should be " +
                        std::to_string(p.get_param().size()) + ".");
                p.set_param(param);
            },
            "Parameter vector :math:`p` of the problem");

    py::class_<pa::PolymorphicPANOCDirectionBase,
               std::shared_ptr<pa::PolymorphicPANOCDirectionBase>,
               pa::PolymorphicPANOCDirectionTrampoline>(
        m, "PANOCDirection",
        "Class that provides fast directions for the PANOC algorithm (e.g. "
        "L-BFGS)")
        .def(py::init<>())
        .def("initialize", &pa::PolymorphicPANOCDirectionBase::initialize)
        .def("update", &pa::PolymorphicPANOCDirectionBase::update)
        .def("apply", &pa::PolymorphicPANOCDirectionBase::apply_ret)
        .def("changed_γ", &pa::PolymorphicPANOCDirectionBase::changed_γ)
        .def("reset", &pa::PolymorphicPANOCDirectionBase::reset)
        .def("get_name", &pa::PolymorphicPANOCDirectionBase::get_name)
        .def("__str__", &pa::PolymorphicPANOCDirectionBase::get_name);

    using paLBFGSParamCBFGS = decltype(pa::LBFGSParams::cbfgs);
    py::class_<paLBFGSParamCBFGS>(
        m, "LBFGSParamsCBFGS",
        "C++ documentation: :cpp:member:`pa::LBFGSParams::cbfgs`")
        .def(py::init())
        .def(py::init(&kwargs_to_struct<paLBFGSParamCBFGS>))
        .def("to_dict", &struct_to_dict<paLBFGSParamCBFGS>)
        .def_readwrite("α", &paLBFGSParamCBFGS::α)
        .def_readwrite("ϵ", &paLBFGSParamCBFGS::ϵ);

    py::class_<pa::LBFGSParams>(
        m, "LBFGSParams", "C++ documentation: :cpp:class:`pa::LBFGSParams`")
        .def(py::init())
        .def(py::init(&kwargs_to_struct<pa::LBFGSParams>))
        .def("to_dict", &struct_to_dict<pa::LBFGSParams>)
        .def_readwrite("memory", &pa::LBFGSParams::memory)
        .def_readwrite("cbfgs", &pa::LBFGSParams::cbfgs)
        .def_readwrite("rescale_when_γ_changes",
                       &pa::LBFGSParams::rescale_when_γ_changes);

    py::class_<pa::PolymorphicLBFGSDirection,
               std::shared_ptr<pa::PolymorphicLBFGSDirection>,
               pa::PolymorphicPANOCDirectionBase>(
        m, "LBFGSDirection",
        "C++ documentation: :cpp:class:`pa::LBFGSDirection`")
        .def(py::init<pa::LBFGSParams>(), "params"_a)
        .def("initialize", &pa::PolymorphicLBFGSDirection::initialize)
        .def("update", &pa::PolymorphicLBFGSDirection::update)
        .def("apply", &pa::PolymorphicLBFGSDirection::apply_ret)
        .def("changed_γ", &pa::PolymorphicLBFGSDirection::changed_γ)
        .def("reset", &pa::PolymorphicLBFGSDirection::reset)
        .def("get_name", &pa::PolymorphicLBFGSDirection::get_name)
        .def("__str__", &pa::PolymorphicLBFGSDirection::get_name)
        .def_property_readonly("params",
                               &pa::PolymorphicLBFGSDirection::get_params);

    py::enum_<pa::LBFGSStepSize>(
        m, "LBFGSStepsize", "C++ documentation: :cpp:enum:`pa::LBFGSStepSize`")
        .value("BasedOnGradientStepSize",
               pa::LBFGSStepSize::BasedOnGradientStepSize)
        .value("BasedOnCurvature", pa::LBFGSStepSize::BasedOnCurvature)
        .export_values();

    py::class_<pa::LipschitzEstimateParams>(
        m, "LipschitzEstimateParams",
        "C++ documentation: :cpp:class:`pa::LipschitzEstimateParams`")
        .def(py::init())
        .def(py::init(&kwargs_to_struct<pa::LipschitzEstimateParams>))
        .def("to_dict", &struct_to_dict<pa::LipschitzEstimateParams>)
        .def_readwrite("L_0", &pa::LipschitzEstimateParams::L₀)
        .def_readwrite("ε", &pa::LipschitzEstimateParams::ε)
        .def_readwrite("δ", &pa::LipschitzEstimateParams::δ)
        .def_readwrite("Lγ_factor", &pa::LipschitzEstimateParams::Lγ_factor);

    py::class_<pa::PANOCParams>(
        m, "PANOCParams", "C++ documentation: :cpp:class:`pa::PANOCParams`")
        .def(py::init())
        .def(py::init(&kwargs_to_struct<pa::PANOCParams>))
        .def("to_dict", &struct_to_dict<pa::PANOCParams>)
        .def_readwrite("Lipschitz", &pa::PANOCParams::Lipschitz)
        .def_readwrite("max_iter", &pa::PANOCParams::max_iter)
        .def_readwrite("max_time", &pa::PANOCParams::max_time)
        .def_readwrite("τ_min", &pa::PANOCParams::τ_min)
        .def_readwrite("L_min", &pa::PANOCParams::L_min)
        .def_readwrite("L_max", &pa::PANOCParams::L_max)
        .def_readwrite("max_no_progress", &pa::PANOCParams::max_no_progress)
        .def_readwrite("print_interval", &pa::PANOCParams::print_interval)
        .def_readwrite("quadratic_upperbound_tolerance_factor",
                       &pa::PANOCParams::quadratic_upperbound_tolerance_factor)
        .def_readwrite("update_lipschitz_in_linesearch",
                       &pa::PANOCParams::update_lipschitz_in_linesearch)
        .def_readwrite("alternative_linesearch_cond",
                       &pa::PANOCParams::alternative_linesearch_cond)
        .def_readwrite("lbfgs_stepsize", &pa::PANOCParams::lbfgs_stepsize);

    py::enum_<pa::SolverStatus>(
        m, "SolverStatus", py::arithmetic(),
        "C++ documentation: :cpp:enum:`pa::SolverStatus`")
        .value("Unknown", pa::SolverStatus::Unknown, "Initial value")
        .value("Converged", pa::SolverStatus::Converged,
               "Converged and reached given tolerance")
        .value("MaxTime", pa::SolverStatus::MaxTime,
               "Maximum allowed execution time exceeded")
        .value("MaxIter", pa::SolverStatus::MaxIter,
               "Maximum number of iterations exceeded")
        .value("NotFinite", pa::SolverStatus::NotFinite,
               "Intermediate results were infinite or not-a-number")
        .value("NoProgress", pa::SolverStatus::NoProgress,
               "No progress was made in the last iteration")
        .value("Interrupted", pa::SolverStatus::Interrupted,
               "Solver was interrupted by the user")
        .export_values();

    py::class_<pa::PolymorphicInnerSolverBase::Stats>(m, "InnerSolverStats")
        .def(py::init(&pa::PolymorphicInnerSolverBase::Stats::from_dict));

    py::class_<pa::PolymorphicInnerSolverBase,
               std::shared_ptr<pa::PolymorphicInnerSolverBase>,
               pa::PolymorphicInnerSolverTrampoline>(m, "InnerSolver")
        .def(py::init<>())
        .def("__call__", &pa::PolymorphicInnerSolverBase::operator())
        .def("stop", &pa::PolymorphicInnerSolverBase::stop)
        .def("get_name", &pa::PolymorphicInnerSolverBase::get_name)
        .def("get_params", &pa::PolymorphicInnerSolverBase::get_params);

    py::class_<pa::PGAParams>(m, "PGAParams",
                              "C++ documentation: :cpp:class:`pa::PGAParams`")
        .def(py::init())
        .def(py::init(&kwargs_to_struct<pa::PGAParams>))
        .def("to_dict", &struct_to_dict<pa::PGAParams>)
        .def_readwrite("Lipschitz", &pa::PGAParams::Lipschitz)
        .def_readwrite("max_iter", &pa::PGAParams::max_iter)
        .def_readwrite("max_time", &pa::PGAParams::max_time)
        .def_readwrite("L_min", &pa::PGAParams::L_min)
        .def_readwrite("L_max", &pa::PGAParams::L_max)
        .def_readwrite("stop_crit", &pa::PGAParams::stop_crit)
        .def_readwrite("print_interval", &pa::PGAParams::print_interval)
        .def_readwrite("quadratic_upperbound_tolerance_factor",
                       &pa::PGAParams::quadratic_upperbound_tolerance_factor);

    py::class_<pa::PGAProgressInfo>(
        m, "PGAProgressInfo",
        "C++ documentation: :cpp:class:`pa::PGAProgressInfo`")
        .def_readonly("k", &pa::PGAProgressInfo::k)
        .def_readonly("x", &pa::PGAProgressInfo::x)
        .def_readonly("p", &pa::PGAProgressInfo::p)
        .def_readonly("norm_sq_p", &pa::PGAProgressInfo::norm_sq_p)
        .def_readonly("x_hat", &pa::PGAProgressInfo::x_hat)
        .def_readonly("ψ", &pa::PGAProgressInfo::ψ)
        .def_readonly("grad_ψ", &pa::PGAProgressInfo::grad_ψ)
        .def_readonly("ψ_hat", &pa::PGAProgressInfo::ψ_hat)
        .def_readonly("grad_ψ_hat", &pa::PGAProgressInfo::grad_ψ_hat)
        .def_readonly("L", &pa::PGAProgressInfo::L)
        .def_readonly("γ", &pa::PGAProgressInfo::γ)
        .def_readonly("ε", &pa::PGAProgressInfo::ε)
        .def_readonly("Σ", &pa::PGAProgressInfo::Σ)
        .def_readonly("y", &pa::PGAProgressInfo::y)
        .def_property_readonly("fpr", [](const pa::PGAProgressInfo &p) {
            return std::sqrt(p.norm_sq_p) / p.γ;
        });

    py::class_<pa::GAAPGAParams>(
        m, "GAAPGAParams", "C++ documentation: :cpp:class:`pa::GAAPGAParams`")
        .def(py::init())
        .def(py::init(&kwargs_to_struct<pa::GAAPGAParams>))
        .def("to_dict", &struct_to_dict<pa::GAAPGAParams>)
        .def_readwrite("Lipschitz", &pa::GAAPGAParams::Lipschitz)
        .def_readwrite("limitedqr_mem", &pa::GAAPGAParams::limitedqr_mem)
        .def_readwrite("max_iter", &pa::GAAPGAParams::max_iter)
        .def_readwrite("max_time", &pa::GAAPGAParams::max_time)
        .def_readwrite("L_min", &pa::GAAPGAParams::L_min)
        .def_readwrite("L_max", &pa::GAAPGAParams::L_max)
        .def_readwrite("stop_crit", &pa::GAAPGAParams::stop_crit)
        .def_readwrite("print_interval", &pa::GAAPGAParams::print_interval)
        .def_readwrite("quadratic_upperbound_tolerance_factor",
                       &pa::GAAPGAParams::quadratic_upperbound_tolerance_factor)
        .def_readwrite("max_no_progress", &pa::GAAPGAParams::max_no_progress)
        .def_readwrite("full_flush_on_γ_change",
                       &pa::GAAPGAParams::full_flush_on_γ_change);

    py::class_<pa::GAAPGAProgressInfo>(
        m, "GAAPGAProgressInfo",
        "C++ documentation: :cpp:class:`pa::GAAPGAProgressInfo`")
        .def_readonly("k", &pa::GAAPGAProgressInfo::k)
        .def_readonly("x", &pa::GAAPGAProgressInfo::x)
        .def_readonly("p", &pa::GAAPGAProgressInfo::p)
        .def_readonly("norm_sq_p", &pa::GAAPGAProgressInfo::norm_sq_p)
        .def_readonly("x_hat", &pa::GAAPGAProgressInfo::x_hat)
        .def_readonly("ψ", &pa::GAAPGAProgressInfo::ψ)
        .def_readonly("grad_ψ", &pa::GAAPGAProgressInfo::grad_ψ)
        .def_readonly("ψ_hat", &pa::GAAPGAProgressInfo::ψ_hat)
        .def_readonly("grad_ψ_hat", &pa::GAAPGAProgressInfo::grad_ψ_hat)
        .def_readonly("L", &pa::GAAPGAProgressInfo::L)
        .def_readonly("γ", &pa::GAAPGAProgressInfo::γ)
        .def_readonly("ε", &pa::GAAPGAProgressInfo::ε)
        .def_readonly("Σ", &pa::GAAPGAProgressInfo::Σ)
        .def_readonly("y", &pa::GAAPGAProgressInfo::y)
        .def_property_readonly("fpr", [](const pa::GAAPGAProgressInfo &p) {
            return std::sqrt(p.norm_sq_p) / p.γ;
        });

    py::class_<pa::PANOCProgressInfo>(
        m, "PANOCProgressInfo",
        "Data passed to the PANOC progress callback.\n\n"
        "C++ documentation: :cpp:class:`pa::PANOCProgressInfo`")
        .def_readonly("k", &pa::PANOCProgressInfo::k, //
                      "Iteration")
        .def_readonly("x", &pa::PANOCProgressInfo::x, //
                      "Decision variable :math:`x`")
        .def_readonly("p", &pa::PANOCProgressInfo::p, //
                      "Projected gradient step :math:`p`")
        .def_readonly("norm_sq_p", &pa::PANOCProgressInfo::norm_sq_p, //
                      ":math:`\\left\\|p\\right\\|^2`")
        .def_readonly(
            "x_hat", &pa::PANOCProgressInfo::x_hat, //
            "Decision variable after projected gradient step :math:`\\hat x`")
        .def_readonly("φγ", &pa::PANOCProgressInfo::φγ, //
                      "Forward-backward envelope :math:`\\varphi_\\gamma(x)`")
        .def_readonly("ψ", &pa::PANOCProgressInfo::ψ, //
                      "Objective value :math:`\\psi(x)`")
        .def_readonly("grad_ψ", &pa::PANOCProgressInfo::grad_ψ, //
                      "Gradient of objective :math:`\\nabla\\psi(x)`")
        .def_readonly("ψ_hat", &pa::PANOCProgressInfo::ψ_hat)
        .def_readonly("grad_ψ_hat", &pa::PANOCProgressInfo::grad_ψ_hat)
        .def_readonly("L", &pa::PANOCProgressInfo::L, //
                      "Estimate of Lipschitz constant of objective :math:`L`")
        .def_readonly("γ", &pa::PANOCProgressInfo::γ,
                      "Step size :math:`\\gamma`")
        .def_readonly("τ", &pa::PANOCProgressInfo::τ, //
                      "Line search parameter :math:`\\tau`")
        .def_readonly("ε", &pa::PANOCProgressInfo::ε, //
                      "Tolerance reached :math:`\\varepsilon_k`")
        .def_readonly("Σ", &pa::PANOCProgressInfo::Σ, //
                      "Penalty factor :math:`\\Sigma`")
        .def_readonly("y", &pa::PANOCProgressInfo::y, //
                      "Lagrange multipliers :math:`y`")
        .def_property_readonly(
            "fpr",
            [](const pa::PANOCProgressInfo &p) {
                return std::sqrt(p.norm_sq_p) / p.γ;
            },
            "Fixed-point residual :math:`\\left\\|p\\right\\| / \\gamma`");

    py::class_<pa::StructuredPANOCLBFGSProgressInfo>(
        m, "StructuredPANOCLBFGSProgressInfo",
        "Data passed to the structured PANOC progress callback.\n\n"
        "C++ documentation: :cpp:class:`pa::StructuredPANOCLBFGSProgressInfo`")
        .def_readonly("k", &pa::StructuredPANOCLBFGSProgressInfo::k)
        .def_readonly("x", &pa::StructuredPANOCLBFGSProgressInfo::x)
        .def_readonly("p", &pa::StructuredPANOCLBFGSProgressInfo::p)
        .def_readonly("norm_sq_p",
                      &pa::StructuredPANOCLBFGSProgressInfo::norm_sq_p)
        .def_readonly("x_hat", &pa::StructuredPANOCLBFGSProgressInfo::x_hat)
        .def_readonly("φγ", &pa::StructuredPANOCLBFGSProgressInfo::φγ)
        .def_readonly("ψ", &pa::StructuredPANOCLBFGSProgressInfo::ψ)
        .def_readonly("grad_ψ", &pa::StructuredPANOCLBFGSProgressInfo::grad_ψ)
        .def_readonly("ψ_hat", &pa::StructuredPANOCLBFGSProgressInfo::ψ_hat)
        .def_readonly("grad_ψ_hat",
                      &pa::StructuredPANOCLBFGSProgressInfo::grad_ψ_hat)
        .def_readonly("L", &pa::StructuredPANOCLBFGSProgressInfo::L)
        .def_readonly("γ", &pa::StructuredPANOCLBFGSProgressInfo::γ)
        .def_readonly("τ", &pa::StructuredPANOCLBFGSProgressInfo::τ)
        .def_readonly("ε", &pa::StructuredPANOCLBFGSProgressInfo::ε)
        .def_readonly("Σ", &pa::StructuredPANOCLBFGSProgressInfo::Σ)
        .def_readonly("y", &pa::StructuredPANOCLBFGSProgressInfo::y)
        .def_property_readonly(
            "fpr", [](const pa::StructuredPANOCLBFGSProgressInfo &p) {
                return std::sqrt(p.norm_sq_p) / p.γ;
            });

    py::class_<pa::PolymorphicPANOCSolver,
               std::shared_ptr<pa::PolymorphicPANOCSolver>,
               pa::PolymorphicInnerSolverBase>(
        m, "PANOCSolver", "C++ documentation: :cpp:class:`pa::PANOCSolver`")
        .def(py::init([] {
            return std::make_shared<pa::PolymorphicPANOCSolver>(
                pa::PANOCSolver<pa::PolymorphicPANOCDirectionBase>{
                    pa::PANOCParams{},
                    std::static_pointer_cast<pa::PolymorphicPANOCDirectionBase>(
                        std::make_shared<pa::PolymorphicLBFGSDirection>(
                            pa::LBFGSParams{}))});
        }))
        .def(py::init(PolymorphicPANOCConstructor< //
                      pa::PolymorphicLBFGSDirection>()),
             "panoc_params"_a, "lbfgs_direction"_a)
        .def(py::init(PolymorphicPANOCConversion< //
                      pa::PolymorphicLBFGSDirection, pa::LBFGSParams>()),
             "panoc_params"_a, "lbfgs_params"_a)
        .def(py::init(PolymorphicPANOCConstructor< //
                      pa::PolymorphicPANOCDirectionTrampoline>()),
             "panoc_params"_a, "direction"_a)
        .def(
            "set_progress_callback",
            &pa::PolymorphicPANOCSolver::set_progress_callback, "callback"_a,
            "Attach a callback that is called on each iteration of the solver.")
        .def("__call__",
             pa::InnerSolverCallWrapper<pa::PolymorphicPANOCSolver>(),
             py::call_guard<py::scoped_ostream_redirect,
                            py::scoped_estream_redirect>(),
             "problem"_a, "Σ"_a, "ε"_a, "x"_a,
             "y"_a, //
             "Solve.\n\n"
             ":param problem: Problem to solve\n"
             ":param Σ: Penalty factor\n"
             ":param ε: Desired tolerance\n"
             ":param x: Initial guess\n"
             ":param y: Initial Lagrange multipliers\n\n"
             ":return: * Solution :math:`x`\n"
             "         * Updated Lagrange multipliers :math:`y`\n"
             "         * Slack variable error :math:`g(x) - z`\n"
             "         * Statistics\n\n")
        .def("__str__", &pa::PolymorphicPANOCSolver::get_name)
        .def_property_readonly("params",
                               &pa::PolymorphicPANOCSolver::get_params)
        .def_property_readonly(
            "direction", [](const pa::PolymorphicPANOCSolver &s) {
                return s.innersolver.direction_provider.direction;
            });

    py::class_<pa::PolymorphicPGASolver,
               std::shared_ptr<pa::PolymorphicPGASolver>,
               pa::PolymorphicInnerSolverBase>(
        m, "PGASolver", "C++ documentation: :cpp:class:`pa::PGASolver`")
        .def(py::init<pa::PGAParams>())
        .def(
            "set_progress_callback",
            &pa::PolymorphicPGASolver::set_progress_callback, "callback"_a,
            "Attach a callback that is called on each iteration of the solver.")
        .def("__call__", pa::InnerSolverCallWrapper<pa::PolymorphicPGASolver>(),
             py::call_guard<py::scoped_ostream_redirect,
                            py::scoped_estream_redirect>(),
             "problem"_a, "Σ"_a, "ε"_a, "x"_a,
             "y"_a, //
             "Solve.\n\n"
             ":param problem: Problem to solve\n"
             ":param Σ: Penalty factor\n"
             ":param ε: Desired tolerance\n"
             ":param x: Initial guess\n"
             ":param y: Initial Lagrange multipliers\n\n"
             ":return: * Solution :math:`x`\n"
             "         * Updated Lagrange multipliers :math:`y`\n"
             "         * Slack variable error :math:`g(x) - z`\n"
             "         * Statistics\n\n")
        .def("__str__", &pa::PolymorphicPGASolver::get_name)
        .def_property_readonly("params", &pa::PolymorphicPGASolver::get_params);

    py::class_<pa::PolymorphicGAAPGASolver,
               std::shared_ptr<pa::PolymorphicGAAPGASolver>,
               pa::PolymorphicInnerSolverBase>(
        m, "GAAPGASolver", "C++ documentation: :cpp:class:`pa::GAAPGASolver`")
        .def(py::init<pa::GAAPGAParams>())
        .def(
            "set_progress_callback",
            &pa::PolymorphicGAAPGASolver::set_progress_callback, "callback"_a,
            "Attach a callback that is called on each iteration of the solver.")
        .def("__call__",
             pa::InnerSolverCallWrapper<pa::PolymorphicGAAPGASolver>(),
             py::call_guard<py::scoped_ostream_redirect,
                            py::scoped_estream_redirect>(),
             "problem"_a, "Σ"_a, "ε"_a, "x"_a,
             "y"_a, //
             "Solve.\n\n"
             ":param problem: Problem to solve\n"
             ":param Σ: Penalty factor\n"
             ":param ε: Desired tolerance\n"
             ":param x: Initial guess\n"
             ":param y: Initial Lagrange multipliers\n\n"
             ":return: * Solution :math:`x`\n"
             "         * Updated Lagrange multipliers :math:`y`\n"
             "         * Slack variable error :math:`g(x) - z`\n"
             "         * Statistics\n\n")
        .def("__str__", &pa::PolymorphicGAAPGASolver::get_name)
        .def_property_readonly("params",
                               &pa::PolymorphicGAAPGASolver::get_params);

    py::enum_<pa::PANOCStopCrit>(
        m, "PANOCStopCrit", "C++ documentation: :cpp:enum:`pa::PANOCStopCrit`")
        .value("ApproxKKT", pa::PANOCStopCrit::ApproxKKT)
        .value("ApproxKKT2", pa::PANOCStopCrit::ApproxKKT2)
        .value("ProjGradNorm", pa::PANOCStopCrit::ProjGradNorm)
        .value("ProjGradNorm2", pa::PANOCStopCrit::ProjGradNorm2)
        .value("ProjGradUnitNorm", pa::PANOCStopCrit::ProjGradUnitNorm)
        .value("ProjGradUnitNorm2", pa::PANOCStopCrit::ProjGradUnitNorm2)
        .value("FPRNorm", pa::PANOCStopCrit::FPRNorm)
        .value("FPRNorm2", pa::PANOCStopCrit::FPRNorm2)
        .export_values();

    py::class_<pa::StructuredPANOCLBFGSParams>(
        m, "StructuredPANOCLBFGSParams",
        "C++ documentation: :cpp:class:`pa::StructuredPANOCLBFGSParams`")
        .def(py::init(&kwargs_to_struct<pa::StructuredPANOCLBFGSParams>))
        .def("to_dict", &struct_to_dict<pa::StructuredPANOCLBFGSParams>)

        .def_readwrite("Lipschitz", &pa::StructuredPANOCLBFGSParams::Lipschitz)
        .def_readwrite("max_iter", &pa::StructuredPANOCLBFGSParams::max_iter)
        .def_readwrite("max_time", &pa::StructuredPANOCLBFGSParams::max_time)
        .def_readwrite("τ_min", &pa::StructuredPANOCLBFGSParams::τ_min)
        .def_readwrite("L_min", &pa::StructuredPANOCLBFGSParams::L_min)
        .def_readwrite("L_max", &pa::StructuredPANOCLBFGSParams::L_max)
        .def_readwrite("nonmonotone_linesearch",
                       &pa::StructuredPANOCLBFGSParams::nonmonotone_linesearch)
        .def_readwrite("stop_crit", &pa::StructuredPANOCLBFGSParams::stop_crit)
        .def_readwrite("max_no_progress",
                       &pa::StructuredPANOCLBFGSParams::max_no_progress)
        .def_readwrite("print_interval",
                       &pa::StructuredPANOCLBFGSParams::print_interval)
        .def_readwrite("quadratic_upperbound_tolerance_factor",
                       &pa::StructuredPANOCLBFGSParams::
                           quadratic_upperbound_tolerance_factor)
        .def_readwrite(
            "update_lipschitz_in_linesearch",
            &pa::StructuredPANOCLBFGSParams::update_lipschitz_in_linesearch)
        .def_readwrite(
            "alternative_linesearch_cond",
            &pa::StructuredPANOCLBFGSParams::alternative_linesearch_cond)
        .def_readwrite(
            "hessian_vec_finited_differences",
            &pa::StructuredPANOCLBFGSParams::hessian_vec_finited_differences)
        .def_readwrite("full_augmented_hessian",
                       &pa::StructuredPANOCLBFGSParams::full_augmented_hessian)
        .def_readwrite("lbfgs_stepsize",
                       &pa::StructuredPANOCLBFGSParams::lbfgs_stepsize);

    py::class_<pa::PolymorphicStructuredPANOCLBFGSSolver,
               std::shared_ptr<pa::PolymorphicStructuredPANOCLBFGSSolver>,
               pa::PolymorphicInnerSolverBase>(
        m, "StructuredPANOCLBFGSSolver",
        "C++ documentation: :cpp:class:`pa::StructuredPANOCLBFGSSolver`")
        .def(py::init([] {
            return std::make_shared<pa::PolymorphicStructuredPANOCLBFGSSolver>(
                pa::StructuredPANOCLBFGSSolver{
                    pa::StructuredPANOCLBFGSParams{},
                    pa::LBFGSParams{},
                });
        }))
        .def(py::init<pa::StructuredPANOCLBFGSParams, pa::LBFGSParams>(),
             "panoc_params"_a, "lbfgs_params"_a)
        .def(
            "set_progress_callback",
            &pa::PolymorphicStructuredPANOCLBFGSSolver::set_progress_callback,
            "callback"_a,
            "Attach a callback that is called on each iteration of the solver.")
        .def("__call__",
             pa::InnerSolverCallWrapper<
                 pa::PolymorphicStructuredPANOCLBFGSSolver>(),
             py::call_guard<py::scoped_ostream_redirect,
                            py::scoped_estream_redirect>(),
             "problem"_a, "Σ"_a, "ε"_a, "x"_a,
             "y"_a, //
             "Solve.\n\n"
             ":param problem: Problem to solve\n"
             ":param Σ: Penalty factor\n"
             ":param ε: Desired tolerance\n"
             ":param x: Initial guess\n"
             ":param y: Initial Lagrange multipliers\n\n"
             ":return: * Solution :math:`x`\n"
             "         * Updated Lagrange multipliers :math:`y`\n"
             "         * Slack variable error :math:`g(x) - z`\n"
             "         * Statistics\n\n")
        .def("__str__", &pa::PolymorphicStructuredPANOCLBFGSSolver::get_name)
        .def_property_readonly(
            "params", &pa::PolymorphicStructuredPANOCLBFGSSolver::get_params);

    py::class_<pa::ALMParams>(m, "ALMParams",
                              "C++ documentation: :cpp:class:`pa::ALMParams`")
        .def(py::init())
        .def(py::init(&kwargs_to_struct<pa::ALMParams>))
        .def("to_dict", &struct_to_dict<pa::ALMParams>)
        .def_readwrite("ε", &pa::ALMParams::ε)
        .def_readwrite("δ", &pa::ALMParams::δ)
        .def_readwrite("Δ", &pa::ALMParams::Δ)
        .def_readwrite("Δ_lower", &pa::ALMParams::Δ_lower)
        .def_readwrite("Σ_0", &pa::ALMParams::Σ₀)
        .def_readwrite("σ_0", &pa::ALMParams::σ₀)
        .def_readwrite("Σ_0_lower", &pa::ALMParams::Σ₀_lower)
        .def_readwrite("ε_0", &pa::ALMParams::ε₀)
        .def_readwrite("ε_0_increase", &pa::ALMParams::ε₀_increase)
        .def_readwrite("ρ", &pa::ALMParams::ρ)
        .def_readwrite("ρ_increase", &pa::ALMParams::ρ_increase)
        .def_readwrite("θ", &pa::ALMParams::θ)
        .def_readwrite("M", &pa::ALMParams::M)
        .def_readwrite("Σ_max", &pa::ALMParams::Σ_max)
        .def_readwrite("Σ_min", &pa::ALMParams::Σ_min)
        .def_readwrite("max_iter", &pa::ALMParams::max_iter)
        .def_readwrite("max_time", &pa::ALMParams::max_time)
        .def_readwrite("max_num_initial_retries",
                       &pa::ALMParams::max_num_initial_retries)
        .def_readwrite("max_num_retries", &pa::ALMParams::max_num_retries)
        .def_readwrite("max_total_num_retries",
                       &pa::ALMParams::max_total_num_retries)
        .def_readwrite("print_interval", &pa::ALMParams::print_interval)
        .def_readwrite("preconditioning", &pa::ALMParams::preconditioning)
        .def_readwrite("single_penalty_factor",
                       &pa::ALMParams::single_penalty_factor);

    py::class_<pa::PolymorphicALMSolver>(
        m, "ALMSolver",
        "Main augmented Lagrangian solver.\n\n"
        "C++ documentation: :cpp:class:`pa::ALMSolver`")
        // Default constructor
        .def(py::init([] {
                 return pa::PolymorphicALMSolver{
                     pa::ALMParams{},
                     std::static_pointer_cast<pa::PolymorphicInnerSolverBase>(
                         std::make_shared<
                             pa::PolymorphicStructuredPANOCLBFGSSolver>(
                             pa::StructuredPANOCLBFGSParams{},
                             pa::LBFGSParams{})),
                 };
             }),
             "Build an ALM solver using Structured PANOC as inner solver.")
        // Params and solver
        .def(py::init(PolymorphicALMConstructor<pa::PolymorphicPANOCSolver>()),
             "alm_params"_a, "panoc_solver"_a,
             "Build an ALM solver using PANOC as inner solver.")
        .def(py::init(PolymorphicALMConstructor<pa::PolymorphicPGASolver>()),
             "alm_params"_a, "pga_solver"_a,
             "Build an ALM solver using the projected gradient algorithm as "
             "inner solver.")
        .def(py::init(PolymorphicALMConstructor<
                      pa::PolymorphicStructuredPANOCLBFGSSolver>()),
             "alm_params"_a, "structuredpanoc_solver"_a,
             "Build an ALM solver using Structured PANOC as inner solver.")
        .def(py::init(PolymorphicALMConstructor<
                      pa::PolymorphicInnerSolverTrampoline>()),
             "alm_params"_a, "inner_solver"_a,
             "Build an ALM solver using a custom inner solver.")
        // Only solver (default params)
        .def(py::init(PolymorphicALMConstructorDefaultParams<
                      pa::PolymorphicPANOCSolver>()),
             "panoc_solver"_a,
             "Build an ALM solver using PANOC as inner solver.")
        .def(py::init(PolymorphicALMConstructorDefaultParams<
                      pa::PolymorphicPGASolver>()),
             "pga_solver"_a,
             "Build an ALM solver using the projected gradient algorithm as "
             "inner solver.")
        .def(py::init(PolymorphicALMConstructorDefaultParams<
                      pa::PolymorphicStructuredPANOCLBFGSSolver>()),
             "structuredpanoc_solver"_a,
             "Build an ALM solver using Structured PANOC as inner solver.")
        .def(py::init(PolymorphicALMConstructorDefaultParams<
                      pa::PolymorphicInnerSolverTrampoline>()),
             "inner_solver"_a,
             "Build an ALM solver using a custom inner solver.")
        // Only params (default solver)
        .def(py::init([](const pa::ALMParams &params) {
                 return pa::PolymorphicALMSolver{
                     params,
                     std::static_pointer_cast<pa::PolymorphicInnerSolverBase>(
                         std::make_shared<
                             pa::PolymorphicStructuredPANOCLBFGSSolver>(
                             pa::StructuredPANOCLBFGSParams{},
                             pa::LBFGSParams{})),
                 };
             }),
             "alm_params"_a,
             "Build an ALM solver using Structured PANOC as inner solver.")
        // Other functions and properties
        .def_property_readonly("inner_solver",
                               [](const pa::PolymorphicALMSolver &s) {
                                   return s.inner_solver.solver;
                               })
        .def(
            "__call__",
            [](pa::PolymorphicALMSolver &solver, const pa::Problem &p,
               std::optional<pa::vec> x, std::optional<pa::vec> y)
                -> std::tuple<pa::vec, pa::vec, py::dict> {
                if (!x)
                    x = pa::vec::Zero(p.n);
                else if (x->size() != p.n)
                    throw std::invalid_argument(
                        "Length of x does not match problem size problem.n");
                if (!y)
                    y = pa::vec::Zero(p.m);
                else if (y->size() != p.m)
                    throw std::invalid_argument(
                        "Length of y does not match problem size problem.m");
                if (p.C.lowerbound.size() != p.n)
                    throw std::invalid_argument(
                        "Length of problem.C.lowerbound does not match problem "
                        "size problem.n");
                if (p.C.upperbound.size() != p.n)
                    throw std::invalid_argument(
                        "Length of problem.C.upperbound does not match problem "
                        "size problem.n");
                if (p.D.lowerbound.size() != p.m)
                    throw std::invalid_argument(
                        "Length of problem.D.lowerbound does not match problem "
                        "size problem.m");
                if (p.D.upperbound.size() != p.m)
                    throw std::invalid_argument(
                        "Length of problem.D.upperbound does not match problem "
                        "size problem.m");

                auto stats = solver(p, *y, *x);
                return std::make_tuple(std::move(*x), std::move(*y),
                                       stats_to_dict(stats));
            },
            "problem"_a, "x"_a = std::nullopt, "y"_a = std::nullopt,
            py::call_guard<py::scoped_ostream_redirect,
                           py::scoped_estream_redirect>(),
            "Solve.\n\n"
            ":param problem: Problem to solve.\n"
            ":param y: Initial guess for Lagrange multipliers :math:`y`\n"
            ":param x: Initial guess for decision variables :math:`x`\n\n"
            ":return: * Lagrange multipliers :math:`y` at the solution\n"
            "         * Solution :math:`x`\n"
            "         * Statistics\n\n")
        .def("__str__", &pa::PolymorphicALMSolver::get_name)
        .def_property_readonly("params", &pa::PolymorphicALMSolver::get_params);

    constexpr auto panoc = [](std::function<pa::real_t(pa::crvec)> ψ,
                              std::function<pa::vec(pa::crvec)> grad_ψ,
                              const pa::Box &C, std::optional<pa::vec> x0,
                              pa::real_t ε, const pa::PANOCParams &params,
                              const pa::LBFGSParams &lbfgs_params) {
        auto n = C.lowerbound.size();
        if (C.upperbound.size() != n)
            throw std::invalid_argument("Length of C.upperbound does not "
                                        "match length of C.lowerbound");
        if (!x0)
            x0 = pa::vec::Zero(n);
        else if (x0->size() != n)
            throw std::invalid_argument(
                "Length of x does not match length of bounds C");
        auto grad_ψ_ = [&](pa::crvec x, pa::rvec gr) {
            auto &&t = grad_ψ(x);
            if (t.size() != x.size())
                throw std::runtime_error("Invalid grad_ψ dimension");
            gr = std::move(t);
        };
        auto stats =
            pa::panoc<pa::LBFGS>(ψ, grad_ψ_, C, *x0, ε, params, {lbfgs_params});
        return std::make_tuple(std::move(*x0), stats_to_dict(stats));
    };

    m.def("panoc", panoc, "ψ"_a, "grad_ψ"_a, "C"_a, "x0"_a = std::nullopt,
          "ε"_a = 1e-8, "params"_a = pa::PANOCParams{},
          "lbfgs_params"_a = pa::LBFGSParams{});

#if !PANOCPY_HAVE_CASADI
    auto load_CasADi_problem = [](const char *, unsigned, unsigned,
                                  bool) -> pa::Problem {
        throw std::runtime_error(
            "This version of panocpy was compiled without CasADi support");
    };
    auto load_CasADi_problem_with_param = [](const char *, unsigned, unsigned,
                                             unsigned,
                                             bool) -> pa::ProblemWithParam {
        throw std::runtime_error(
            "This version of panocpy was compiled without CasADi support");
    };
#else
    using pa::load_CasADi_problem;
    using pa::load_CasADi_problem_with_param;
#endif

    m.def("load_casadi_problem", load_CasADi_problem, "so_name"_a, "n"_a = 0,
          "m"_a = 0, "second_order"_a = false,
          "Load a compiled CasADi problem without parameters.\n\n"
          "C++ documentation: :cpp:func:`pa::load_CasADi_problem`");
    m.def("load_casadi_problem_with_param", load_CasADi_problem_with_param,
          "so_name"_a, "n"_a = 0, "m"_a = 0, "p"_a = 0,
          "second_order"_a = false,
          "Load a compiled CasADi problem with parameters.\n\n"
          "C++ documentation: :cpp:func:`pa::load_CasADi_problem_with_param`");
}
