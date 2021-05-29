/**
 * @file
 * This file defines all Python bindings.
 */

#include <panoc-alm/inner/decl/panoc-stop-crit.hpp>
#include <panoc-alm/inner/directions/lbfgs.hpp>
#include <panoc-alm/inner/panoc.hpp>
#include <panoc-alm/inner/second-order-panoc-lbfgs.hpp>
#include <panoc-alm/interop/casadi/CasADiLoader.hpp>
#include <panoc-alm/util/solverstatus.hpp>

#include <pybind11/attr.h>
#include <pybind11/cast.h>
#include <pybind11/chrono.h>
#include <pybind11/detail/common.h>
#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include <pybind11/iostream.h>
#include <pybind11/pybind11.h>
#include <pybind11/pytypes.h>

#include <memory>
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

struct ostream_redirect {
    py::scoped_ostream_redirect stdout_stream{
        std::cout,                                 // std::ostream&
        py::module_::import("sys").attr("stdout"), // Python output
    };
    py::scoped_ostream_redirect stderr_stream{
        std::cerr,                                 // std::ostream&
        py::module_::import("sys").attr("stderr"), // Python output
    };
};

PYBIND11_MODULE(PANOCPY_MODULE_NAME, m) {
    using py::operator""_a;

    m.doc() = "PANOC+ALM solvers"; // TODO

    py::class_<pa::Box>(m, "Box")
        .def_readwrite("upperbound", &pa::Box::upperbound)
        .def_readwrite("lowerbound", &pa::Box::lowerbound);

    py::class_<pa::Problem>(m, "Problem")
        // .def(py::init())
        .def(py::init<unsigned, unsigned>(), "n"_a, "m"_a)
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
                      prob_setter_hess_L_prod());

    py::class_<ProblemWithParam, pa::Problem>(m, "ProblemWithParam")
        .def(py::init<unsigned, unsigned>(), "n"_a, "m"_a)
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
            "param", py::overload_cast<>(&ProblemWithParam::get_param),
            py::overload_cast<pa::crvec>(&ProblemWithParam::set_param));

    py::class_<pa::PolymorphicPANOCDirectionBase,
               std::shared_ptr<pa::PolymorphicPANOCDirectionBase>,
               pa::PolymorphicPANOCDirectionTrampoline>(m, "PANOCDirection")
        .def(py::init<>())
        .def("initialize", &pa::PolymorphicPANOCDirectionBase::initialize)
        .def("update", &pa::PolymorphicPANOCDirectionBase::update)
        .def("apply", &pa::PolymorphicPANOCDirectionBase::apply_ret)
        .def("changed_γ", &pa::PolymorphicPANOCDirectionBase::changed_γ)
        .def("reset", &pa::PolymorphicPANOCDirectionBase::reset)
        .def("get_name", &pa::PolymorphicPANOCDirectionBase::get_name)
        .def("__str__", &pa::PolymorphicPANOCDirectionBase::get_name);

    py::class_<pa::PolymorphicLBFGSDirection,
               std::shared_ptr<pa::PolymorphicLBFGSDirection>,
               pa::PolymorphicPANOCDirectionBase>(m, "LBFGSDirection")
        .def(py::init<pa::LBFGSParams>())
        .def("initialize", &pa::PolymorphicLBFGSDirection::initialize)
        .def("update", &pa::PolymorphicLBFGSDirection::update)
        .def("apply", &pa::PolymorphicLBFGSDirection::apply_ret)
        .def("changed_γ", &pa::PolymorphicLBFGSDirection::changed_γ)
        .def("reset", &pa::PolymorphicLBFGSDirection::reset)
        .def("get_name", &pa::PolymorphicLBFGSDirection::get_name)
        .def("__str__", &pa::PolymorphicLBFGSDirection::get_name);

    using paLBFGSParamCBFGS = decltype(pa::LBFGSParams::cbfgs);
    py::class_<paLBFGSParamCBFGS>(m, "LBFGSParamsCBFGS")
        .def(py::init())
        .def(py::init(&kwargs_to_struct<paLBFGSParamCBFGS>))
        .def_readwrite("α", &paLBFGSParamCBFGS::α)
        .def_readwrite("ϵ", &paLBFGSParamCBFGS::ϵ);

    py::class_<pa::LBFGSParams>(m, "LBFGSParams")
        .def(py::init())
        .def(py::init(&kwargs_to_struct<pa::LBFGSParams>))
        .def_readwrite("memory", &pa::LBFGSParams::memory)
        .def_readwrite("cbfgs", &pa::LBFGSParams::cbfgs)
        .def_readwrite("rescale_when_γ_changes",
                       &pa::LBFGSParams::rescale_when_γ_changes);

    py::enum_<pa::LBFGSStepSize>(m, "LBFGSStepsize")
        .value("BasedOnGradientStepSize",
               pa::LBFGSStepSize::BasedOnGradientStepSize)
        .value("BasedOnCurvature", pa::LBFGSStepSize::BasedOnCurvature)
        .export_values();

    using paPANOCParamsLipschitz = decltype(pa::PANOCParams::Lipschitz);
    py::class_<paPANOCParamsLipschitz>(m, "PANOCParamsLipschitz")
        .def(py::init())
        .def(py::init(&kwargs_to_struct<paPANOCParamsLipschitz>))
        .def_readwrite("L_0", &paPANOCParamsLipschitz::L₀)
        .def_readwrite("ε", &paPANOCParamsLipschitz::ε)
        .def_readwrite("δ", &paPANOCParamsLipschitz::δ)
        .def_readwrite("Lγ_factor", &paPANOCParamsLipschitz::Lγ_factor);

    py::class_<pa::PANOCParams>(m, "PANOCParams")
        .def(py::init())
        .def(py::init(&kwargs_to_struct<pa::PANOCParams>))
        .def_readwrite("Lipschitz", &pa::PANOCParams::Lipschitz)
        .def_readwrite("max_iter", &pa::PANOCParams::max_iter)
        .def_readwrite("max_time", &pa::PANOCParams::max_time)
        .def_readwrite("τ_min", &pa::PANOCParams::τ_min)
        .def_readwrite("γ_min", &pa::PANOCParams::γ_min)
        .def_readwrite("max_no_progress", &pa::PANOCParams::max_no_progress)
        .def_readwrite("print_interval", &pa::PANOCParams::print_interval)
        .def_readwrite("quadratic_upperbound_tolerance_factor",
                       &pa::PANOCParams::quadratic_upperbound_tolerance_factor)
        .def_readwrite("update_lipschitz_in_linesearch",
                       &pa::PANOCParams::update_lipschitz_in_linesearch)
        .def_readwrite("alternative_linesearch_cond",
                       &pa::PANOCParams::alternative_linesearch_cond)
        .def_readwrite("lbfgs_stepsize", &pa::PANOCParams::lbfgs_stepsize);

    py::enum_<pa::SolverStatus>(m, "SolverStatus", py::arithmetic(),
                                "Solver status")
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
        .def("get_name", &pa::PolymorphicInnerSolverBase::get_name);

    py::class_<pa::PolymorphicPANOCSolver,
               std::shared_ptr<pa::PolymorphicPANOCSolver>,
               pa::PolymorphicInnerSolverBase>(m, "PANOCSolver")
        .def(py::init(PolymorphicPANOCConstructor< //
                      pa::PolymorphicLBFGSDirection>()))
        .def(py::init(PolymorphicPANOCConversion< //
                      pa::PolymorphicLBFGSDirection,
                      pa::LBFGSParams>()))
        .def(py::init(PolymorphicPANOCConstructor< //
                      pa::PolymorphicPANOCDirectionTrampoline>()))
        .def("__call__",
             pa::InnerSolverCallWrapper<pa::PolymorphicPANOCSolver>(),
             py::call_guard<py::scoped_ostream_redirect,
                            py::scoped_estream_redirect>())
        .def("__str__", &pa::PolymorphicPANOCSolver::get_name);

    py::enum_<pa::PANOCStopCrit>(m, "PANOCStopCrit")
        .value("ApproxKKT", pa::PANOCStopCrit::ApproxKKT)
        .value("ProjGradNorm", pa::PANOCStopCrit::ProjGradNorm)
        .value("ProjGradUnitNorm", pa::PANOCStopCrit::ProjGradUnitNorm)
        .value("FPRNorm", pa::PANOCStopCrit::FPRNorm)
        .export_values();

    py::class_<pa::SecondOrderPANOCLBFGSParams>(m,
                                                "SecondOrderPANOCLBFGSParams")
        .def(py::init(&kwargs_to_struct<pa::SecondOrderPANOCLBFGSParams>));

    py::class_<pa::PolymorphicSecondOrderPANOCLBFGSSolver,
               std::shared_ptr<pa::PolymorphicSecondOrderPANOCLBFGSSolver>,
               pa::PolymorphicInnerSolverBase>(m, "SecondOrderPANOCLBFGSSolver")
        .def(py::init<pa::SecondOrderPANOCLBFGSParams, pa::LBFGSParams>())
        .def("__call__",
             pa::InnerSolverCallWrapper<
                 pa::PolymorphicSecondOrderPANOCLBFGSSolver>(),
             py::call_guard<py::scoped_ostream_redirect,
                            py::scoped_estream_redirect>())
        .def("__str__", &pa::PolymorphicSecondOrderPANOCLBFGSSolver::get_name);

    py::class_<pa::ALMParams>(m, "ALMParams")
        .def(py::init())
        .def(py::init(&kwargs_to_struct<pa::ALMParams>))
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

    py::class_<pa::PolymorphicALMSolver>(m, "ALMSolver")
        .def(py::init(PolymorphicALMConstructor<pa::PolymorphicPANOCSolver>()))
        .def(py::init(PolymorphicALMConstructor<
                      pa::PolymorphicSecondOrderPANOCLBFGSSolver>()))
        .def(py::init(
            PolymorphicALMConstructor<pa::PolymorphicInnerSolverTrampoline>()))
        .def(
            "__call__",
            [](pa::PolymorphicALMSolver &solver, const pa::Problem &p,
               pa::vec y, pa::vec x) -> std::tuple<pa::vec, pa::vec, py::dict> {
                auto stats = solver(p, y, x);
                return std::make_tuple(std::move(y), std::move(x),
                                       stats_to_dict(stats));
            },
            py::call_guard<py::scoped_ostream_redirect,
                           py::scoped_estream_redirect>());

    m.def("load_casadi_problem", &load_CasADi_problem, "so_name"_a, "n"_a,
          "m"_a);
    m.def("load_casadi_problem_with_param", &load_CasADi_problem_with_param,
          "so_name"_a, "n"_a, "m"_a);
}
