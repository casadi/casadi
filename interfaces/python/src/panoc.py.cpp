#include <alpaqa/util/quadmath/quadmath.hpp>

#include <pybind11/chrono.h>
#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
using namespace py::literals;

#include <alpaqa/implementation/inner/panoc.tpp>
#include <alpaqa/inner/directions/panoc/lbfgs.hpp>
#include <alpaqa/inner/directions/panoc/structured-lbfgs.hpp>
#include <alpaqa/inner/panoc.hpp>
#include <alpaqa/util/check-dim.hpp>

#include "async.hpp"
#include "copy.hpp"
#include "inner-solve.hpp"
#include "member.hpp"
#include "params/params.hpp"
#include "stats-to-dict.hpp"
#include "type-erased-inner-solver.hpp"
#include "type-erased-panoc-direction.hpp"

template <alpaqa::Config Conf>
void register_panoc(py::module_ &m) {
    USING_ALPAQA_CONFIG(Conf);

    using TEProblem         = alpaqa::TypeErasedProblem<config_t>;
    using InnerSolver       = alpaqa::TypeErasedInnerSolver<config_t, TEProblem>;
    using InnerSolveOptions = alpaqa::InnerSolveOptions<config_t>;
    register_dataclass<InnerSolveOptions>(m, "InnerSolveOptions");
    py::implicitly_convertible<py::dict, InnerSolveOptions>();

    // ----------------------------------------------------------------------------------------- //
    using TypeErasedPANOCDirection = alpaqa::TypeErasedPANOCDirection<Conf>;
    using LBFGSParams              = alpaqa::LBFGSParams<config_t>;
    using StructuredLBFGSDir       = alpaqa::StructuredLBFGSDirection<config_t>;
    using StrucLBFGSDirParams      = alpaqa::StructuredLBFGSDirectionParams<config_t>;

    // ----------------------------------------------------------------------------------------- //
    using LipschitzEstimateParams = alpaqa::LipschitzEstimateParams<config_t>;
    register_dataclass<LipschitzEstimateParams>(
        m, "LipschitzEstimateParams",
        "C++ documentation: :cpp:class:`alpaqa::LipschitzEstimateParams`");

    using PANOCParams = alpaqa::PANOCParams<config_t>;
    register_dataclass<PANOCParams>(m, "PANOCParams",
                                    "C++ documentation: :cpp:class:`alpaqa::PANOCParams`");

    // ----------------------------------------------------------------------------------------- //
    using PANOCProgressInfo = alpaqa::PANOCProgressInfo<config_t>;
    py::class_<PANOCProgressInfo>(m, "PANOCProgressInfo",
                                  "Data passed to the PANOC progress callback.\n\n"
                                  "C++ documentation: :cpp:class:`alpaqa::PANOCProgressInfo`")
        // clang-format off
        .def_readonly("k", &PANOCProgressInfo::k, "Iteration")
        .def_readonly("x", &PANOCProgressInfo::x, "Decision variable :math:`x`")
        .def_readonly("p", &PANOCProgressInfo::p, "Projected gradient step :math:`p`")
        .def_readonly("norm_sq_p", &PANOCProgressInfo::norm_sq_p, ":math:`\\left\\|p\\right\\|^2`")
        .def_readonly("x̂", &PANOCProgressInfo::x̂, "Decision variable after projected gradient step :math:`\\hat x`")
        .def_readonly("φγ", &PANOCProgressInfo::φγ, "Forward-backward envelope :math:`\\varphi_\\gamma(x)`")
        .def_readonly("ψ", &PANOCProgressInfo::ψ, "Objective value :math:`\\psi(x)`")
        .def_readonly("grad_ψ", &PANOCProgressInfo::grad_ψ, "Gradient of objective :math:`\\nabla\\psi(x)`")
        .def_readonly("ψ_hat", &PANOCProgressInfo::ψ_hat, "Objective at x̂ :math:`\\psi(\\hat x)`")
        .def_readonly("grad_ψ_hat", &PANOCProgressInfo::grad_ψ_hat, "Gradient of objective at x̂ :math:`\\nabla\\psi(\\hat x)`")
        .def_readonly("q", &PANOCProgressInfo::q, "Previous quasi-Newton step :math:`\\nabla\\psi(\\hat x)`")
        .def_readonly("L", &PANOCProgressInfo::L, "Estimate of Lipschitz constant of objective :math:`L`")
        .def_readonly("γ", &PANOCProgressInfo::γ, "Step size :math:`\\gamma`")
        .def_readonly("τ", &PANOCProgressInfo::τ, "Previous line search parameter :math:`\\tau`")
        .def_readonly("ε", &PANOCProgressInfo::ε, "Tolerance reached :math:`\\varepsilon_k`")
        .def_readonly("Σ", &PANOCProgressInfo::Σ, "Penalty factor :math:`\\Sigma`")
        .def_readonly("y", &PANOCProgressInfo::y, "Lagrange multipliers :math:`y`")
        .def_property_readonly("problem", member_ptr<&PANOCProgressInfo::problem>(), "Problem being solved")
        .def_property_readonly("params", member_ptr<&PANOCProgressInfo::params>(), "Solver parameters")
        // clang-format on
        .def_property_readonly(
            "fpr", [](const PANOCProgressInfo &p) { return std::sqrt(p.norm_sq_p) / p.γ; },
            "Fixed-point residual :math:`\\left\\|p\\right\\| / \\gamma`");

    // Solve without ALM
    using PANOCSolver = alpaqa::PANOCSolver<TypeErasedPANOCDirection>;
    using Problem     = typename PANOCSolver::Problem;

    py::class_<PANOCSolver> panoc_solver(m, "PANOCSolver",
                                         "C++ documentation: :cpp:class:`alpaqa::PANOCSolver`");
    default_copy_methods(panoc_solver);
    panoc_solver
        // Constructors
        .def(py::init([](params_or_dict<PANOCParams> params,
                         params_or_dict<LBFGSParams> lbfgs_params,
                         params_or_dict<StrucLBFGSDirParams> direction_params) {
                 return PANOCSolver{var_kwargs_to_struct(params),
                                    alpaqa::erase_direction_with_params_dict<StructuredLBFGSDir>(
                                        var_kwargs_to_struct(lbfgs_params),
                                        var_kwargs_to_struct(direction_params))};
             }),
             "panoc_params"_a = py::dict{}, "lbfgs_params"_a = py::dict{},
             "direction_params"_a = py::dict{},
             "Create a PANOC solver using structured L-BFGS directions.")
        .def(py::init(
                 [](params_or_dict<PANOCParams> params, const TypeErasedPANOCDirection &direction) {
                     return PANOCSolver{var_kwargs_to_struct(params),
                                        typename PANOCSolver::Direction{direction}};
                 }),
             "panoc_params"_a, "direction"_a, "Create a PANOC solver using a custom direction.")
        .def_property_readonly("direction", member_ref<&PANOCSolver::direction>());
    register_inner_solver_methods<PANOCSolver, Problem, InnerSolver>(panoc_solver);
}

template void register_panoc<alpaqa::EigenConfigd>(py::module_ &);
template void register_panoc<alpaqa::EigenConfigf>(py::module_ &);
template void register_panoc<alpaqa::EigenConfigl>(py::module_ &);
#ifdef ALPAQA_WITH_QUAD_PRECISION
template void register_panoc<alpaqa::EigenConfigq>(py::module_ &);
#endif

// TODO: Why is this needed?
template class alpaqa::PANOCSolver<alpaqa::TypeErasedPANOCDirection<alpaqa::EigenConfigf>>;
template class alpaqa::PANOCSolver<alpaqa::TypeErasedPANOCDirection<alpaqa::EigenConfigd>>;
template class alpaqa::PANOCSolver<alpaqa::TypeErasedPANOCDirection<alpaqa::EigenConfigl>>;
#ifdef ALPAQA_WITH_QUAD_PRECISION
template class alpaqa::PANOCSolver<alpaqa::TypeErasedPANOCDirection<alpaqa::EigenConfigq>>;
#endif