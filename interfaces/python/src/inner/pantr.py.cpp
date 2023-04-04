#include <alpaqa/util/quadmath/quadmath.hpp>

#include <pybind11/chrono.h>
#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
using namespace py::literals;

#include <alpaqa/implementation/inner/pantr.tpp>
#include <alpaqa/inner/directions/pantr/newton-tr.hpp>
#include <alpaqa/util/check-dim.hpp>

#include "type-erased-tr-direction.hpp"
#include <dict/stats-to-dict.hpp>
#include <inner/inner-solve.hpp>
#include <inner/type-erased-inner-solver.hpp>
#include <params/params.hpp>
#include <util/async.hpp>
#include <util/copy.hpp>
#include <util/member.hpp>

template <alpaqa::Config Conf>
void register_pantr(py::module_ &m) {
    USING_ALPAQA_CONFIG(Conf);

    using TEProblem   = alpaqa::TypeErasedProblem<config_t>;
    using InnerSolver = alpaqa::TypeErasedInnerSolver<config_t, TEProblem>;

    // ----------------------------------------------------------------------------------------- //
    using TypeErasedTRDirection = alpaqa::TypeErasedTRDirection<Conf>;
    using NewtonTRDir           = alpaqa::NewtonTRDirection<config_t>;
    using SteihaugCGParams      = alpaqa::SteihaugCGParams<config_t>;
    using NewtonTRDirParams     = alpaqa::NewtonTRDirectionParams<config_t>;

    // ----------------------------------------------------------------------------------------- //
    using PANTRParams = alpaqa::PANTRParams<config_t>;
    register_dataclass<PANTRParams>(m, "PANTRParams",
                                    "C++ documentation: :cpp:class:`alpaqa::PANTRParams`");

    // ----------------------------------------------------------------------------------------- //
    using PANTRProgressInfo = alpaqa::PANTRProgressInfo<config_t>;
    py::class_<PANTRProgressInfo>(m, "PANTRProgressInfo",
                                  "Data passed to the PANTR progress callback.\n\n"
                                  "C++ documentation: :cpp:class:`alpaqa::PANTRProgressInfo`")
        // clang-format off
        .def_readonly("k", &PANTRProgressInfo::k, "Iteration")
        .def_readonly("status", &PANTRProgressInfo::status, "Current solver status")
        .def_readonly("x", &PANTRProgressInfo::x, "Decision variable :math:`x`")
        .def_readonly("p", &PANTRProgressInfo::p, "Projected gradient step :math:`p`")
        .def_readonly("norm_sq_p", &PANTRProgressInfo::norm_sq_p, ":math:`\\left\\|p\\right\\|^2`")
        .def_readonly("x̂", &PANTRProgressInfo::x̂, "Decision variable after projected gradient step :math:`\\hat x`")
        .def_readonly("φγ", &PANTRProgressInfo::φγ, "Forward-backward envelope :math:`\\varphi_\\gamma(x)`")
        .def_readonly("ψ", &PANTRProgressInfo::ψ, "Objective value :math:`\\psi(x)`")
        .def_readonly("grad_ψ", &PANTRProgressInfo::grad_ψ, "Gradient of objective :math:`\\nabla\\psi(x)`")
        .def_readonly("ψ_hat", &PANTRProgressInfo::ψ_hat, "Objective at x̂ :math:`\\psi(\\hat x)`")
        .def_readonly("grad_ψ_hat", &PANTRProgressInfo::grad_ψ_hat, "Gradient of objective at x̂ :math:`\\nabla\\psi(\\hat x)`")
        .def_readonly("q", &PANTRProgressInfo::q, "Previous quasi-Newton step :math:`\\nabla\\psi(\\hat x)`")
        .def_readonly("L", &PANTRProgressInfo::L, "Estimate of Lipschitz constant of objective :math:`L`")
        .def_readonly("γ", &PANTRProgressInfo::γ, "Step size :math:`\\gamma`")
        .def_readonly("Δ", &PANTRProgressInfo::Δ, "Previous trust radius :math:`\\Delta`")
        .def_readonly("ρ", &PANTRProgressInfo::ρ, "Previous decrease ratio :math:`\\rho`")
        .def_readonly("ε", &PANTRProgressInfo::ε, "Tolerance reached :math:`\\varepsilon_k`")
        .def_readonly("Σ", &PANTRProgressInfo::Σ, "Penalty factor :math:`\\Sigma`")
        .def_readonly("y", &PANTRProgressInfo::y, "Lagrange multipliers :math:`y`")
        .def_property_readonly("problem", member_ptr<&PANTRProgressInfo::problem>(), "Problem being solved")
        .def_property_readonly("params", member_ptr<&PANTRProgressInfo::params>(), "Solver parameters")
        // clang-format on
        .def_property_readonly(
            "fpr", [](const PANTRProgressInfo &p) { return std::sqrt(p.norm_sq_p) / p.γ; },
            "Fixed-point residual :math:`\\left\\|p\\right\\| / \\gamma`");

    // Solve without ALM
    using PANTRSolver = alpaqa::PANTRSolver<TypeErasedTRDirection>;
    using Problem     = typename PANTRSolver::Problem;

    py::class_<PANTRSolver> pantr_solver(m, "PANTRSolver",
                                         "C++ documentation: :cpp:class:`alpaqa::PANTRSolver`");
    default_copy_methods(pantr_solver);
    pantr_solver
        // Constructors
        .def(py::init([](params_or_dict<PANTRParams> params,
                         params_or_dict<SteihaugCGParams> steihaug_params,
                         params_or_dict<NewtonTRDirParams> direction_params) {
                 return PANTRSolver{var_kwargs_to_struct(params),
                                    alpaqa::erase_tr_direction_with_params_dict<NewtonTRDir>(
                                        var_kwargs_to_struct(steihaug_params),
                                        var_kwargs_to_struct(direction_params))};
             }),
             "pantr_params"_a = py::dict{}, "steihaug_params"_a = py::dict{},
             "direction_params"_a = py::dict{},
             "Create a PANTR solver using a structured Newton CG subproblem solver.")
        .def(py::init(
                 [](params_or_dict<PANTRParams> params, const TypeErasedTRDirection &direction) {
                     return PANTRSolver{var_kwargs_to_struct(params),
                                        typename PANTRSolver::Direction{direction}};
                 }),
             "pantr_params"_a, "direction"_a, "Create a PANTR solver using a custom direction.")
        .def_property_readonly("direction", member_ref<&PANTRSolver::direction>());
    register_inner_solver_methods<PANTRSolver, Problem, InnerSolver>(pantr_solver);
}

template void register_pantr<alpaqa::EigenConfigd>(py::module_ &);
template void register_pantr<alpaqa::EigenConfigf>(py::module_ &);
template void register_pantr<alpaqa::EigenConfigl>(py::module_ &);
#ifdef ALPAQA_WITH_QUAD_PRECISION
template void register_pantr<alpaqa::EigenConfigq>(py::module_ &);
#endif

// TODO: Why is this needed?
template class alpaqa::PANTRSolver<alpaqa::TypeErasedTRDirection<alpaqa::EigenConfigf>>;
template class alpaqa::PANTRSolver<alpaqa::TypeErasedTRDirection<alpaqa::EigenConfigd>>;
template class alpaqa::PANTRSolver<alpaqa::TypeErasedTRDirection<alpaqa::EigenConfigl>>;
#ifdef ALPAQA_WITH_QUAD_PRECISION
template class alpaqa::PANTRSolver<alpaqa::TypeErasedTRDirection<alpaqa::EigenConfigq>>;
#endif