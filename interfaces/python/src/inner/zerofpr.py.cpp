#include <alpaqa/util/quadmath/quadmath.hpp>

#include <pybind11/chrono.h>
#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
using namespace py::literals;

#include <alpaqa/implementation/inner/zerofpr.tpp>
#include <alpaqa/inner/directions/panoc/structured-lbfgs.hpp>
#include <alpaqa/util/check-dim.hpp>

#include "type-erased-panoc-direction.hpp"
#include <dict/stats-to-dict.hpp>
#include <inner/inner-solve.hpp>
#include <inner/type-erased-inner-solver.hpp>
#include <params/params.hpp>
#include <util/async.hpp>
#include <util/copy.hpp>
#include <util/member.hpp>

template <alpaqa::Config Conf>
void register_zerofpr(py::module_ &m) {
    USING_ALPAQA_CONFIG(Conf);

    // ----------------------------------------------------------------------------------------- //
    using TypeErasedPANOCDirection = alpaqa::TypeErasedPANOCDirection<Conf>;
    using LBFGSParams              = alpaqa::LBFGSParams<config_t>;
    using StructuredLBFGSDir       = alpaqa::StructuredLBFGSDirection<config_t>;
    using StrucLBFGSDirParams      = alpaqa::StructuredLBFGSDirectionParams<config_t>;

    // ----------------------------------------------------------------------------------------- //
    using ZeroFPRParams = alpaqa::ZeroFPRParams<config_t>;
    register_dataclass<ZeroFPRParams>(m, "ZeroFPRParams",
                                      "C++ documentation: :cpp:class:`alpaqa::ZeroFPRParams`");

    // ----------------------------------------------------------------------------------------- //
    using ZeroFPRProgressInfo = alpaqa::ZeroFPRProgressInfo<config_t>;
    py::class_<ZeroFPRProgressInfo>(m, "ZeroFPRProgressInfo",
                                    "Data passed to the ZeroFPR progress callback.\n\n"
                                    "C++ documentation: :cpp:class:`alpaqa::ZeroFPRProgressInfo`")
        // clang-format off
        .def_readonly("k", &ZeroFPRProgressInfo::k, "Iteration")
        .def_readonly("status", &ZeroFPRProgressInfo::status, "Current solver status")
        .def_readonly("x", &ZeroFPRProgressInfo::x, "Decision variable :math:`x`")
        .def_readonly("p", &ZeroFPRProgressInfo::p, "Projected gradient step :math:`p`")
        .def_readonly("norm_sq_p", &ZeroFPRProgressInfo::norm_sq_p, ":math:`\\left\\|p\\right\\|^2`")
        .def_readonly("x̂", &ZeroFPRProgressInfo::x̂, "Decision variable after projected gradient step :math:`\\hat x`")
        .def_readonly("φγ", &ZeroFPRProgressInfo::φγ, "Forward-backward envelope :math:`\\varphi_\\gamma(x)`")
        .def_readonly("ψ", &ZeroFPRProgressInfo::ψ, "Objective value :math:`\\psi(x)`")
        .def_readonly("grad_ψ", &ZeroFPRProgressInfo::grad_ψ, "Gradient of objective :math:`\\nabla\\psi(x)`")
        .def_readonly("ψ_hat", &ZeroFPRProgressInfo::ψ_hat, "Objective at x̂ :math:`\\psi(\\hat x)`")
        .def_readonly("grad_ψ_hat", &ZeroFPRProgressInfo::grad_ψ_hat, "Gradient of objective at x̂ :math:`\\nabla\\psi(\\hat x)`")
        .def_readonly("q", &ZeroFPRProgressInfo::q, "Previous quasi-Newton step :math:`\\nabla\\psi(\\hat x)`")
        .def_readonly("L", &ZeroFPRProgressInfo::L, "Estimate of Lipschitz constant of objective :math:`L`")
        .def_readonly("γ", &ZeroFPRProgressInfo::γ, "Step size :math:`\\gamma`")
        .def_readonly("τ", &ZeroFPRProgressInfo::τ, "Previous line search parameter :math:`\\tau`")
        .def_readonly("ε", &ZeroFPRProgressInfo::ε, "Tolerance reached :math:`\\varepsilon_k`")
        .def_readonly("Σ", &ZeroFPRProgressInfo::Σ, "Penalty factor :math:`\\Sigma`")
        .def_readonly("y", &ZeroFPRProgressInfo::y, "Lagrange multipliers :math:`y`")
        .def_property_readonly("problem", member_ptr<&ZeroFPRProgressInfo::problem>(), "Problem being solved")
        .def_property_readonly("params", member_ptr<&ZeroFPRProgressInfo::params>(), "Solver parameters")
        // clang-format on
        .def_property_readonly(
            "fpr", [](const ZeroFPRProgressInfo &p) { return std::sqrt(p.norm_sq_p) / p.γ; },
            "Fixed-point residual :math:`\\left\\|p\\right\\| / \\gamma`");

    // Solve without ALM
    using ZeroFPRSolver = alpaqa::ZeroFPRSolver<TypeErasedPANOCDirection>;
    using Problem       = typename ZeroFPRSolver::Problem;
    using InnerSolver   = alpaqa::TypeErasedInnerSolver<config_t, Problem>;

    py::class_<ZeroFPRSolver> zfpr_solver(m, "ZeroFPRSolver",
                                          "C++ documentation: :cpp:class:`alpaqa::ZeroFPRSolver`");
    default_copy_methods(zfpr_solver);
    zfpr_solver
        // Constructors
        .def(py::init([](params_or_dict<ZeroFPRParams> params,
                         params_or_dict<LBFGSParams> lbfgs_params,
                         params_or_dict<StrucLBFGSDirParams> direction_params) {
                 return ZeroFPRSolver{var_kwargs_to_struct(params),
                                      alpaqa::erase_direction_with_params_dict<StructuredLBFGSDir>(
                                          var_kwargs_to_struct(lbfgs_params),
                                          var_kwargs_to_struct(direction_params))};
             }),
             "zerofpr_params"_a = py::dict{}, "lbfgs_params"_a = py::dict{},
             "direction_params"_a = py::dict{},
             "Create a ZeroFPR solver using structured L-BFGS directions.")
        .def(py::init([](params_or_dict<ZeroFPRParams> params,
                         const TypeErasedPANOCDirection &direction) {
                 return ZeroFPRSolver{var_kwargs_to_struct(params),
                                      typename ZeroFPRSolver::Direction{direction}};
             }),
             "zerofpr_params"_a, "direction"_a, "Create a ZeroFPR solver using a custom direction.")
        .def_property_readonly("direction", member_ref<&ZeroFPRSolver::direction>());
    register_inner_solver_methods<ZeroFPRSolver, Problem, InnerSolver>(zfpr_solver);
}

template void register_zerofpr<alpaqa::EigenConfigd>(py::module_ &);
template void register_zerofpr<alpaqa::EigenConfigf>(py::module_ &);
template void register_zerofpr<alpaqa::EigenConfigl>(py::module_ &);
#ifdef ALPAQA_WITH_QUAD_PRECISION
template void register_zerofpr<alpaqa::EigenConfigq>(py::module_ &);
#endif

// TODO: Why is this needed?
template class alpaqa::ZeroFPRSolver<alpaqa::TypeErasedPANOCDirection<alpaqa::EigenConfigf>>;
template class alpaqa::ZeroFPRSolver<alpaqa::TypeErasedPANOCDirection<alpaqa::EigenConfigd>>;
template class alpaqa::ZeroFPRSolver<alpaqa::TypeErasedPANOCDirection<alpaqa::EigenConfigl>>;
#ifdef ALPAQA_WITH_QUAD_PRECISION
template class alpaqa::ZeroFPRSolver<alpaqa::TypeErasedPANOCDirection<alpaqa::EigenConfigq>>;
#endif