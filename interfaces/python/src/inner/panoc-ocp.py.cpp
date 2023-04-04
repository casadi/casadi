#include <alpaqa/util/quadmath/quadmath.hpp>

#include <pybind11/chrono.h>
#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
using namespace py::literals;

#include <alpaqa/inner/panoc-ocp.hpp>
#include <alpaqa/util/check-dim.hpp>

#include <dict/stats-to-dict.hpp>
#include <inner/inner-solve.hpp>
#include <inner/type-erased-inner-solver.hpp>
#include <params/params.hpp>
#include <util/async.hpp>
#include <util/copy.hpp>
#include <util/member.hpp>

template <alpaqa::Config Conf>
void register_panoc_ocp(py::module_ &m) {
    USING_ALPAQA_CONFIG(Conf);

    using TEOCProblem    = alpaqa::TypeErasedControlProblem<config_t>;
    using InnerOCPSolver = alpaqa::TypeErasedInnerSolver<config_t, TEOCProblem>;

    // ----------------------------------------------------------------------------------------- //
    using PANOCOCPParams = alpaqa::PANOCOCPParams<config_t>;
    register_dataclass<PANOCOCPParams>(m, "PANOCOCPParams",
                                       "C++ documentation: :cpp:class:`alpaqa::PANOCOCPParams`");

    using PANOCOCPProgressInfo = alpaqa::PANOCOCPProgressInfo<config_t>;
    py::class_<PANOCOCPProgressInfo>(m, "PANOCOCPProgressInfo",
                                     "Data passed to the PANOC progress callback.\n\n"
                                     "C++ documentation: :cpp:class:`alpaqa::PANOCOCPProgressInfo`")
        // clang-format off
        .def_readonly("k", &PANOCOCPProgressInfo::k, "Iteration")
        .def_readonly("status", &PANOCOCPProgressInfo::status, "Current solver status")
        .def_readonly("xu", &PANOCOCPProgressInfo::xu, "States :math:`x` and inputs :math:`u`")
        .def_readonly("p", &PANOCOCPProgressInfo::p, "Projected gradient step :math:`p`")
        .def_readonly("norm_sq_p", &PANOCOCPProgressInfo::norm_sq_p, ":math:`\\left\\|p\\right\\|^2`")
        .def_readonly("x̂u", &PANOCOCPProgressInfo::x̂u, "Variables after projected gradient step :math:`\\hat u`")
        .def_readonly("φγ", &PANOCOCPProgressInfo::φγ, "Forward-backward envelope :math:`\\varphi_\\gamma(u)`")
        .def_readonly("ψ", &PANOCOCPProgressInfo::ψ, "Objective value :math:`\\psi(u)`")
        .def_readonly("grad_ψ", &PANOCOCPProgressInfo::grad_ψ, "Gradient of objective :math:`\\nabla\\psi(u)`")
        .def_readonly("ψ_hat", &PANOCOCPProgressInfo::ψ_hat, "Objective at x̂ :math:`\\psi(\\hat u)`")
        .def_readonly("q", &PANOCOCPProgressInfo::q, "Previous accelerated step :math:`q`")
        .def_readonly("gn", &PANOCOCPProgressInfo::gn, "Was :math:`q` a Gauss-Newton or L-BFGS step?")
        .def_readonly("nJ", &PANOCOCPProgressInfo::nJ, "Number of inactive constraints :math:`\\#\\mathcal J`")
        .def_readonly("lqr_min_rcond", &PANOCOCPProgressInfo::lqr_min_rcond, "Minimum reciprocal condition number encountered in LQR factorization")
        .def_readonly("L", &PANOCOCPProgressInfo::L, "Estimate of Lipschitz constant of objective :math:`L`")
        .def_readonly("γ", &PANOCOCPProgressInfo::γ, "Step size :math:`\\gamma`")
        .def_readonly("τ", &PANOCOCPProgressInfo::τ, "Line search parameter :math:`\\tau`")
        .def_readonly("ε", &PANOCOCPProgressInfo::ε, "Tolerance reached :math:`\\varepsilon_k`")
        .def_property_readonly("problem", member_ptr<&PANOCOCPProgressInfo::problem>(), "Problem being solved")
        .def_property_readonly("params", member_ptr<&PANOCOCPProgressInfo::params>(), "Solver parameters")
        .def_property_readonly("u", &PANOCOCPProgressInfo::u, "Inputs")
        .def_property_readonly("û", &PANOCOCPProgressInfo::û, "Inputs after projected gradient step")
        .def_property_readonly("x", &PANOCOCPProgressInfo::x, "States")
        .def_property_readonly("x̂", &PANOCOCPProgressInfo::x̂, "States after projected gradient step")
        // clang-format on
        .def_property_readonly(
            "fpr", [](const PANOCOCPProgressInfo &p) { return std::sqrt(p.norm_sq_p) / p.γ; },
            "Fixed-point residual :math:`\\left\\|p\\right\\| / \\gamma`");

    using PANOCOCPSolver = alpaqa::PANOCOCPSolver<config_t>;
    using ControlProblem = typename PANOCOCPSolver::Problem;

    py::class_<PANOCOCPSolver> panoc_ocp_solver(
        m, "PANOCOCPSolver", "C++ documentation: :cpp:class:`alpaqa::PANOCOCPSolver`");
    default_copy_methods(panoc_ocp_solver);
    panoc_ocp_solver
        // Constructor
        .def(py::init([](params_or_dict<PANOCOCPParams> params) {
                 return PANOCOCPSolver{var_kwargs_to_struct(params)};
             }),
             "panoc_params"_a, "Create a PANOC solver.");
    register_inner_solver_methods<PANOCOCPSolver, ControlProblem, InnerOCPSolver>(panoc_ocp_solver);
}

template void register_panoc_ocp<alpaqa::EigenConfigd>(py::module_ &);
template void register_panoc_ocp<alpaqa::EigenConfigf>(py::module_ &);
template void register_panoc_ocp<alpaqa::EigenConfigl>(py::module_ &);
#ifdef ALPAQA_WITH_QUAD_PRECISION
template void register_panoc_ocp<alpaqa::EigenConfigq>(py::module_ &);
#endif
