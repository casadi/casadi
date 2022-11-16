#include <alpaqa/util/quadmath/quadmath.hpp>

#include <pybind11/chrono.h>
#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
using namespace py::literals;
constexpr auto ret_ref_internal = py::return_value_policy::reference_internal;

#include <alpaqa/inner/src/structured-panoc.tpp>
#include <alpaqa/inner/structured-panoc.hpp>
#include <alpaqa/util/check-dim.hpp>

#include "kwargs-to-struct.hpp"
#include "lbfgs-params.hpp"
#include "type-erased-panoc-direction.hpp"

template <alpaqa::Config Conf>
struct kwargs_to_struct_table<alpaqa::StructuredPANOCLBFGSParams<Conf>> {
    inline static const kwargs_to_struct_table_t<alpaqa::StructuredPANOCLBFGSParams<Conf>> table{
        // clang-format off
        {"Lipschitz", &alpaqa::StructuredPANOCLBFGSParams<Conf>::Lipschitz},
        {"max_iter", &alpaqa::StructuredPANOCLBFGSParams<Conf>::max_iter},
        {"max_time", &alpaqa::StructuredPANOCLBFGSParams<Conf>::max_time},
        {"τ_min", &alpaqa::StructuredPANOCLBFGSParams<Conf>::τ_min},
        {"L_min", &alpaqa::StructuredPANOCLBFGSParams<Conf>::L_min},
        {"L_max", &alpaqa::StructuredPANOCLBFGSParams<Conf>::L_max},
        {"nonmonotone_linesearch", &alpaqa::StructuredPANOCLBFGSParams<Conf>::nonmonotone_linesearch},
        {"fpr_shortcut_accept_factor", &alpaqa::StructuredPANOCLBFGSParams<Conf>::fpr_shortcut_accept_factor},
        {"fpr_shortcut_history", &alpaqa::StructuredPANOCLBFGSParams<Conf>::fpr_shortcut_history},
        {"L_max", &alpaqa::StructuredPANOCLBFGSParams<Conf>::L_max},
        {"stop_crit", &alpaqa::StructuredPANOCLBFGSParams<Conf>::stop_crit},
        {"max_no_progress", &alpaqa::StructuredPANOCLBFGSParams<Conf>::max_no_progress},
        {"print_interval", &alpaqa::StructuredPANOCLBFGSParams<Conf>::print_interval},
        {"print_precision", &alpaqa::StructuredPANOCLBFGSParams<Conf>::print_precision},
        {"quadratic_upperbound_tolerance_factor", &alpaqa::StructuredPANOCLBFGSParams<Conf>::quadratic_upperbound_tolerance_factor},
        {"linesearch_tolerance_factor", &alpaqa::StructuredPANOCLBFGSParams<Conf>::linesearch_tolerance_factor},
        {"update_lipschitz_in_linesearch", &alpaqa::StructuredPANOCLBFGSParams<Conf>::update_lipschitz_in_linesearch},
        {"alternative_linesearch_cond", &alpaqa::StructuredPANOCLBFGSParams<Conf>::alternative_linesearch_cond},
        {"hessian_vec", &alpaqa::StructuredPANOCLBFGSParams<Conf>::hessian_vec},
        {"hessian_vec_finite_differences", &alpaqa::StructuredPANOCLBFGSParams<Conf>::hessian_vec_finite_differences},
        {"full_augmented_hessian", &alpaqa::StructuredPANOCLBFGSParams<Conf>::full_augmented_hessian},
        {"hessian_step_size_heuristic", &alpaqa::StructuredPANOCLBFGSParams<Conf>::hessian_step_size_heuristic},
        // clang-format on
    };
};

template <alpaqa::Config Conf>
void register_structured_panoc(py::module_ &m) {
    USING_ALPAQA_CONFIG(Conf);

    using LBFGSParams                = alpaqa::LBFGSParams<config_t>;
    using StructuredPANOCLBFGSParams = alpaqa::StructuredPANOCLBFGSParams<config_t>;
    py::class_<StructuredPANOCLBFGSParams>(
        m, "StructuredPANOCLBFGSParams",
        "C++ documentation: :cpp:class:`alpaqa::StructuredPANOCLBFGSParams`")
        .def(py::init())
        .def(py::init(&kwargs_to_struct<StructuredPANOCLBFGSParams>))
        .def("to_dict", &struct_to_dict<StructuredPANOCLBFGSParams>)
        // clang-format off
        .def_readwrite("Lipschitz", &StructuredPANOCLBFGSParams::Lipschitz)
        .def_readwrite("max_iter", &StructuredPANOCLBFGSParams::max_iter)
        .def_readwrite("max_time", &StructuredPANOCLBFGSParams::max_time)
        .def_readwrite("τ_min", &StructuredPANOCLBFGSParams::τ_min)
        .def_readwrite("L_min", &StructuredPANOCLBFGSParams::L_min)
        .def_readwrite("L_max", &StructuredPANOCLBFGSParams::L_max)
        .def_readwrite("nonmonotone_linesearch", &StructuredPANOCLBFGSParams::nonmonotone_linesearch)
        .def_readwrite("fpr_shortcut_accept_factor", &StructuredPANOCLBFGSParams::fpr_shortcut_accept_factor)
        .def_readwrite("fpr_shortcut_history", &StructuredPANOCLBFGSParams::fpr_shortcut_history)
        .def_readwrite("L_max", &StructuredPANOCLBFGSParams::L_max)
        .def_readwrite("stop_crit", &StructuredPANOCLBFGSParams::stop_crit)
        .def_readwrite("max_no_progress", &StructuredPANOCLBFGSParams::max_no_progress)
        .def_readwrite("print_interval", &StructuredPANOCLBFGSParams::print_interval)
        .def_readwrite("print_precision", &StructuredPANOCLBFGSParams::print_precision)
        .def_readwrite("quadratic_upperbound_tolerance_factor", &StructuredPANOCLBFGSParams::quadratic_upperbound_tolerance_factor)
        .def_readwrite("linesearch_tolerance_factor", &StructuredPANOCLBFGSParams::linesearch_tolerance_factor)
        .def_readwrite("update_lipschitz_in_linesearch", &StructuredPANOCLBFGSParams::update_lipschitz_in_linesearch)
        .def_readwrite("alternative_linesearch_cond", &StructuredPANOCLBFGSParams::alternative_linesearch_cond)
        .def_readwrite("hessian_vec", &StructuredPANOCLBFGSParams::hessian_vec)
        .def_readwrite("hessian_vec_finite_differences", &StructuredPANOCLBFGSParams::hessian_vec_finite_differences)
        .def_readwrite("full_augmented_hessian", &StructuredPANOCLBFGSParams::full_augmented_hessian)
        .def_readwrite("hessian_step_size_heuristic", &StructuredPANOCLBFGSParams::hessian_step_size_heuristic)
        // clang-format on
        ;

    // ----------------------------------------------------------------------------------------- //
    using StructuredPANOCLBFGSProgressInfo = alpaqa::StructuredPANOCLBFGSProgressInfo<config_t>;
    py::class_<StructuredPANOCLBFGSProgressInfo>(
        m, "StructuredPANOCLBFGSProgressInfo",
        "Data passed to the PANOC progress callback.\n\n"
        "C++ documentation: :cpp:class:`alpaqa::StructuredPANOCLBFGSProgressInfo`")
        // clang-format off
        .def_readonly("k", &StructuredPANOCLBFGSProgressInfo::k, "Iteration")
        .def_readonly("x", &StructuredPANOCLBFGSProgressInfo::x, "Decision variable :math:`x`")
        .def_readonly("p", &StructuredPANOCLBFGSProgressInfo::p, "Projected gradient step :math:`p`")
        .def_readonly("norm_sq_p", &StructuredPANOCLBFGSProgressInfo::norm_sq_p, ":math:`\\left\\|p\\right\\|^2`")
        .def_readonly("x̂", &StructuredPANOCLBFGSProgressInfo::x̂, "Decision variable after projected gradient step :math:`\\hat x`")
        .def_readonly("q", &StructuredPANOCLBFGSProgressInfo::q, "Previous quasi-Newton step :math:`q`")
        .def_readonly("J", &StructuredPANOCLBFGSProgressInfo::J, "Index set :math:`\\mathcal J`")
        .def_readonly("φγ", &StructuredPANOCLBFGSProgressInfo::φγ, "Forward-backward envelope :math:`\\varphi_\\gamma(x)`")
        .def_readonly("ψ", &StructuredPANOCLBFGSProgressInfo::ψ, "Objective value :math:`\\psi(x)`")
        .def_readonly("grad_ψ", &StructuredPANOCLBFGSProgressInfo::grad_ψ, "Gradient of objective :math:`\\nabla\\psi(x)`")
        .def_readonly("ψ_hat", &StructuredPANOCLBFGSProgressInfo::ψ_hat)
        .def_readonly("grad_ψ_hat", &StructuredPANOCLBFGSProgressInfo::grad_ψ_hat)
        .def_readonly("L", &StructuredPANOCLBFGSProgressInfo::L, "Estimate of Lipschitz constant of objective :math:`L`")
        .def_readonly("γ", &StructuredPANOCLBFGSProgressInfo::γ, "Step size :math:`\\gamma`")
        .def_readonly("τ", &StructuredPANOCLBFGSProgressInfo::τ, "Line search parameter :math:`\\tau`")
        .def_readonly("ε", &StructuredPANOCLBFGSProgressInfo::ε, "Tolerance reached :math:`\\varepsilon_k`")
        .def_readonly("Σ", &StructuredPANOCLBFGSProgressInfo::Σ, "Penalty factor :math:`\\Sigma`")
        .def_readonly("y", &StructuredPANOCLBFGSProgressInfo::y, "Lagrange multipliers :math:`y`")
        .def_property_readonly("problem", [](const StructuredPANOCLBFGSProgressInfo &p) -> auto & { return p.problem; }, "Problem being solved")
        .def_property_readonly("params", [](const StructuredPANOCLBFGSProgressInfo &p) -> auto & { return p.params; }, "Solver parameters")
        // clang-format on
        .def_property_readonly(
            "fpr",
            [](const StructuredPANOCLBFGSProgressInfo &p) { return std::sqrt(p.norm_sq_p) / p.γ; },
            "Fixed-point residual :math:`\\left\\|p\\right\\| / \\gamma`");

    using StructuredPANOCLBFGSSolver = alpaqa::StructuredPANOCLBFGSSolver<config_t>;
    py::class_<StructuredPANOCLBFGSSolver>(
        m, "StructuredPANOCLBFGSSolver",
        "C++ documentation: :cpp:class:`alpaqa::StructuredPANOCLBFGSSolver`")
        .def(py::init([](params_or_dict<StructuredPANOCLBFGSParams> params,
                         params_or_dict<LBFGSParams> lbfgs_params) {
                 return StructuredPANOCLBFGSSolver{var_kwargs_to_struct(params),
                                                   var_kwargs_to_struct(lbfgs_params)};
             }),
             "panoc_params"_a = py::dict{}, "lbfgs_params"_a = py::dict{},
             "Create a PANOC solver using L-BFGS directions.")
        .def("set_progress_callback", &StructuredPANOCLBFGSSolver::set_progress_callback,
             "callback"_a,
             "Specify a callable that is invoked with some intermediate results on each iteration "
             "of the algorithm.");
}

template void register_structured_panoc<alpaqa::EigenConfigd>(py::module_ &);
template void register_structured_panoc<alpaqa::EigenConfigf>(py::module_ &);
template void register_structured_panoc<alpaqa::EigenConfigl>(py::module_ &);
#ifdef ALPAQA_WITH_QUAD_PRECISION
template void register_structured_panoc<alpaqa::EigenConfigq>(py::module_ &);
#endif
