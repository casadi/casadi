#include <alpaqa/util/quadmath/quadmath.hpp>

#include <pybind11/chrono.h>
#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
using namespace py::literals;

#include <alpaqa/inner/directions/panoc/lbfgs.hpp>
#include <alpaqa/inner/directions/panoc/structured-lbfgs.hpp>
#include <alpaqa/inner/panoc.hpp>
#include <alpaqa/implementation/inner/panoc.tpp>
#include <alpaqa/util/check-dim.hpp>

#include "async.hpp"
#include "params/params.hpp"
#include "stats-to-dict.hpp"
#include "type-erased-panoc-direction.hpp"

template <alpaqa::Config Conf>
void register_panoc(py::module_ &m) {
    USING_ALPAQA_CONFIG(Conf);

    // ----------------------------------------------------------------------------------------- //
    using TypeErasedPANOCDirection = alpaqa::TypeErasedPANOCDirection<Conf>;
    using LBFGSParams              = alpaqa::LBFGSParams<config_t>;
    using StructuredLBFGSDir       = alpaqa::StructuredLBFGSDirection<config_t>;
    using StrucLBFGSDirParams      = alpaqa::StructuredLBFGSDirectionParams<config_t>;

    // ----------------------------------------------------------------------------------------- //
    using LipschitzEstimateParams = alpaqa::LipschitzEstimateParams<config_t>;
    py::class_<LipschitzEstimateParams>(
        m, "LipschitzEstimateParams",
        "C++ documentation: :cpp:class:`alpaqa::LipschitzEstimateParams`")
        .def(py::init(&dict_to_struct<LipschitzEstimateParams>))
        .def(py::init(&kwargs_to_struct<LipschitzEstimateParams>))
        .def("to_dict", &struct_to_dict<LipschitzEstimateParams>)
        .def_readwrite("L_0", &LipschitzEstimateParams::L_0)
        .def_readwrite("ε", &LipschitzEstimateParams::ε)
        .def_readwrite("δ", &LipschitzEstimateParams::δ)
        .def_readwrite("Lγ_factor", &LipschitzEstimateParams::Lγ_factor);

    using PANOCParams = alpaqa::PANOCParams<config_t>;
    py::class_<PANOCParams>(m, "PANOCParams", "C++ documentation: :cpp:class:`alpaqa::PANOCParams`")
        .def(py::init(&dict_to_struct<PANOCParams>))
        .def(py::init(&kwargs_to_struct<PANOCParams>))
        .def("to_dict", &struct_to_dict<PANOCParams>)
        // clang-format off
        .def_readwrite("Lipschitz", &PANOCParams::Lipschitz)
        .def_readwrite("max_iter", &PANOCParams::max_iter)
        .def_readwrite("max_time", &PANOCParams::max_time)
        .def_readwrite("τ_min", &PANOCParams::τ_min)
        .def_readwrite("β", &PANOCParams::β)
        .def_readwrite("L_min", &PANOCParams::L_min)
        .def_readwrite("L_max", &PANOCParams::L_max)
        .def_readwrite("stop_crit", &PANOCParams::stop_crit)
        .def_readwrite("max_no_progress", &PANOCParams::max_no_progress)
        .def_readwrite("print_interval", &PANOCParams::print_interval)
        .def_readwrite("print_precision", &PANOCParams::print_precision)
        .def_readwrite("quadratic_upperbound_tolerance_factor", &PANOCParams::quadratic_upperbound_tolerance_factor)
        .def_readwrite("linesearch_tolerance_factor", &PANOCParams::linesearch_tolerance_factor)
        // clang-format on
        ;

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
        .def_property_readonly("problem", [](const PANOCProgressInfo &p) -> auto & { return p.problem; }, "Problem being solved")
        .def_property_readonly("params", [](const PANOCProgressInfo &p) -> auto & { return p.params; }, "Solver parameters")
        // clang-format on
        .def_property_readonly(
            "fpr", [](const PANOCProgressInfo &p) { return std::sqrt(p.norm_sq_p) / p.γ; },
            "Fixed-point residual :math:`\\left\\|p\\right\\| / \\gamma`");

    // Solve without ALM
    using PANOCSolver            = alpaqa::PANOCSolver<TypeErasedPANOCDirection>;
    using Problem                = typename PANOCSolver::Problem;
    auto panoc_independent_solve = [](PANOCSolver &solver, const Problem &problem, crvec Σ,
                                      real_t ε, std::optional<vec> x, std::optional<vec> y,
                                      bool async) {
        bool always_overwrite_results = true;
        auto n                        = problem.get_n();
        auto m                        = problem.get_m();
        alpaqa::util::check_dim<config_t>("Σ", Σ, m);
        if (x)
            alpaqa::util::check_dim<config_t>("x", *x, n);
        else
            x = vec::Zero(n);
        if (y)
            alpaqa::util::check_dim<config_t>("y", *y, m);
        else
            y = vec::Zero(m);
        vec err_z          = vec::Zero(m);
        auto invoke_solver = [&] {
            return solver(problem, Σ, ε, always_overwrite_results, *x, *y, err_z);
        };
        auto &&stats = async_solve(async, solver, invoke_solver, problem);
        return std::make_tuple(std::move(*x), std::move(*y), std::move(err_z),
                               alpaqa::conv::stats_to_dict(stats));
    };
    // Solve without ALM and without constraints
    auto panoc_independent_solve_unconstr = [](PANOCSolver &solver, const Problem &problem,
                                               real_t ε, std::optional<vec> x, bool async) {
        bool always_overwrite_results = true;
        auto n                        = problem.get_n();
        auto m                        = problem.get_m();
        if (x)
            alpaqa::util::check_dim<config_t>("x", *x, n);
        else
            x = vec::Zero(n);
        vec Σ(0), y(0), err_z(0);
        auto invoke_solver = [&] {
            return solver(problem, Σ, ε, always_overwrite_results, *x, y, err_z);
        };
        auto &&stats = async_solve(async, solver, invoke_solver, problem);
        return std::make_tuple(std::move(*x), alpaqa::conv::stats_to_dict(stats));
    };

    py::class_<PANOCSolver>(m, "PANOCSolver", "C++ documentation: :cpp:class:`alpaqa::PANOCSolver`")
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
        // Copy
        .def("__copy__", [](const PANOCSolver &self) { return PANOCSolver{self}; })
        .def(
            "__deepcopy__", [](const PANOCSolver &self, py::dict) { return PANOCSolver{self}; },
            "memo"_a)
        // Call
        .def("__call__", panoc_independent_solve, "problem"_a, "Σ"_a, "ε"_a, "x"_a = py::none(),
             "y"_a = py::none(), "asynchronous"_a = true, //
             "Solve.\n\n"
             ":param problem: Problem to solve\n"
             ":param Σ: Penalty factor\n"
             ":param ε: Desired tolerance\n"
             ":param x: Initial guess\n"
             ":param y: Initial Lagrange multipliers\n\n"
             ":param asynchronous: Release the GIL and run the solver on a separate thread\n"
             ":return: * Solution :math:`x`\n"
             "         * Updated Lagrange multipliers :math:`y`\n"
             "         * Slack variable error :math:`g(x) - z`\n"
             "         * Statistics\n\n")
        .def("__call__", panoc_independent_solve_unconstr, "problem"_a, "ε"_a, "x"_a = py::none(),
             "asynchronous"_a = true, //
             "Solve.\n\n"
             ":param problem: Problem to solve\n"
             ":param ε: Desired tolerance\n"
             ":param x: Initial guess\n"
             ":param asynchronous: Release the GIL and run the solver on a separate thread\n"
             ":return: * Solution :math:`x`\n"
             "         * Statistics\n\n")
        .def("__str__", &PANOCSolver::get_name)
        .def("set_progress_callback", &PANOCSolver::set_progress_callback, "callback"_a,
             "Specify a callable that is invoked with some intermediate results on each iteration "
             "of the algorithm.")
        .def_property_readonly("direction",
                               py::cpp_function(
                                   [](const PANOCSolver &s) -> const TypeErasedPANOCDirection & {
                                       return s.direction;
                                   },
                                   py::return_value_policy::reference_internal));
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