#include <alpaqa/util/quadmath/quadmath.hpp>

#include <pybind11/chrono.h>
#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
using namespace py::literals;
constexpr auto ret_ref_internal = py::return_value_policy::reference_internal;

#include <alpaqa/inner/panoc-ocp.hpp>
#include <alpaqa/util/check-dim.hpp>

#include "async.hpp"
#include "params/params.hpp"
#include "stats-to-dict.hpp"

template <alpaqa::Config Conf>
void register_panoc_ocp(py::module_ &m) {
    USING_ALPAQA_CONFIG(Conf);

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
        .def_property_readonly("problem", [](const PANOCOCPProgressInfo &p) -> auto & { return p.problem; }, "Problem being solved")
        .def_property_readonly("params", [](const PANOCOCPProgressInfo &p) -> auto & { return p.params; }, "Solver parameters")
        .def_property_readonly("u", &PANOCOCPProgressInfo::u, "Inputs")
        .def_property_readonly("û", &PANOCOCPProgressInfo::û, "Inputs after projected gradient step")
        .def_property_readonly("x", &PANOCOCPProgressInfo::x, "States")
        .def_property_readonly("x̂", &PANOCOCPProgressInfo::x̂, "States after projected gradient step")
        // clang-format on
        .def_property_readonly(
            "fpr", [](const PANOCOCPProgressInfo &p) { return std::sqrt(p.norm_sq_p) / p.γ; },
            "Fixed-point residual :math:`\\left\\|p\\right\\| / \\gamma`");

    using PANOCOCPSolver             = alpaqa::PANOCOCPSolver<config_t>;
    using ControlProblem             = typename PANOCOCPSolver::Problem;
    using InnerSolveOptions          = alpaqa::InnerSolveOptions<config_t>;
    auto panoc_ocp_independent_solve = [](PANOCOCPSolver &solver, const ControlProblem &problem,
                                          const InnerSolveOptions &opts, std::optional<vec> u,
                                          std::optional<vec> y, std::optional<vec> μ, bool async) {
        auto N     = problem.get_N();
        auto nu    = problem.get_nu();
        auto nc    = problem.get_nc();
        auto nc_N  = problem.get_nc_N();
        auto m     = nc * N + nc_N;
        bool ret_y = y.has_value();
        if (u)
            alpaqa::util::check_dim<config_t>("u", *u, nu * N);
        else
            u = vec::Zero(nu * N);
        if (y)
            alpaqa::util::check_dim<config_t>("y", *y, m);
        else if (m == 0)
            y = vec{};
        else
            throw std::invalid_argument("Missing argument y");
        if (μ)
            alpaqa::util::check_dim<config_t>("μ", *μ, m);
        else if (m == 0)
            μ = vec{};
        else
            throw std::invalid_argument("Missing argument μ");
        vec err_z          = vec::Zero(m);
        auto invoke_solver = [&] { return solver(problem, opts, *u, *y, *μ, err_z); };
        auto &&stats       = async_solve(async, solver, invoke_solver, problem);
        if (ret_y)
            return py::make_tuple(std::move(*u), std::move(*y), std::move(err_z),
                                  alpaqa::conv::stats_to_dict(stats));
        else
            return py::make_tuple(std::move(*u), alpaqa::conv::stats_to_dict(stats));
    };

    py::class_<PANOCOCPSolver>(m, "PANOCOCPSolver",
                               "C++ documentation: :cpp:class:`alpaqa::PANOCOCPSolver`")
        // Constructor
        .def(py::init([](params_or_dict<PANOCOCPParams> params) {
                 return PANOCOCPSolver{var_kwargs_to_struct(params)};
             }),
             "panoc_params"_a, "Create a PANOC solver.")
        // Copy
        .def("__copy__", [](const PANOCOCPSolver &self) { return PANOCOCPSolver{self}; })
        .def(
            "__deepcopy__",
            [](const PANOCOCPSolver &self, py::dict) { return PANOCOCPSolver{self}; }, "memo"_a)
        // Call
        .def("__call__", panoc_ocp_independent_solve, "problem"_a, "opts"_a = py::dict(),
             "u"_a = py::none(), "y"_a = py::none(), "μ"_a = py::none(), "asynchronous"_a = true,
             "Solve.\n\n"
             ":param problem: Problem to solve\n"
             ":param opts: Options\n"
             ":param u: Initial guess\n"
             ":param y: Lagrange multipliers\n"
             ":param μ: Penalty factors\n"
             ":param asynchronous: Release the GIL and run the solver on a separate thread\n"
             ":return: * Solution :math:`u`\n"
             "         * Statistics\n\n")
        .def("__str__", &PANOCOCPSolver::get_name)
        .def("set_progress_callback", &PANOCOCPSolver::set_progress_callback, "callback"_a,
             "Specify a callable that is invoked with some intermediate results on each iteration "
             "of the algorithm.");
}

template void register_panoc_ocp<alpaqa::EigenConfigd>(py::module_ &);
template void register_panoc_ocp<alpaqa::EigenConfigf>(py::module_ &);
template void register_panoc_ocp<alpaqa::EigenConfigl>(py::module_ &);
#ifdef ALPAQA_WITH_QUAD_PRECISION
template void register_panoc_ocp<alpaqa::EigenConfigq>(py::module_ &);
#endif
