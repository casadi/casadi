#include <alpaqa/util/quadmath/quadmath.hpp>

#include <pybind11/chrono.h>
#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include <pybind11/gil.h>
#include <pybind11/iostream.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
using namespace py::literals;
constexpr auto ret_ref_internal = py::return_value_policy::reference_internal;

#include <chrono>
#include <exception>
#include <future>
using namespace std::chrono_literals;

#include <alpaqa/accelerators/lbfgs.hpp>
#include <alpaqa/inner/directions/panoc/lbfgs.hpp>
#include <alpaqa/inner/directions/panoc/structured-lbfgs.hpp>
#include <alpaqa/inner/panoc-ocp.hpp>
#include <alpaqa/inner/panoc.hpp>
#include <alpaqa/inner/src/panoc.tpp>
#include <alpaqa/util/check-dim.hpp>

#include "params/params.hpp"
#include "stats-to-dict.hpp"
#include "stream-replacer.hpp"
#include "thread-checker.hpp"
#include "type-erased-panoc-direction.hpp"

template <alpaqa::Config Conf>
void register_panoc(py::module_ &m) {
    USING_ALPAQA_CONFIG(Conf);

    // ----------------------------------------------------------------------------------------- //
    using Box                      = alpaqa::Box<config_t>;
    using TypeErasedPANOCDirection = alpaqa::TypeErasedPANOCDirection<Conf>;
    py::class_<TypeErasedPANOCDirection> te_direction(m, "PANOCDirection");
    te_direction //
        .def_property_readonly("params", &TypeErasedPANOCDirection::template get_params<>)
        .def("__str__", &TypeErasedPANOCDirection::template get_name<>);

    // ----------------------------------------------------------------------------------------- //
    using LBFGSDir       = alpaqa::LBFGSDirection<config_t>;
    using LBFGSParams    = alpaqa::LBFGSParams<config_t>;
    using LBFGSDirParams = alpaqa::LBFGSDirectionParams<config_t>;

    py::class_<LBFGSDir> lbfgs(m, "LBFGSDirection",
                               "C++ documentation: :cpp:class:`alpaqa::LBFGSDirection`");
    py::class_<LBFGSDirParams> lbfgs_params(
        lbfgs, "DirectionParams",
        "C++ documentation: :cpp:class:`alpaqa::LBFGSDirection::DirectionParams`");
    lbfgs_params //
        .def(py::init())
        .def(py::init(&kwargs_to_struct<LBFGSDirParams>))
        .def("to_dict", &struct_to_dict<LBFGSDirParams>)
        .def_readwrite("rescale_when_γ_changes", &LBFGSDirParams::rescale_when_γ_changes);
    lbfgs //
        .def(py::init([](params_or_dict<LBFGSParams> lbfgs_params,
                         params_or_dict<LBFGSDirParams> direction_params) {
                 return LBFGSDir{var_kwargs_to_struct(lbfgs_params),
                                 var_kwargs_to_struct(direction_params)};
             }),
             "lbfgs_params"_a = py::dict{}, "direction_params"_a = py::dict{})
        .def_property_readonly(
            "params",
            py::cpp_function(&LBFGSDir::get_params, py::return_value_policy::reference_internal))
        .def("__str__", &LBFGSDir::get_name);

    te_direction.def(
        py::init(&alpaqa::erase_direction_with_params_dict<LBFGSDir, const LBFGSDir &>));
    py::implicitly_convertible<LBFGSDir, TypeErasedPANOCDirection>();

    // ----------------------------------------------------------------------------------------- //
    using StructuredLBFGSDir  = alpaqa::StructuredLBFGSDirection<config_t>;
    using StrucLBFGSDirParams = alpaqa::StructuredLBFGSDirectionParams<config_t>;

    py::class_<StructuredLBFGSDir> struc_lbfgs(
        m, "StructuredLBFGSDirection",
        "C++ documentation: :cpp:class:`alpaqa::StructuredLBFGSDirection`");
    py::class_<StrucLBFGSDirParams> struc_lbfgs_params(
        struc_lbfgs, "DirectionParams",
        "C++ documentation: :cpp:class:`alpaqa::StructuredLBFGSDirection::DirectionParams`");
    struc_lbfgs_params //
        .def(py::init())
        .def(py::init(&kwargs_to_struct<StrucLBFGSDirParams>))
        .def("to_dict", &struct_to_dict<StrucLBFGSDirParams>)
        .def_readwrite("hessian_vec", &StrucLBFGSDirParams::hessian_vec)
        .def_readwrite("hessian_vec_finite_differences",
                       &StrucLBFGSDirParams::hessian_vec_finite_differences)
        .def_readwrite("full_augmented_hessian", &StrucLBFGSDirParams::full_augmented_hessian);
    struc_lbfgs //
        .def(py::init([](params_or_dict<LBFGSParams> lbfgs_params,
                         params_or_dict<StrucLBFGSDirParams> direction_params) {
                 return StructuredLBFGSDir{var_kwargs_to_struct(lbfgs_params),
                                           var_kwargs_to_struct(direction_params)};
             }),
             "lbfgs_params"_a = py::dict{}, "direction_params"_a = py::dict{})
        .def_property_readonly("params",
                               py::cpp_function(&StructuredLBFGSDir::get_params,
                                                py::return_value_policy::reference_internal))
        .def("__str__", &StructuredLBFGSDir::get_name);

    te_direction.def(py::init(
        &alpaqa::erase_direction_with_params_dict<StructuredLBFGSDir, const StructuredLBFGSDir &>));
    py::implicitly_convertible<StructuredLBFGSDir, TypeErasedPANOCDirection>();

    // ----------------------------------------------------------------------------------------- //
    using LipschitzEstimateParams = alpaqa::LipschitzEstimateParams<config_t>;
    py::class_<LipschitzEstimateParams>(
        m, "LipschitzEstimateParams",
        "C++ documentation: :cpp:class:`alpaqa::LipschitzEstimateParams`")
        .def(py::init())
        .def(py::init(&kwargs_to_struct<LipschitzEstimateParams>))
        .def("to_dict", &struct_to_dict<LipschitzEstimateParams>)
        .def_readwrite("L_0", &LipschitzEstimateParams::L_0)
        .def_readwrite("ε", &LipschitzEstimateParams::ε)
        .def_readwrite("δ", &LipschitzEstimateParams::δ)
        .def_readwrite("Lγ_factor", &LipschitzEstimateParams::Lγ_factor);

    using PANOCParams = alpaqa::PANOCParams<config_t>;
    py::class_<PANOCParams>(m, "PANOCParams", "C++ documentation: :cpp:class:`alpaqa::PANOCParams`")
        .def(py::init())
        .def(py::init(&kwargs_to_struct<PANOCParams>))
        .def("to_dict", &struct_to_dict<PANOCParams>)
        // clang-format off
        .def_readwrite("Lipschitz", &PANOCParams::Lipschitz)
        .def_readwrite("max_iter", &PANOCParams::max_iter)
        .def_readwrite("max_time", &PANOCParams::max_time)
        .def_readwrite("τ_min", &PANOCParams::τ_min)
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

    using PANOCSolver            = alpaqa::PANOCSolver<TypeErasedPANOCDirection>;
    using Problem                = typename PANOCSolver::Problem;
    auto panoc_independent_solve = [](PANOCSolver &solver, const Problem &problem, crvec Σ,
                                      real_t ε, std::optional<vec> x, std::optional<vec> y) {
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
        vec err_z  = vec::Zero(m);
        auto stats = solver(problem, Σ, ε, always_overwrite_results, *x, *y, err_z);
        return std::make_tuple(std::move(*x), std::move(*y), std::move(err_z),
                               alpaqa::conv::stats_to_dict(stats));
    };
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
        if (!async) {
            // Replace the output stream
            StreamReplacer stream{&solver};
            // Invoke the solver synchronously
            auto stats = invoke_solver();
            return std::make_tuple(std::move(*x), alpaqa::conv::stats_to_dict(stats));
        } else {
            // Check that the user doesn't use the same solver/problem in multiple threads
            ThreadChecker solver_checker{&solver};
            ThreadChecker problem_checker{&problem};
            // Replace the output stream
            StreamReplacer stream{&solver};
            // Invoke the solver asynchronously
            auto stats = std::async(std::launch::async, invoke_solver);
            {
                py::gil_scoped_release gil;
                while (stats.wait_for(50ms) != std::future_status::ready) {
                    py::gil_scoped_acquire gil;
                    // Check if Python received a signal (e.g. Ctrl+C)
                    if (PyErr_CheckSignals() != 0) {
                        // Nicely ask the solver to stop
                        solver.stop();
                        // It should return a result soon
                        if (py::gil_scoped_release gil;
                            stats.wait_for(15s) != std::future_status::ready) {
                            // If it doesn't, we terminate the entire program,
                            // because the solver uses variables local to this
                            // function, so we cannot safely return without
                            // waiting for the solver to finish.
                            std::terminate();
                        }
                        if (PyErr_Occurred())
                            throw py::error_already_set();
                        break;
                    }
                }
            }
            return std::make_tuple(std::move(*x), alpaqa::conv::stats_to_dict(stats.get()));
        }
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
             "y"_a = py::none(), //
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

    using PANOCOCPParams = alpaqa::PANOCOCPParams<config_t>;
    py::class_<PANOCOCPParams>(m, "PANOCOCPParams",
                               "C++ documentation: :cpp:class:`alpaqa::PANOCOCPParams`")
        .def(py::init())
        .def(py::init(&kwargs_to_struct<PANOCOCPParams>))
        .def("to_dict", &struct_to_dict<PANOCOCPParams>)
        // clang-format off
        .def_readwrite("Lipschitz", &PANOCOCPParams::Lipschitz)
        .def_readwrite("max_iter", &PANOCOCPParams::max_iter)
        .def_readwrite("max_time", &PANOCOCPParams::max_time)
        .def_readwrite("τ_min", &PANOCOCPParams::τ_min)
        .def_readwrite("β", &PANOCOCPParams::β)
        .def_readwrite("L_min", &PANOCOCPParams::L_min)
        .def_readwrite("L_max", &PANOCOCPParams::L_max)
        .def_readwrite("L_max_inc", &PANOCOCPParams::L_max_inc)
        .def_readwrite("stop_crit", &PANOCOCPParams::stop_crit)
        .def_readwrite("max_no_progress", &PANOCOCPParams::max_no_progress)
        .def_readwrite("gn_interval", &PANOCOCPParams::gn_interval)
        .def_readwrite("gn_sticky", &PANOCOCPParams::gn_sticky)
        .def_readwrite("reset_lbfgs_on_gn_step", &PANOCOCPParams::reset_lbfgs_on_gn_step)
        .def_readwrite("lqr_factor_cholesky", &PANOCOCPParams::lqr_factor_cholesky)
        .def_readwrite("print_interval", &PANOCOCPParams::print_interval)
        .def_readwrite("print_precision", &PANOCOCPParams::print_precision)
        .def_readwrite("quadratic_upperbound_tolerance_factor", &PANOCOCPParams::quadratic_upperbound_tolerance_factor)
        .def_readwrite("linesearch_tolerance_factor", &PANOCOCPParams::linesearch_tolerance_factor)
        .def_readwrite("disable_acceleration", &PANOCOCPParams::disable_acceleration)
        // clang-format on
        ;

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
        .def_property_readonly("qu", &PANOCOCPProgressInfo::qu, "Accelerated step on inputs")
        .def_property_readonly("qx", &PANOCOCPProgressInfo::qx, "Accelerated step on states")
        // clang-format on
        .def_property_readonly(
            "fpr", [](const PANOCOCPProgressInfo &p) { return std::sqrt(p.norm_sq_p) / p.γ; },
            "Fixed-point residual :math:`\\left\\|p\\right\\| / \\gamma`");

    using PANOCOCPSolver                      = alpaqa::PANOCOCPSolver<config_t>;
    using ControlProblem                      = typename PANOCOCPSolver::Problem;
    auto panoc_ocp_independent_solve_unconstr = [](PANOCOCPSolver &solver,
                                                   const ControlProblem &problem, real_t ε,
                                                   std::optional<vec> u, bool async) {
        auto N  = problem.get_N();
        auto nu = problem.get_nu();
        if (u)
            alpaqa::util::check_dim<config_t>("u", *u, nu * N);
        else
            u = vec::Zero(nu * N);
        auto invoke_solver = [&] { return solver(problem, ε, *u); };
        if (!async) {
            auto stats = invoke_solver();
            return std::make_tuple(std::move(*u), alpaqa::conv::stats_to_dict(stats));
        } else {
            // Check that the user doesn't use the same solver/problem in multiple threads
            ThreadChecker solver_checker{&solver};
            ThreadChecker problem_checker{&problem};
            // Replace the output stream
            StreamReplacer stream{&solver};
            // Invoke the solver asynchronously
            auto stats = std::async(std::launch::async, invoke_solver);
            {
                py::gil_scoped_release gil{};
                while (stats.wait_for(50ms) != std::future_status::ready) {
                    py::gil_scoped_acquire gil{};
                    // Check if Python received a signal (e.g. Ctrl+C)
                    if (PyErr_CheckSignals() != 0) {
                        // Nicely ask the solver to stop
                        solver.stop();
                        // It should return a result soon
                        if (py::gil_scoped_release gil{};
                            stats.wait_for(15s) != std::future_status::ready) {
                            // If it doesn't, we terminate the entire program,
                            // because the solver uses variables local to this
                            // function, so we cannot safely return without
                            // waiting for the solver to finish.
                            std::terminate();
                        }
                        if (PyErr_Occurred())
                            throw py::error_already_set();
                        break;
                    }
                }
            }
            return std::make_tuple(std::move(*u), alpaqa::conv::stats_to_dict(stats.get()));
        }
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
        .def("__call__", panoc_ocp_independent_solve_unconstr, "problem"_a, "ε"_a,
             "u"_a = py::none(), "asynchronous"_a = true, //
             "Solve.\n\n"
             ":param problem: Problem to solve\n"
             ":param ε: Desired tolerance\n"
             ":param u: Initial guess\n"
             ":param asynchronous: Release the GIL and run the solver on a separate thread\n"
             ":return: * Solution :math:`u`\n"
             "         * Statistics\n\n")
        .def("__str__", &PANOCOCPSolver::get_name)
        .def("set_progress_callback", &PANOCOCPSolver::set_progress_callback, "callback"_a,
             "Specify a callable that is invoked with some intermediate results on each iteration "
             "of the algorithm.");

    // Catch-all, must be last
    te_direction //
        .def(py::init([](py::object o) {
            struct {
                using Problem = alpaqa::TypeErasedProblem<Conf>;
                void initialize(const Problem &problem, crvec y, crvec Σ, real_t γ_0, crvec x_0,
                                crvec x̂_0, crvec p_0, crvec grad_ψx_0) {
                    alpaqa::ScopedMallocAllower ma;
                    py::gil_scoped_acquire gil;
                    o.attr("initialize")(problem, y, Σ, γ_0, x_0, x̂_0, p_0, grad_ψx_0);
                }
                bool update(real_t γₖ, real_t γₙₑₓₜ, crvec xₖ, crvec xₙₑₓₜ, crvec pₖ, crvec pₙₑₓₜ,
                            crvec grad_ψxₖ, crvec grad_ψxₙₑₓₜ) {
                    alpaqa::ScopedMallocAllower ma;
                    py::gil_scoped_acquire gil;
                    return py::cast<bool>(
                        o.attr("update")(γₖ, γₙₑₓₜ, xₖ, xₙₑₓₜ, pₖ, pₙₑₓₜ, grad_ψxₖ, grad_ψxₙₑₓₜ));
                }
                bool has_initial_direction() const {
                    alpaqa::ScopedMallocAllower ma;
                    py::gil_scoped_acquire gil;
                    return py::cast<bool>(o.attr("has_initial_direction")());
                }
                bool apply(real_t γₖ, crvec xₖ, crvec x̂ₖ, crvec pₖ, crvec grad_ψxₖ, rvec qₖ) const {
                    alpaqa::ScopedMallocAllower ma;
                    py::gil_scoped_acquire gil;
                    return py::cast<bool>(o.attr("apply")(γₖ, xₖ, x̂ₖ, pₖ, grad_ψxₖ, qₖ));
                }
                void changed_γ(real_t γₖ, real_t old_γₖ) {
                    alpaqa::ScopedMallocAllower ma;
                    py::gil_scoped_acquire gil;
                    o.attr("changed_γ")(γₖ, old_γₖ);
                }
                void reset() {
                    alpaqa::ScopedMallocAllower ma;
                    py::gil_scoped_acquire gil;
                    o.attr("reset")();
                }
                std::string get_name() const {
                    py::gil_scoped_acquire gil;
                    return py::cast<std::string>(py::str(o));
                }
                py::object get_params() const {
                    py::gil_scoped_acquire gil;
                    return py::getattr(o, "params");
                }

                py::object o;
            } s{std::move(o)};
            return TypeErasedPANOCDirection{std::move(s)};
        }));
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