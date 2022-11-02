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
#include <alpaqa/inner/panoc.hpp>
#include <alpaqa/inner/src/panoc.tpp>
#include <alpaqa/util/check-dim.hpp>

#include "kwargs-to-struct.hpp"
#include "lbfgs-params.hpp"
#include "stats-to-dict.hpp"
#include "type-erased-panoc-direction.hpp"

template <alpaqa::Config Conf>
struct kwargs_to_struct_table<alpaqa::PANOCParams<Conf>> {
    inline static const kwargs_to_struct_table_t<alpaqa::PANOCParams<Conf>> table{
        // clang-format off
        {"Lipschitz", &alpaqa::PANOCParams<Conf>::Lipschitz},
        {"max_iter", &alpaqa::PANOCParams<Conf>::max_iter},
        {"max_time", &alpaqa::PANOCParams<Conf>::max_time},
        {"τ_min", &alpaqa::PANOCParams<Conf>::τ_min},
        {"L_min", &alpaqa::PANOCParams<Conf>::L_min},
        {"L_max", &alpaqa::PANOCParams<Conf>::L_max},
        {"stop_crit", &alpaqa::PANOCParams<Conf>::stop_crit},
        {"max_no_progress", &alpaqa::PANOCParams<Conf>::max_no_progress},
        {"print_interval", &alpaqa::PANOCParams<Conf>::print_interval},
        {"print_precision", &alpaqa::PANOCParams<Conf>::print_precision},
        {"quadratic_upperbound_tolerance_factor", &alpaqa::PANOCParams<Conf>::quadratic_upperbound_tolerance_factor},
        {"update_lipschitz_in_linesearch", &alpaqa::PANOCParams<Conf>::update_lipschitz_in_linesearch},
        {"alternative_linesearch_cond", &alpaqa::PANOCParams<Conf>::alternative_linesearch_cond},
        // clang-format on
    };
};

template <alpaqa::Config Conf>
struct kwargs_to_struct_table<alpaqa::LipschitzEstimateParams<Conf>> {
    inline static const kwargs_to_struct_table_t<alpaqa::LipschitzEstimateParams<Conf>> table{
        {"L_0", &alpaqa::LipschitzEstimateParams<Conf>::L_0},
        {"δ", &alpaqa::LipschitzEstimateParams<Conf>::δ},
        {"ε", &alpaqa::LipschitzEstimateParams<Conf>::ε},
        {"Lγ_factor", &alpaqa::LipschitzEstimateParams<Conf>::Lγ_factor},
    };
};

template <alpaqa::Config Conf>
struct kwargs_to_struct_table<alpaqa::CBFGSParams<Conf>> {
    inline static const kwargs_to_struct_table_t<alpaqa::CBFGSParams<Conf>> table{
        {"α", &alpaqa::CBFGSParams<Conf>::α},
        {"ϵ", &alpaqa::CBFGSParams<Conf>::ϵ},
    };
};

template <alpaqa::Config Conf>
void register_panoc(py::module_ &m) {
    USING_ALPAQA_CONFIG(Conf);

    // ----------------------------------------------------------------------------------------- //
    using Box                      = alpaqa::Box<config_t>;
    using TypeErasedPANOCDirection = alpaqa::TypeErasedPANOCDirection<Conf>;
    py::class_<TypeErasedPANOCDirection>(m, "PANOCDirection")
        .def(py::init([](py::object o) {
            struct {
                using Problem = alpaqa::TypeErasedProblem<Conf>;
                void initialize(const Problem &problem, crvec y, crvec Σ, real_t γ_0, crvec x_0,
                                crvec x̂_0, crvec p_0, crvec grad_ψx_0) {
                    alpaqa::ScopedMallocAllower ma;
                    o.attr("initialize")(problem, y, Σ, γ_0, x_0, x̂_0, p_0, grad_ψx_0);
                }
                bool update(real_t γₖ, real_t γₙₑₓₜ, crvec xₖ, crvec xₙₑₓₜ, crvec pₖ, crvec pₙₑₓₜ,
                            crvec grad_ψxₖ, crvec grad_ψxₙₑₓₜ) {
                    alpaqa::ScopedMallocAllower ma;
                    return py::cast<bool>(
                        o.attr("update")(γₖ, γₙₑₓₜ, xₖ, xₙₑₓₜ, pₖ, pₙₑₓₜ, grad_ψxₖ, grad_ψxₙₑₓₜ));
                }
                bool apply(real_t γₖ, crvec xₖ, crvec x̂ₖ, crvec pₖ, crvec grad_ψxₖ, rvec qₖ) const {
                    alpaqa::ScopedMallocAllower ma;
                    return py::cast<bool>(o.attr("apply")(γₖ, xₖ, x̂ₖ, pₖ, grad_ψxₖ, qₖ));
                }
                void changed_γ(real_t γₖ, real_t old_γₖ) {
                    alpaqa::ScopedMallocAllower ma;
                    o.attr("changed_γ")(γₖ, old_γₖ);
                }
                void reset() {
                    alpaqa::ScopedMallocAllower ma;
                    o.attr("reset")();
                }
                std::string get_name() const { return py::cast<std::string>(py::str(o)); }

                py::object o;
            } s{std::move(o)};
            return TypeErasedPANOCDirection{std::move(s)};
        }))
        .def("__str__", &TypeErasedPANOCDirection::template get_name<>);

    // ----------------------------------------------------------------------------------------- //
    using LBFGS = alpaqa::LBFGS<config_t>;
    py::class_<LBFGS> lbfgs(m, "LBFGS", "C++ documentation :cpp:class:`alpaqa::LBFGS`");
    using LBFGSParams = typename LBFGS::Params;
    py::class_<LBFGSParams> lbfgsparams(lbfgs, "Params",
                                        "C++ documentation :cpp:class:`alpaqa::LBFGSParams`");
    using CBFGS = alpaqa::CBFGSParams<config_t>;
    py::class_<CBFGS> cbfgs(lbfgsparams, "CBFGS",
                            "C++ documentation :cpp:class:`alpaqa::CBFGSParams`");
    using LBFGSSign = typename LBFGS::Sign;
    py::enum_<LBFGSSign> lbfgssign(lbfgs, "Sign",
                                   "C++ documentation :cpp:enum:`alpaqa::LBFGS::Sign`");
    cbfgs //
        .def(py::init())
        .def(py::init(&kwargs_to_struct<CBFGS>))
        .def("to_dict", &struct_to_dict<CBFGS>)
        .def_readwrite("α", &CBFGS::α)
        .def_readwrite("ϵ", &CBFGS::ϵ);
    lbfgsparams //
        .def(py::init())
        .def(py::init(&kwargs_to_struct<LBFGSParams>))
        .def("to_dict", &struct_to_dict<LBFGSParams>)
        .def_readwrite("memory", &LBFGSParams::memory)
        .def_readwrite("min_div_fac", &LBFGSParams::min_div_fac)
        .def_readwrite("min_abs_s", &LBFGSParams::min_abs_s)
        .def_readwrite("cbfgs", &LBFGSParams::cbfgs)
        .def_readwrite("force_pos_def", &LBFGSParams::force_pos_def)
        .def_readwrite("stepsize", &LBFGSParams::stepsize);
    lbfgssign //
        .value("Positive", LBFGSSign::Positive)
        .value("Negative", LBFGSSign::Negative)
        .export_values();

    auto safe_lbfgs_update = [](LBFGS &self, crvec xk, crvec xkp1, crvec pk, crvec pkp1,
                                LBFGSSign sign, bool forced) {
        check_dim("xk", xk, self.n());
        check_dim("xkp1", xkp1, self.n());
        check_dim("pk", pk, self.n());
        check_dim("pkp1", pkp1, self.n());
        return self.update(xk, xkp1, pk, pkp1, sign, forced);
    };
    auto safe_lbfgs_update_sy = [](LBFGS &self, crvec sk, crvec yk, real_t pkp1Tpkp1, bool forced) {
        check_dim("sk", sk, self.n());
        check_dim("yk", yk, self.n());
        return self.update_sy(sk, yk, pkp1Tpkp1, forced);
    };
    auto safe_lbfgs_apply = [](LBFGS &self, rvec q, real_t γ) {
        check_dim("q", q, self.n());
        return self.apply(q, γ);
    };

    lbfgs //
        .def(py::init([](params_or_dict<LBFGSParams> params) {
                 return LBFGS{var_kwargs_to_struct(params)};
             }),
             "params"_a)
        .def(py::init([](params_or_dict<LBFGSParams> params, length_t n) {
                 return LBFGS{var_kwargs_to_struct(params), n};
             }),
             "params"_a, "n"_a)
        .def_static("update_valid", LBFGS::update_valid, "params"_a, "yᵀs"_a, "sᵀs"_a, "pᵀp"_a)
        .def("update", safe_lbfgs_update, "xk"_a, "xkp1"_a, "pk"_a, "pkp1"_a,
             "sign"_a = LBFGS::Sign::Positive, "forced"_a = false)
        .def("update_sy", safe_lbfgs_update_sy, "sk"_a, "yk"_a, "pkp1Tpkp1"_a, "forced"_a = false)
        .def("apply", safe_lbfgs_apply, "q"_a, "γ"_a)
        .def("apply_masked",
             py::overload_cast<rvec, real_t, const std::vector<index_t> &>(&LBFGS::apply_masked),
             // [](LBFGS &self, rvec q, real_t γ, const std::vector<index_t> &J) {
             //     return self.apply_masked(q, γ, J);
             // },
             "q"_a, "γ"_a, "J"_a)
        .def("reset", &LBFGS::reset)
        .def("current_history", &LBFGS::current_history)
        .def("resize", &LBFGS::resize, "n"_a)
        .def("scale_y", &LBFGS::scale_y, "factor"_a)
        .def_property_readonly("n", &LBFGS::n)
        .def(
            "s", [](LBFGS &self, index_t i) -> rvec { return self.s(i); }, ret_ref_internal, "i"_a)
        .def(
            "y", [](LBFGS &self, index_t i) -> rvec { return self.y(i); }, ret_ref_internal, "i"_a)
        .def(
            "ρ", [](LBFGS &self, index_t i) -> real_t & { return self.ρ(i); }, ret_ref_internal,
            "i"_a)
        .def(
            "α", [](LBFGS &self, index_t i) -> real_t & { return self.α(i); }, ret_ref_internal,
            "i"_a)
        .def_property_readonly("params", &LBFGS::get_params)
        .def("__str__", &LBFGS::get_name);

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
        .def_readwrite("update_lipschitz_in_linesearch", &PANOCParams::update_lipschitz_in_linesearch)
        .def_readwrite("alternative_linesearch_cond", &PANOCParams::alternative_linesearch_cond)
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
        .def_readonly("L", &PANOCProgressInfo::L, "Estimate of Lipschitz constant of objective :math:`L`")
        .def_readonly("γ", &PANOCProgressInfo::γ, "Step size :math:`\\gamma`")
        .def_readonly("τ", &PANOCProgressInfo::τ, "Line search parameter :math:`\\tau`")
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
        check_dim("Σ", Σ, m);
        if (x)
            check_dim("x", *x, n);
        else
            x = vec::Zero(n);
        if (y)
            check_dim("y", *y, m);
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
            check_dim("x", *x, n);
        else
            x = vec::Zero(n);
        vec Σ(0), y(0), err_z(0);
        auto invoke_solver = [&] {
            return solver(problem, Σ, ε, always_overwrite_results, *x, y, err_z);
        };
        if (!async) {
            auto stats = invoke_solver();
            return std::make_tuple(std::move(*x), alpaqa::conv::stats_to_dict(stats));
        } else {
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
            return std::make_tuple(std::move(*x), alpaqa::conv::stats_to_dict(stats.get()));
        }
    };

    py::class_<PANOCSolver>(m, "PANOCSolver", "C++ documentation: :cpp:class:`alpaqa::PANOCSolver`")
        .def(py::init([](params_or_dict<PANOCParams> params, const LBFGS &lbfgs) {
                 return PANOCSolver{var_kwargs_to_struct(params),
                                    alpaqa::erase_direction<LBFGS>(lbfgs)};
             }),
             "panoc_params"_a, "LBFGS"_a, "Create a PANOC solver using L-BFGS directions.")
        .def(py::init(
                 [](params_or_dict<PANOCParams> params, params_or_dict<LBFGSParams> lbfgs_params) {
                     return PANOCSolver{
                         var_kwargs_to_struct(params),
                         alpaqa::erase_direction<LBFGS>(LBFGS{var_kwargs_to_struct(lbfgs_params)})};
                 }),
             "panoc_params"_a = py::dict{}, "lbfgs_params"_a = py::dict{},
             "Create a PANOC solver using L-BFGS directions.")
        .def(py::init(
                 [](params_or_dict<PANOCParams> params, const TypeErasedPANOCDirection &direction) {
                     return PANOCSolver{var_kwargs_to_struct(params),
                                        typename PANOCSolver::Direction{direction}};
                 }),
             "panoc_params"_a, "direction"_a, "Create a PANOC solver using a custom direction.")
        .def("__call__", panoc_independent_solve,
             py::call_guard<py::scoped_ostream_redirect, py::scoped_estream_redirect>(),
             "problem"_a, "Σ"_a, "ε"_a, "x"_a = py::none(),
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
        .def("__call__", panoc_independent_solve_unconstr,
             py::call_guard<py::scoped_ostream_redirect, py::scoped_estream_redirect>(),
             "problem"_a, "ε"_a, "x"_a = py::none(), "async_"_a = false, //
             "Solve.\n\n"
             ":param problem: Problem to solve\n"
             ":param ε: Desired tolerance\n"
             ":param x: Initial guess\n"
             ":param async_: Release the GIL and run the solver on a separate thread\n"
             ":return: * Solution :math:`x`\n"
             "         * Statistics\n\n")
        .def("__str__", &PANOCSolver::get_name)
        .def("set_progress_callback", &PANOCSolver::set_progress_callback, "callback"_a,
             "Specify a callable that is invoked with some intermediate results on each iteration "
             "of the algorithm.");
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