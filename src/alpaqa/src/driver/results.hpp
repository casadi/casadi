#pragma once

#include <alpaqa/config/config.hpp>
#include <alpaqa/problem/kkt-error.hpp>
#include <alpaqa/problem/ocproblem-counters.hpp>
#include <alpaqa/problem/problem-counters.hpp>
#include <alpaqa/util/print.hpp>

#include <bit>
#include <charconv>
#include <chrono>
#include <iomanip>
#include <map>
#include <numeric>
#include <string_view>
#include <variant>

#include "problem.hpp"

struct SolverResults {
    USING_ALPAQA_CONFIG(alpaqa::DefaultConfig);
    static constexpr real_t NaN = alpaqa::NaN<config_t>;

    std::string_view status;
    bool success = false;
    alpaqa::EvalCounter evals;
    std::chrono::nanoseconds duration{};
    std::string solver;
    real_t h = NaN, δ = NaN, ε = NaN, γ = NaN, Σ = NaN;
    vec solution{};
    vec multipliers{};
    vec multipliers_bounds{};
    index_t outer_iter = -1, inner_iter = -1;
    using any_stat_t = std::variant<index_t, real_t, std::string, bool, vec,
                                    std::vector<real_t>>;
    std::vector<std::pair<std::string, any_stat_t>> extra{};
};

struct BenchmarkResults {
    USING_ALPAQA_CONFIG(alpaqa::DefaultConfig);
    static constexpr real_t NaN = alpaqa::NaN<config_t>;

    LoadedProblem &problem;
    SolverResults solver_results;
    real_t objective = NaN, smooth_objective = NaN;
    alpaqa::KKTError<config_t> error{};
    std::span<const std::string_view> options;
    int64_t timestamp = 0;
};

inline std::string random_hex_string(auto &&rng) {
    auto rnd     = std::uniform_int_distribution<uint32_t>()(rng);
    auto rnd_str = std::string(8, '0');
    std::to_chars(rnd_str.data() + std::countl_zero(rnd) / 4,
                  rnd_str.data() + rnd_str.size(), rnd, 16);
    return rnd_str;
}

template <class Clk>
auto timestamp_ms() {
    using ms    = std::chrono::milliseconds;
    auto now    = Clk::now();
    auto now_ms = std::chrono::duration_cast<ms>(now.time_since_epoch());
    return now_ms;
}

inline void write_evaluations(std::ostream &os,
                              const alpaqa::EvalCounter &evals) {
    auto dict_elem = [&os](std::string_view name, const auto &value) {
        os << name << ": " << value << '\n';
    };
#define EVAL(name) dict_elem(#name, evals.name)
    EVAL(proj_diff_g);
    EVAL(proj_multipliers);
    EVAL(prox_grad_step);
    EVAL(f);
    EVAL(grad_f);
    EVAL(f_grad_f);
    EVAL(f_g);
    EVAL(grad_f_grad_g_prod);
    EVAL(g);
    EVAL(grad_g_prod);
    EVAL(grad_gi);
    EVAL(grad_L);
    EVAL(hess_L_prod);
    EVAL(hess_L);
    EVAL(hess_ψ_prod);
    EVAL(hess_ψ);
    EVAL(ψ);
    EVAL(grad_ψ);
    EVAL(ψ_grad_ψ);
#undef EVAL
}

inline void write_evaluations(std::ostream &os,
                              const alpaqa::OCPEvalCounter &evals) {
    auto dict_elem = [&os](std::string_view name, const auto &value) {
        os << name << ": " << value << '\n';
    };
#define EVAL(name) dict_elem(#name, evals.name)
    EVAL(f);
    EVAL(jac_f);
    EVAL(grad_f_prod);
    EVAL(h);
    EVAL(h_N);
    EVAL(l);
    EVAL(l_N);
    EVAL(qr);
    EVAL(q_N);
    EVAL(add_Q);
    EVAL(add_Q_N);
    EVAL(add_R_masked);
    EVAL(add_S_masked);
    EVAL(add_R_prod_masked);
    EVAL(add_S_prod_masked);
    EVAL(constr);
    EVAL(constr_N);
    EVAL(grad_constr_prod);
    EVAL(grad_constr_prod_N);
    EVAL(add_gn_hess_constr);
    EVAL(add_gn_hess_constr_N);
#undef EVAL
}

namespace detail {
template <class... Ts>
struct overloaded : Ts... {
    using Ts::operator()...;
};
template <class... Ts>
overloaded(Ts...) -> overloaded<Ts...>;
} // namespace detail

inline void print_results(std::ostream &os, const BenchmarkResults &results) {
    USING_ALPAQA_CONFIG(BenchmarkResults::config_t);
    using alpaqa::float_to_str;
    const auto &solstats = results.solver_results;
    const auto &kkterr   = results.error;
    auto time_s          = std::chrono::duration<double>(solstats.duration);
    os << '\n'
       << solstats.evals << '\n'
       << "solver:  " << solstats.solver << '\n'
       << "problem: " << results.problem.path.filename().c_str() << " (from "
       << results.problem.path.parent_path() << ")" << '\n'
       << "status:  " << (solstats.success ? "\033[0;32m" : "\033[0;31m")
       << solstats.status << "\033[0m" << '\n'
       << "num var: " << results.problem.problem.get_n() << '\n'
       << "num con: " << results.problem.problem.get_m() << '\n'
       << "ε = " << float_to_str(solstats.ε) << '\n'
       << "δ = " << float_to_str(solstats.δ) << '\n'
       << "final step size  = " << float_to_str(solstats.γ) << '\n'
       << "penalty norm     = " << float_to_str(solstats.Σ) << '\n'
       << "nonsmooth objective = " << float_to_str(solstats.h) << '\n'
       << "smooth objective    = " << float_to_str(results.objective) << '\n'
       << "objective           = " << float_to_str(results.objective) << '\n'
       << "stationarity    = " << float_to_str(kkterr.stationarity) << '\n'
       << "violation       = " << float_to_str(kkterr.constr_violation) << '\n'
       << "complementarity = " << float_to_str(kkterr.complementarity) << '\n'
       << "time: " << float_to_str(time_s.count(), 3) << " s\n"
       << "outer iter: " << std::setw(6) << solstats.outer_iter << '\n'
       << "iter:       " << std::setw(6) << solstats.inner_iter << '\n'
       << std::endl;
    for (const auto &[key, value] : solstats.extra) {
        auto print_key = [&os, k{key}](const auto &v) {
            os << k << ": " << v << '\n';
        };
        auto print = detail::overloaded{
            [&print_key](const auto &v) { print_key(v); },
            [&print_key](real_t v) { print_key(float_to_str(v)); },
            [&print_key](bool v) { print_key(v ? "yes" : "no"); },
            [](const std::vector<real_t> &) {},
            [](const vec &) {},
        };
        std::visit(print, value);
    }
    os << std::endl;
}

inline void write_results(std::ostream &os, const BenchmarkResults &results) {
    // TODO
    (void)os;
    (void)results;
}
