#pragma once

#include <alpaqa/config/config.hpp>
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

#include "output.hpp"
#include "problem.hpp"

struct BenchmarkResults {
    USING_ALPAQA_CONFIG(alpaqa::DefaultConfig);
    static constexpr real_t NaN = alpaqa::NaN<config_t>;

    LoadedProblem &problem;
    alpaqa::EvalCounter evals;
    std::chrono::nanoseconds duration{};
    std::string_view status;
    bool success = false;
    std::string solver;
    real_t f = NaN, δ = NaN, ε = NaN, γ = NaN, Σ = NaN;
    real_t stationarity = NaN, constr_violation = NaN, complementarity = NaN;
    vec solution;
    vec multipliers;
    index_t outer_iter = -1, inner_iter = -1;
    using any_stat_t = std::variant<index_t, real_t, std::string, bool, vec,
                                    std::vector<real_t>>;
    std::vector<std::pair<std::string, any_stat_t>> extra;
    std::span<const std::string_view> options;
    int64_t timestamp = 0;
};

std::string random_hex_string(auto &&rng) {
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

template <class... A>
void dict_elem(std::ostream &os, std::string_view k, A &&...a) {
    os << "    " << std::quoted(k) << ": ";
    python_literal(os, std::forward<A>(a)...);
    os << ",\n";
};

void write_evaluations(std::ostream &os, const alpaqa::EvalCounter &evals) {
    os << std::quoted("evaluations") << ": {\n";
#define EVAL(name) dict_elem(os, #name, evals.name)
    EVAL(proj_diff_g);
    EVAL(proj_multipliers);
    EVAL(prox_grad_step);
    EVAL(f);
    EVAL(grad_f);
    EVAL(f_grad_f);
    EVAL(f_g);
    EVAL(f_grad_f_g);
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
    EVAL(grad_ψ_from_ŷ);
    EVAL(ψ_grad_ψ);
#undef EVAL
    os << "},\n";
}

void write_evaluations(std::ostream &os, const alpaqa::OCPEvalCounter &evals) {
    os << std::quoted("evaluations") << ": {\n";
#define EVAL(name) dict_elem(os, #name, evals.name)
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
    os << "},\n";
}

namespace detail {
template <class... Ts>
struct overloaded : Ts... {
    using Ts::operator()...;
};
template <class... Ts>
overloaded(Ts...) -> overloaded<Ts...>;
} // namespace detail

void print_results(std::ostream &os, const BenchmarkResults &results) {
    USING_ALPAQA_CONFIG(BenchmarkResults::config_t);
    using alpaqa::float_to_str;
    auto time_s = std::chrono::duration<double>(results.duration);
    os << '\n'
       << results.evals << '\n'
       << "solver:  " << results.solver << '\n'
       << "problem: " << results.problem.path.filename().c_str() << " (from "
       << results.problem.path.parent_path() << ")" << '\n'
       << "status:  " << (results.success ? "\033[0;32m" : "\033[0;31m")
       << results.status << "\033[0m" << '\n'
       << "num var: " << results.problem.problem.get_n() << '\n'
       << "num con: " << results.problem.problem.get_m() << '\n'
       << "f = " << float_to_str(results.f) << '\n'
       << "ε = " << float_to_str(results.ε) << '\n'
       << "δ = " << float_to_str(results.δ) << '\n'
       << "γ = " << float_to_str(results.γ) << '\n'
       << "Σ = " << float_to_str(results.Σ) << '\n'
       << "stationarity    = " << float_to_str(results.stationarity) << '\n'
       << "violation       = " << float_to_str(results.constr_violation) << '\n'
       << "complementarity = " << float_to_str(results.complementarity) << '\n'
       << "time: " << float_to_str(time_s.count(), 3) << " s\n"
       << "outer iter: " << std::setw(6) << results.outer_iter << '\n'
       << "iter:       " << std::setw(6) << results.inner_iter << '\n'
       << std::endl;
    for (const auto &[key, value] : results.extra) {
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

void write_results(std::ostream &os, const BenchmarkResults &results) {
    os << "from numpy import nan, inf\n"
          "import numpy as np\n"
          "__all__ = ['results']\n"
          "results = {\n";
    dict_elem(os, "opts", results.options);
    dict_elem(os, "time_utc_ms", results.timestamp);
    dict_elem(os, "runtime", results.duration.count());
    dict_elem(os, "iter", results.inner_iter);
    dict_elem(os, "outer_iter", results.outer_iter);
    dict_elem(os, "status", results.status);
    dict_elem(os, "success", results.success);
    dict_elem(os, "solver", results.solver);
    dict_elem(os, "solution", results.solution);
    dict_elem(os, "multipliers", results.multipliers);
    dict_elem(os, "f", results.f);
    dict_elem(os, "δ", results.δ);
    dict_elem(os, "ε", results.ε);
    dict_elem(os, "γ", results.γ);
    dict_elem(os, "Σ", results.Σ);
    dict_elem(os, "stationarity", results.stationarity);
    dict_elem(os, "constr_violation", results.constr_violation);
    dict_elem(os, "complementarity", results.complementarity);
    write_evaluations(os, results.evals);
    for (const auto &[key, value] : results.extra)
        std::visit([&, k{std::string_view{key}}](
                       const auto &v) { dict_elem(os, k, v); },
                   value);
    dict_elem(os, "nvar", results.problem.problem.get_n());
    dict_elem(os, "ncon", results.problem.problem.get_m());
    dict_elem(os, "path", results.problem.path);
    dict_elem(os, "abs_path", results.problem.abs_path);
    os << "}\n";
}

struct KKTError {
    USING_ALPAQA_CONFIG(alpaqa::DefaultConfig);
    real_t stationarity, constr_violation, complementarity;
};

KKTError compute_kkt_error(
    const alpaqa::TypeErasedProblem<alpaqa::DefaultConfig> &problem,
    alpaqa::crvec<alpaqa::DefaultConfig> x,
    alpaqa::crvec<alpaqa::DefaultConfig> y) {

    USING_ALPAQA_CONFIG(alpaqa::DefaultConfig);
    const auto n = x.size(), m = y.size();
    vec z(n), grad_Lx(n), work(n), g(m), e(m);
    // Gradient of Lagrangian, ∇ℒ(x,y) = ∇f(x) + ∇g(x) y
    problem.eval_grad_L(x, y, grad_Lx, work);
    // Eliminate normal cone of bound constraints, z = Π(x - ∇ℒ(x,y)) - x
    problem.eval_prox_grad_step(1, x, grad_Lx, work, z);
    // Stationarity, ‖Π(x - ∇ℒ(x,y)) - x‖
    auto stationarity = alpaqa::vec_util::norm_inf(z);
    // Constraints, g(x)
    problem.eval_g(x, g);
    // Distance to feasible set, e = g(x) - Π(g(x))
    problem.eval_proj_diff_g(g, e);
    // Constraint violation, ‖g(x) - Π(g(x))‖
    auto constr_violation = alpaqa::vec_util::norm_inf(e);
    // Complementary slackness
    real_t complementarity = std::inner_product(
        y.begin(), y.end(), e.begin(), real_t(0),
        [](real_t acc, real_t yv) { return std::fmax(acc, std::abs(yv)); },
        std::multiplies<>{});
    return {.stationarity     = stationarity,
            .constr_violation = constr_violation,
            .complementarity  = complementarity};
}
