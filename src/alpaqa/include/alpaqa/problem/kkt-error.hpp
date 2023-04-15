#pragma once

#include <alpaqa/problem/type-erased-problem.hpp>
#include <numeric>

namespace alpaqa {

template <Config Conf>
struct KKTError {
    USING_ALPAQA_CONFIG(Conf);
    static constexpr auto NaN = alpaqa::NaN<config_t>;
    real_t stationarity = NaN, constr_violation = NaN, complementarity = NaN;
};

template <Config Conf>
KKTError<Conf> compute_kkt_error(const alpaqa::TypeErasedProblem<Conf> &problem,
                                 alpaqa::crvec<Conf> x, alpaqa::crvec<Conf> y) {

    USING_ALPAQA_CONFIG(Conf);
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
        [](real_t acc, real_t ye) { return std::fmax(acc, std::abs(ye)); }, std::multiplies<>{});
    return {.stationarity     = stationarity,
            .constr_violation = constr_violation,
            .complementarity  = complementarity};
}

} // namespace alpaqa