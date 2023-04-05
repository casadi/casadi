#pragma once

#include <alpaqa/problem/problem-counters.hpp>
#include <alpaqa/problem/type-erased-problem.hpp>
#include <alpaqa/util/timed.hpp>

#include <functional>
#include <type_traits>

namespace alpaqa {

/// @addtogroup grp_Problems
/// @{

/// Problem wrapper that keeps track of the number of evaluations and the run
/// time of each function.
/// You probably want to use @ref problem_with_counters or
/// @ref problem_with_counters_ref instead of instantiating this class directly.
/// @note   The evaluation counters are stored using a `std::shared_pointers`,
///         which means that different copies of a @ref ProblemWithCounters
///         instance all share the same counters. To opt out of this behavior,
///         you can use the @ref decouple_evaluations function.
template <class Problem>
struct ProblemWithCounters {
    USING_ALPAQA_CONFIG_TEMPLATE(std::remove_cvref_t<Problem>::config_t);
    using Box = typename TypeErasedProblem<config_t>::Box;

    // clang-format off
    void eval_proj_diff_g(crvec z, rvec e) const { ++evaluations->proj_diff_g; return timed(evaluations->time.proj_diff_g, std::bind(&std::remove_cvref_t<Problem>::eval_proj_diff_g, &problem, z, e)); }
    void eval_proj_multipliers(rvec y, real_t M) const { ++evaluations->proj_multipliers; return timed(evaluations->time.proj_multipliers, std::bind(&std::remove_cvref_t<Problem>::eval_proj_multipliers, &problem, y, M)); }
    real_t eval_prox_grad_step(real_t γ, crvec x, crvec grad_ψ, rvec x̂, rvec p) const { ++evaluations->prox_grad_step; return timed(evaluations->time.prox_grad_step, std::bind(&std::remove_cvref_t<Problem>::eval_prox_grad_step, &problem, γ, x, grad_ψ, x̂, p)); }
    index_t eval_inactive_indices_res_lna(real_t γ, crvec x, crvec grad_ψ, rindexvec J) const requires requires { &std::remove_cvref_t<Problem>::eval_inactive_indices_res_lna; } { ++evaluations->inactive_indices_res_lna; return timed(evaluations->time.inactive_indices_res_lna, std::bind(&std::remove_cvref_t<Problem>::eval_inactive_indices_res_lna, &problem, γ, x, grad_ψ, J)); }
    real_t eval_f(crvec x) const { ++evaluations->f; return timed(evaluations->time.f, std::bind(&std::remove_cvref_t<Problem>::eval_f, &problem, x)); }
    void eval_grad_f(crvec x, rvec grad_fx) const { ++evaluations->grad_f; return timed(evaluations->time.grad_f, std::bind(&std::remove_cvref_t<Problem>::eval_grad_f, &problem, x, grad_fx)); }
    void eval_g(crvec x, rvec gx) const { ++evaluations->g; return timed(evaluations->time.g, std::bind(&std::remove_cvref_t<Problem>::eval_g, &problem, x, gx)); }
    void eval_grad_g_prod(crvec x, crvec y, rvec grad_gxy) const { ++evaluations->grad_g_prod; return timed(evaluations->time.grad_g_prod, std::bind(&std::remove_cvref_t<Problem>::eval_grad_g_prod, &problem, x, y, grad_gxy)); }
    void eval_grad_gi(crvec x, index_t i, rvec grad_gi) const requires requires { &std::remove_cvref_t<Problem>::eval_grad_gi; } { ++evaluations->grad_gi; return timed(evaluations->time.grad_gi, std::bind(&std::remove_cvref_t<Problem>::eval_grad_gi, &problem, x, i, grad_gi)); }
    void eval_jac_g(crvec x, rindexvec inner_idx, rindexvec outer_ptr, rvec J_values) const requires requires { &std::remove_cvref_t<Problem>::eval_jac_g; } { if (J_values.size() == 0) return problem.eval_jac_g(x, inner_idx, outer_ptr, J_values); ++evaluations->jac_g; return timed(evaluations->time.jac_g, std::bind(&std::remove_cvref_t<Problem>::eval_jac_g, &problem, x, inner_idx, outer_ptr, J_values)); }
    length_t get_jac_g_num_nonzeros() const requires requires { &std::remove_cvref_t<Problem>::get_jac_g_num_nonzeros; } { return problem.get_jac_g_num_nonzeros(); }
    void eval_hess_L_prod(crvec x, crvec y, real_t scale, crvec v, rvec Hv) const requires requires { &std::remove_cvref_t<Problem>::eval_hess_L_prod; } { ++evaluations->hess_L_prod; return timed(evaluations->time.hess_L_prod, std::bind(&std::remove_cvref_t<Problem>::eval_hess_L_prod, &problem, x, y, scale, v, Hv)); }
    void eval_hess_L(crvec x, crvec y, real_t scale, rindexvec inner_idx, rindexvec outer_ptr, rvec H_values) const requires requires { &std::remove_cvref_t<Problem>::eval_hess_L; } { if (H_values.size() == 0) return problem.eval_hess_L(x, y, scale, inner_idx, outer_ptr, H_values); ++evaluations->hess_L; return timed(evaluations->time.hess_L, std::bind(&std::remove_cvref_t<Problem>::eval_hess_L, &problem, x, y, scale, inner_idx, outer_ptr, H_values)); }
    length_t get_hess_L_num_nonzeros() const requires requires { &std::remove_cvref_t<Problem>::get_hess_L_num_nonzeros; } { return problem.get_hess_L_num_nonzeros(); }
    void eval_hess_ψ_prod(crvec x, crvec y, crvec Σ, real_t scale, crvec v, rvec Hv) const requires requires { &std::remove_cvref_t<Problem>::eval_hess_ψ_prod; } { ++evaluations->hess_ψ_prod; return timed(evaluations->time.hess_ψ_prod, std::bind(&std::remove_cvref_t<Problem>::eval_hess_ψ_prod, &problem, x, y, Σ, scale, v, Hv)); }
    void eval_hess_ψ(crvec x, crvec y, crvec Σ, real_t scale, rindexvec inner_idx, rindexvec outer_ptr, rvec H_values) const requires requires { &std::remove_cvref_t<Problem>::eval_hess_ψ; } { if (H_values.size() == 0) return problem.eval_hess_ψ(x, y, Σ, scale, inner_idx, outer_ptr, H_values); ++evaluations->hess_ψ; return timed(evaluations->time.hess_ψ, std::bind(&std::remove_cvref_t<Problem>::eval_hess_ψ, &problem, x, y, Σ, scale, inner_idx, outer_ptr, H_values)); }
    length_t get_hess_ψ_num_nonzeros() const requires requires { &std::remove_cvref_t<Problem>::get_hess_ψ_num_nonzeros; } { return problem.get_hess_ψ_num_nonzeros(); }
    real_t eval_f_grad_f(crvec x, rvec grad_fx) const requires requires { &std::remove_cvref_t<Problem>::eval_f_grad_f; } { ++evaluations->f_grad_f; return timed(evaluations->time.f_grad_f, std::bind(&std::remove_cvref_t<Problem>::eval_f_grad_f, &problem, x, grad_fx)); }
    real_t eval_f_g(crvec x, rvec g) const requires requires { &std::remove_cvref_t<Problem>::eval_f_g; } { ++evaluations->f_g; return timed(evaluations->time.f_g, std::bind(&std::remove_cvref_t<Problem>::eval_f_g, &problem, x, g)); }
    void eval_grad_f_grad_g_prod(crvec x, crvec y, rvec grad_f, rvec grad_gxy) const requires requires { &std::remove_cvref_t<Problem>::eval_grad_f_grad_g_prod; } { ++evaluations->grad_f_grad_g_prod; return timed(evaluations->time.grad_f_grad_g_prod, std::bind(&std::remove_cvref_t<Problem>::eval_grad_f_grad_g_prod, &problem, x, y, grad_f, grad_gxy)); }
    void eval_grad_L(crvec x, crvec y, rvec grad_L, rvec work_n) const requires requires { &std::remove_cvref_t<Problem>::eval_grad_L; } { ++evaluations->grad_L; return timed(evaluations->time.grad_L, std::bind(&std::remove_cvref_t<Problem>::eval_grad_L, &problem, x, y, grad_L, work_n)); }
    real_t eval_ψ(crvec x, crvec y, crvec Σ, rvec ŷ) const requires requires { &std::remove_cvref_t<Problem>::eval_ψ; } { ++evaluations->ψ; return timed(evaluations->time.ψ, std::bind(&std::remove_cvref_t<Problem>::eval_ψ, &problem, x, y, Σ, ŷ)); }
    void eval_grad_ψ(crvec x, crvec y, crvec Σ, rvec grad_ψ, rvec work_n, rvec work_m) const requires requires { &std::remove_cvref_t<Problem>::eval_grad_ψ; } { ++evaluations->grad_ψ; return timed(evaluations->time.grad_ψ, std::bind(&std::remove_cvref_t<Problem>::eval_grad_ψ, &problem, x, y, Σ, grad_ψ, work_n, work_m)); }
    real_t eval_ψ_grad_ψ(crvec x, crvec y, crvec Σ, rvec grad_ψ, rvec work_n, rvec work_m) const requires requires { &std::remove_cvref_t<Problem>::eval_ψ_grad_ψ; } { ++evaluations->ψ_grad_ψ; return timed(evaluations->time.ψ_grad_ψ, std::bind(&std::remove_cvref_t<Problem>::eval_ψ_grad_ψ, &problem, x, y, Σ, grad_ψ, work_n, work_m)); }
    const Box &get_box_C() const requires requires { &std::remove_cvref_t<Problem>::get_box_C; } { return problem.get_box_C(); }
    const Box &get_box_D() const requires requires { &std::remove_cvref_t<Problem>::get_box_D; } { return problem.get_box_D(); }
    void check() const requires requires { &std::remove_cvref_t<Problem>::check; } { return problem.check(); }

    [[nodiscard]] bool provides_eval_grad_gi() const requires requires (Problem p) { { p.provides_eval_grad_gi() } -> std::convertible_to<bool>; } { return problem.provides_eval_grad_gi(); }
    [[nodiscard]] bool provides_eval_inactive_indices_res_lna() const requires requires (Problem p) { { p.provides_eval_inactive_indices_res_lna() } -> std::convertible_to<bool>; } { return problem.provides_eval_inactive_indices_res_lna(); }
    [[nodiscard]] bool provides_eval_jac_g() const requires requires (Problem p) { { p.provides_eval_jac_g() } -> std::convertible_to<bool>; } { return problem.provides_eval_jac_g(); }
    [[nodiscard]] bool provides_get_jac_g_num_nonzeros() const requires requires (Problem p) { { p.provides_get_jac_g_num_nonzeros() } -> std::convertible_to<bool>; } { return problem.provides_get_jac_g_num_nonzeros(); }
    [[nodiscard]] bool provides_eval_hess_L_prod() const requires requires (Problem p) { { p.provides_eval_hess_L_prod() } -> std::convertible_to<bool>; } { return problem.provides_eval_hess_L_prod(); }
    [[nodiscard]] bool provides_eval_hess_L() const requires requires (Problem p) { { p.provides_eval_hess_L() } -> std::convertible_to<bool>; } { return problem.provides_eval_hess_L(); }
    [[nodiscard]] bool provides_get_hess_L_num_nonzeros() const requires requires (Problem p) { { p.provides_get_hess_L_num_nonzeros() } -> std::convertible_to<bool>; } { return problem.provides_get_hess_L_num_nonzeros(); }
    [[nodiscard]] bool provides_eval_hess_ψ_prod() const requires requires (Problem p) { { p.provides_eval_hess_ψ() } -> std::convertible_to<bool>; } { return problem.provides_eval_hess_ψ_prod(); }
    [[nodiscard]] bool provides_eval_hess_ψ() const requires requires (Problem p) { { p.provides_eval_hess_ψ() } -> std::convertible_to<bool>; } { return problem.provides_eval_hess_ψ(); }
    [[nodiscard]] bool provides_get_hess_ψ_num_nonzeros() const requires requires (Problem p) { { p.provides_get_hess_ψ_num_nonzeros() } -> std::convertible_to<bool>; } { return problem.provides_get_hess_ψ_num_nonzeros(); }
    [[nodiscard]] bool provides_eval_f_grad_f() const requires requires (Problem p) { { p.provides_eval_f_grad_f() } -> std::convertible_to<bool>; } { return problem.provides_eval_f_grad_f(); }
    [[nodiscard]] bool provides_eval_f_g() const requires requires (Problem p) { { p.provides_eval_f_g() } -> std::convertible_to<bool>; } { return problem.provides_eval_f_g(); }
    [[nodiscard]] bool provides_eval_grad_f_grad_g_prod() const requires requires (Problem p) { { p.provides_eval_grad_f_grad_g_prod() } -> std::convertible_to<bool>; } { return problem.provides_eval_grad_f_grad_g_prod(); }
    [[nodiscard]] bool provides_eval_grad_L() const requires requires (Problem p) { { p.provides_eval_grad_L() } -> std::convertible_to<bool>; } { return problem.provides_eval_grad_L(); }
    [[nodiscard]] bool provides_eval_ψ() const requires requires (Problem p) { { p.provides_eval_ψ() } -> std::convertible_to<bool>; } { return problem.provides_eval_ψ(); }
    [[nodiscard]] bool provides_eval_grad_ψ() const requires requires (Problem p) { { p.provides_eval_grad_ψ() } -> std::convertible_to<bool>; } { return problem.provides_eval_grad_ψ(); }
    [[nodiscard]] bool provides_eval_ψ_grad_ψ() const requires requires (Problem p) { { p.provides_eval_ψ_grad_ψ() } -> std::convertible_to<bool>; } { return problem.provides_eval_ψ_grad_ψ(); }
    [[nodiscard]] bool provides_get_box_C() const requires requires (Problem p) { { p.provides_get_box_C() } -> std::convertible_to<bool>; } { return problem.provides_get_box_C(); }
    [[nodiscard]] bool provides_get_box_D() const requires requires (Problem p) { { p.provides_get_box_D() } -> std::convertible_to<bool>; } { return problem.provides_get_box_D(); }
    [[nodiscard]] bool provides_check() const requires requires (Problem p) { { p.provides_check() } -> std::convertible_to<bool>; } { return problem.provides_check(); }
    // clang-format on

    [[nodiscard]] length_t get_n() const { return problem.get_n(); }
    [[nodiscard]] length_t get_m() const { return problem.get_m(); }

    std::shared_ptr<EvalCounter> evaluations = std::make_shared<EvalCounter>();
    Problem problem;

    ProblemWithCounters()
        requires std::is_default_constructible_v<Problem>
    = default;
    template <class P>
    explicit ProblemWithCounters(P &&problem)
        requires std::is_same_v<std::remove_cvref_t<P>, std::remove_cvref_t<Problem>>
        : problem{std::forward<P>(problem)} {}
    template <class... Args>
    explicit ProblemWithCounters(std::in_place_t, Args &&...args)
        requires(!std::is_lvalue_reference_v<Problem>)
        : problem{std::forward<Args>(args)...} {}

    /// Reset all evaluation counters and timers to zero. Affects all instances
    /// that share the same evaluations. If you only want to reset the counters
    /// of this instance, use @ref decouple_evaluations first.
    void reset_evaluations() { evaluations.reset(); }
    /// Give this instance its own evaluation counters and timers, decoupling
    /// it from any other instances they might have previously been shared with.
    /// The evaluation counters and timers are preserved (a copy is made).
    void decouple_evaluations() { evaluations = std::make_shared<EvalCounter>(*evaluations); }

  private:
    template <class TimeT, class FunT>
    static decltype(auto) timed(TimeT &time, FunT &&f) {
        util::Timed timed{time};
        return std::forward<FunT>(f)();
    }
};

/// Wraps the given problem into a @ref ProblemWithCounters and keeps track of
/// how many times each function is called, and how long these calls took.
/// The wrapper has its own copy of the given problem. Making copies of the
/// wrapper also copies the underlying problem, but does not copy the evaluation
/// counters, all copies share the same counters.
template <class Problem>
[[nodiscard]] auto problem_with_counters(Problem &&p) {
    using Prob        = std::remove_cvref_t<Problem>;
    using ProbWithCnt = ProblemWithCounters<Prob>;
    return ProbWithCnt{std::forward<Problem>(p)};
}

/// Wraps the given problem into a @ref ProblemWithCounters and keeps track of
/// how many times each function is called, and how long these calls took.
/// The wrapper keeps only a reference to the given problem, it is the
/// responsibility of the caller to make sure that the wrapper does not outlive
/// the original problem. Making copies of the wrapper does not copy the
/// evaluation counters, all copies share the same counters.
template <class Problem>
[[nodiscard]] auto problem_with_counters_ref(Problem &p) {
    using Prob        = std::remove_cvref_t<Problem>;
    using ProbWithCnt = ProblemWithCounters<const Prob &>;
    return ProbWithCnt{p};
}

/// @}

} // namespace alpaqa
