#pragma once

#include <alpaqa/config/config.hpp>
#include <alpaqa/problem/problem-counters.hpp>

#include <memory>

namespace alpaqa {

template <Config Conf  = DefaultConfig,
          class Holder = std::unique_ptr<const ProblemBase<Conf>>>
class WrappedProblemWithCounters : public ProblemBase<Conf> {
  public:
    USING_ALPAQA_CONFIG(Conf);
    using wrapped_type = const ProblemBase<Conf>;
    using holder_type  = Holder;
    using Box          = typename wrapped_type::Box;

    WrappedProblemWithCounters(holder_type p)
        : wrapped_type{p->n, p->m}, problem{std::move(p)} {}

    WrappedProblemWithCounters(const WrappedProblemWithCounters &) = delete;
    WrappedProblemWithCounters &
    operator=(const WrappedProblemWithCounters &)             = delete;
    WrappedProblemWithCounters(WrappedProblemWithCounters &&) = default;
    WrappedProblemWithCounters &
    operator=(WrappedProblemWithCounters &&) = default;

    std::unique_ptr<ProblemBase<config_t>> clone() const & override {
        return std::make_unique<WrappedProblemWithCounters>(
            holder_type{problem->clone().release()});
    }
    std::unique_ptr<ProblemBase<config_t>> clone() && override {
        return std::make_unique<WrappedProblemWithCounters>(std::move(problem));
    }

    real_t eval_f(crvec x) const override;
    void eval_grad_f(crvec x, rvec grad_fx) const override;
    void eval_g(crvec x, rvec gx) const override;
    void eval_grad_g_prod(crvec x, crvec y, rvec grad_gxy) const override;
    void eval_grad_gi(crvec x, index_t i, rvec grad_gi) const override;
    void eval_hess_L_prod(crvec x, crvec y, crvec v, rvec Hv) const override;
    void eval_hess_L(crvec x, crvec y, rmat H) const override;

    real_t eval_f_grad_f(crvec x, rvec grad_fx) const override;
    real_t eval_f_g(crvec x, rvec g) const override;
    real_t eval_f_grad_f_g(crvec x, rvec grad_fx, rvec g) const override;
    void eval_grad_f_grad_g_prod(crvec x, crvec y, rvec grad_f,
                                 rvec grad_gxy) const override;
    void eval_grad_L(crvec x, crvec y, rvec grad_L, rvec work_n) const override;

    real_t eval_ψ_ŷ(crvec x, crvec y, crvec Σ, rvec ŷ) const override;
    void eval_grad_ψ_from_ŷ(crvec x, crvec ŷ, rvec grad_ψ,
                            rvec work_n) const override;
    void eval_grad_ψ(crvec x, crvec y, crvec Σ, rvec grad_ψ, rvec work_n,
                     rvec work_m) const override;
    real_t eval_ψ_grad_ψ(crvec x, crvec y, crvec Σ, rvec grad_ψ, rvec work_n,
                         rvec work_m) const override;

    const Box &get_C() const override;
    const Box &get_D() const override;

    const holder_type &get_problem() const & { return problem; }
    holder_type &&get_problem() && { return std::move(problem); }

    mutable EvalCounter evaluations;

  private:
    template <class TimeT, class FunT>
    static auto timed(TimeT &time, const FunT &f) -> decltype(f()) {
        if constexpr (std::is_same_v<decltype(f()), void>) {
            auto t0 = std::chrono::steady_clock::now();
            f();
            auto t1 = std::chrono::steady_clock::now();
            time += t1 - t0;
        } else {
            auto t0  = std::chrono::steady_clock::now();
            auto res = f();
            auto t1  = std::chrono::steady_clock::now();
            time += t1 - t0;
            return res;
        }
    }
    holder_type problem;
};

template <Config Conf>
WrappedProblemWithCounters<Conf, std::unique_ptr<const ProblemBase<Conf>>>
with_counters(ProblemBase<Conf> &&problem) {
    return std::move(problem).clone();
}

template <Config Conf>
WrappedProblemWithCounters<Conf, const ProblemBase<Conf> *>
with_counters(const ProblemBase<Conf> &problem) {
    return &problem;
}

template <Config Conf>
WrappedProblemWithCounters<Conf, std::unique_ptr<const ProblemBase<Conf>>>
with_counters(std::unique_ptr<ProblemBase<Conf>> problem) {
    return std::move(problem);
}

template <Config Conf>
WrappedProblemWithCounters<Conf, std::shared_ptr<const ProblemBase<Conf>>>
with_counters(std::shared_ptr<ProblemBase<Conf>> problem) {
    return std::move(problem);
}

template <Config Conf, class Holder>
auto WrappedProblemWithCounters<Conf, Holder>::eval_f(crvec x) const -> real_t {
    ++evaluations.f;
    return timed(evaluations.time.f, [&] { return problem->eval_f(x); });
}
template <Config Conf, class Holder>
void WrappedProblemWithCounters<Conf, Holder>::eval_grad_f(crvec x,
                                                           rvec grad_fx) const {
    ++evaluations.grad_f;
    return timed(evaluations.time.grad_f,
                 [&] { return problem->eval_grad_f(x, grad_fx); });
}
template <Config Conf, class Holder>
void WrappedProblemWithCounters<Conf, Holder>::eval_g(crvec x, rvec gx) const {
    ++evaluations.g;
    return timed(evaluations.time.g, [&] { return problem->eval_g(x, gx); });
}
template <Config Conf, class Holder>
void WrappedProblemWithCounters<Conf, Holder>::eval_grad_g_prod(
    crvec x, crvec y, rvec grad_gxy) const {
    ++evaluations.grad_g_prod;
    return timed(evaluations.time.grad_g_prod,
                 [&] { return problem->eval_grad_g_prod(x, y, grad_gxy); });
}
template <Config Conf, class Holder>
void WrappedProblemWithCounters<Conf, Holder>::eval_grad_gi(
    crvec x, index_t i, rvec grad_gi) const {
    ++evaluations.grad_gi;
    return timed(evaluations.time.grad_gi,
                 [&] { return problem->eval_grad_gi(x, i, grad_gi); });
}
template <Config Conf, class Holder>
void WrappedProblemWithCounters<Conf, Holder>::eval_hess_L_prod(crvec x,
                                                                crvec y,
                                                                crvec v,
                                                                rvec Hv) const {
    ++evaluations.hess_L_prod;
    return timed(evaluations.time.hess_L_prod,
                 [&] { return problem->eval_hess_L_prod(x, y, v, Hv); });
}
template <Config Conf, class Holder>
void WrappedProblemWithCounters<Conf, Holder>::eval_hess_L(crvec x, crvec y,
                                                           rmat H) const {
    ++evaluations.hess_L;
    return timed(evaluations.time.hess_L,
                 [&] { return problem->eval_hess_L(x, y, H); });
}

template <Config Conf, class Holder>
typename WrappedProblemWithCounters<Conf, Holder>::real_t
WrappedProblemWithCounters<Conf, Holder>::eval_f_grad_f(crvec x,
                                                        rvec grad_fx) const {
    ++evaluations.f_grad_f;
    return timed(evaluations.time.f_grad_f,
                 [&] { return problem->eval_f_grad_f(x, grad_fx); });
}
template <Config Conf, class Holder>
typename WrappedProblemWithCounters<Conf, Holder>::real_t
WrappedProblemWithCounters<Conf, Holder>::eval_f_g(crvec x, rvec g) const {
    ++evaluations.f_g;
    return timed(evaluations.time.f_g, [&] { return problem->eval_f_g(x, g); });
}
template <Config Conf, class Holder>
typename WrappedProblemWithCounters<Conf, Holder>::real_t
WrappedProblemWithCounters<Conf, Holder>::eval_f_grad_f_g(crvec x, rvec grad_fx,
                                                          rvec g) const {
    ++evaluations.f_grad_f_g;
    return timed(evaluations.time.f_grad_f_g,
                 [&] { return problem->eval_f_grad_f_g(x, grad_fx, g); });
}
template <Config Conf, class Holder>
void WrappedProblemWithCounters<Conf, Holder>::eval_grad_f_grad_g_prod(
    crvec x, crvec y, rvec grad_f, rvec grad_gxy) const {
    ++evaluations.grad_f_grad_g_prod;
    return timed(evaluations.time.grad_f_grad_g_prod, [&] {
        return problem->eval_grad_f_grad_g_prod(x, y, grad_f, grad_gxy);
    });
}
template <Config Conf, class Holder>
void WrappedProblemWithCounters<Conf, Holder>::eval_grad_L(crvec x, crvec y,
                                                           rvec grad_L,
                                                           rvec work_n) const {
    ++evaluations.grad_L;
    return timed(evaluations.time.grad_L,
                 [&] { return problem->eval_grad_L(x, y, grad_L, work_n); });
}
template <Config Conf, class Holder>
auto WrappedProblemWithCounters<Conf, Holder>::eval_ψ_ŷ(crvec x, crvec y,
                                                        crvec Σ, rvec ŷ) const
    -> real_t {
    ++evaluations.ψ;
    return timed(evaluations.time.ψ,
                 [&] { return problem->eval_ψ_ŷ(x, y, Σ, ŷ); });
}
template <Config Conf, class Holder>
void WrappedProblemWithCounters<Conf, Holder>::eval_grad_ψ_from_ŷ(
    crvec x, crvec ŷ, rvec grad_ψ, rvec work_n) const {
    ++evaluations.grad_ψ_from_ŷ;
    return timed(evaluations.time.grad_ψ_from_ŷ, [&] {
        return problem->eval_grad_ψ_from_ŷ(x, ŷ, grad_ψ, work_n);
    });
}
template <Config Conf, class Holder>
void WrappedProblemWithCounters<Conf, Holder>::eval_grad_ψ(crvec x, crvec y,
                                                           crvec Σ, rvec grad_ψ,
                                                           rvec work_n,
                                                           rvec work_m) const {
    ++evaluations.grad_ψ;
    return timed(evaluations.time.grad_ψ, [&] {
        return problem->eval_grad_ψ(x, y, Σ, grad_ψ, work_n, work_m);
    });
}
template <Config Conf, class Holder>
auto WrappedProblemWithCounters<Conf, Holder>::eval_ψ_grad_ψ(
    crvec x, crvec y, crvec Σ, rvec grad_ψ, rvec work_n, rvec work_m) const
    -> real_t {
    ++evaluations.ψ_grad_ψ;
    return timed(evaluations.time.ψ_grad_ψ, [&] {
        return problem->eval_ψ_grad_ψ(x, y, Σ, grad_ψ, work_n, work_m);
    });
}

template <Config Conf, class Holder>
auto WrappedProblemWithCounters<Conf, Holder>::get_C() const -> const Box & {
    return problem->get_C();
}

template <Config Conf, class Holder>
auto WrappedProblemWithCounters<Conf, Holder>::get_D() const -> const Box & {
    return problem->get_D();
}

} // namespace alpaqa
