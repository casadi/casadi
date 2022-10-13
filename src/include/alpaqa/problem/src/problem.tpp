#pragma once

#include <alpaqa/problem/problem.hpp>

#include <iomanip>
#include <ostream>

#include <alpaqa/util/alloc-check.hpp>

namespace alpaqa {

template <Config Conf>
auto ProblemBase<Conf>::eval_f(crvec) const -> real_t {
    throw not_implemented_error("eval_f");
}
template <Config Conf>
void ProblemBase<Conf>::eval_grad_f(crvec, rvec) const {
    throw not_implemented_error("eval_grad_f");
}
template <Config Conf>
void ProblemBase<Conf>::eval_g(crvec, rvec) const {
    throw not_implemented_error("eval_g");
}
template <Config Conf>
void ProblemBase<Conf>::eval_grad_g_prod(crvec, crvec, rvec) const {
    throw not_implemented_error("eval_grad_g_prod");
}
template <Config Conf>
void ProblemBase<Conf>::eval_grad_gi(crvec, index_t, rvec) const {
    throw not_implemented_error("eval_grad_gi");
}
template <Config Conf>
void ProblemBase<Conf>::eval_hess_L_prod(crvec, crvec, crvec, rvec) const {
    throw not_implemented_error("eval_hess_L_prod");
}
template <Config Conf>
void ProblemBase<Conf>::eval_hess_L(crvec, crvec, rmat) const {
    throw not_implemented_error("eval_hess_L");
}

template <Config Conf>
auto ProblemBase<Conf>::eval_f_grad_f(crvec x, rvec grad_fx) const -> real_t {
    eval_grad_f(x, grad_fx);
    return eval_f(x);
}
template <Config Conf>
auto ProblemBase<Conf>::eval_f_g(crvec x, rvec g) const -> real_t {
    eval_g(x, g);
    return eval_f(x);
}

template <Config Conf>
auto ProblemBase<Conf>::eval_f_grad_f_g(crvec x, rvec grad_fx, rvec g) const
    -> real_t {
    eval_g(x, g);
    return eval_f_grad_f(x, grad_fx);
}

template <Config Conf>
void ProblemBase<Conf>::eval_grad_f_grad_g_prod(crvec x, crvec y, rvec grad_f,
                                                rvec grad_gxy) const {
    eval_grad_f(x, grad_f);
    eval_grad_g_prod(x, y, grad_gxy);
}

template <Config Conf>
void ProblemBase<Conf>::eval_grad_L(crvec x, crvec y, rvec grad_L,
                                    rvec work_n) const {
    // ∇L = ∇f(x) + ∇g(x) y
    eval_grad_f_grad_g_prod(x, y, grad_L, work_n);
    grad_L += work_n;
}

template <Config Conf>
auto ProblemBase<Conf>::eval_ψ_ŷ(crvec x, crvec y, crvec Σ, rvec ŷ) const
    -> real_t {
    if (m == 0) /* [[unlikely]] */
        return eval_f(x);

    real_t f   = eval_f_g(x, ŷ);
    real_t dᵀŷ = calc_ŷ_dᵀŷ(ŷ, y, Σ);
    // ψ(x) = f(x) + ½ dᵀŷ
    real_t ψ = f + real_t(0.5) * dᵀŷ;
    return ψ;
}

template <Config Conf>
void ProblemBase<Conf>::eval_grad_ψ_from_ŷ(crvec x, crvec ŷ, rvec grad_ψ,
                                           rvec work_n) const {
    if (m == 0) /* [[unlikely]] */ {
        eval_grad_f(x, grad_ψ);
    } else {
        eval_grad_L(x, ŷ, grad_ψ, work_n);
    }
}

template <Config Conf>
void ProblemBase<Conf>::eval_grad_ψ(crvec x, crvec y, crvec Σ, rvec grad_ψ,
                                    rvec work_n, rvec work_m) const {
    if (m == 0) /* [[unlikely]] */ {
        eval_grad_f(x, grad_ψ);
    } else {
        eval_g(x, work_m);
        (void)calc_ŷ_dᵀŷ(work_m, y, Σ);
        eval_grad_ψ_from_ŷ(x, work_m, grad_ψ, work_n);
    }
}

#if 0 
// Reference implementation. TODO: verify in tests
template <Config Conf>
auto
ProblemBase<Conf>::eval_ψ_grad_ψ(crvec x, crvec y, crvec Σ, rvec grad_ψ,
                             rvec work_n, rvec work_m) const {
    // ψ(x) = f(x) + ½ dᵀŷ
    real_t ψ = eval_ψ_ŷ(x, y, Σ, work_m);
    // ∇ψ = ∇f(x) + ∇g(x) ŷ
    eval_grad_ψ_from_ŷ(x, work_m, grad_ψ, work_n);
    return ψ;
}
#else
template <Config Conf>
auto ProblemBase<Conf>::eval_ψ_grad_ψ(crvec x, crvec y, crvec Σ, rvec grad_ψ,
                                      rvec work_n, rvec work_m) const
    -> real_t {
    if (m == 0) /* [[unlikely]] */
        return eval_f_grad_f(x, grad_ψ);

    auto &ŷ = work_m;
    // ψ(x) = f(x) + ½ dᵀŷ
    real_t f   = eval_f_g(x, ŷ);
    real_t dᵀŷ = calc_ŷ_dᵀŷ(ŷ, y, Σ);
    real_t ψ   = f + real_t(0.5) * dᵀŷ;
    // ∇ψ(x) = ∇f(x) + ∇g(x) ŷ
    eval_grad_L(x, ŷ, grad_ψ, work_n);
    return ψ;
}
#endif

template <Config Conf>
auto ProblemBase<Conf>::calc_ŷ_dᵀŷ(rvec g_ŷ, crvec y, crvec Σ) const -> real_t {
    if (Σ.size() == 1) {
        // ζ = g(x) + Σ⁻¹y
        g_ŷ += (1 / Σ(0)) * y;
        // d = ζ - Π(ζ, D)
        g_ŷ = projecting_difference(g_ŷ, get_D());
        // dᵀŷ, ŷ = Σ d
        real_t dᵀŷ = Σ(0) * g_ŷ.dot(g_ŷ);
        g_ŷ *= Σ(0);
        return dᵀŷ;
    } else {
        // ζ = g(x) + Σ⁻¹y
        g_ŷ += Σ.asDiagonal().inverse() * y;
        // d = ζ - Π(ζ, D)
        g_ŷ = projecting_difference(g_ŷ, get_D());
        // dᵀŷ, ŷ = Σ d
        real_t dᵀŷ = 0;
        for (unsigned i = 0; i < m; ++i) {
            dᵀŷ += g_ŷ(i) * Σ(i) * g_ŷ(i); // TODO: vectorize
            g_ŷ(i) = Σ(i) * g_ŷ(i);
        }
        return dᵀŷ;
    }
}

template <Config Conf>
std::unique_ptr<ProblemBase<Conf>> Problem<Conf>::clone() const & {
    return std::unique_ptr<Problem>(new Problem(*this));
}
template <Config Conf>
std::unique_ptr<ProblemBase<Conf>> Problem<Conf>::clone() && {
    return std::unique_ptr<Problem>(new Problem(std::move(*this)));
}

// LambdaProblem

template <Config Conf>
auto FunctionalProblem<Conf>::eval_f(crvec x) const -> real_t {
    ScopedMallocAllower ma;
    return f(x);
}
template <Config Conf>
void FunctionalProblem<Conf>::eval_grad_f(crvec x, rvec grad_fx) const {
    ScopedMallocAllower ma;
    return grad_f(x, grad_fx);
}
template <Config Conf>
void FunctionalProblem<Conf>::eval_g(crvec x, rvec gx) const {
    ScopedMallocAllower ma;
    return g(x, gx);
}
template <Config Conf>
void FunctionalProblem<Conf>::eval_grad_g_prod(crvec x, crvec y,
                                               rvec grad_gxy) const {
    ScopedMallocAllower ma;
    return grad_g_prod(x, y, grad_gxy);
}
template <Config Conf>
void FunctionalProblem<Conf>::eval_grad_gi(crvec x, index_t i,
                                           rvec gr_gi) const {
    ScopedMallocAllower ma;
    return grad_gi(x, i, gr_gi);
}
template <Config Conf>
void FunctionalProblem<Conf>::eval_hess_L_prod(crvec x, crvec y, crvec v,
                                               rvec Hv) const {
    ScopedMallocAllower ma;
    return hess_L_prod(x, y, v, Hv);
}
template <Config Conf>
void FunctionalProblem<Conf>::eval_hess_L(crvec x, crvec y, rmat H) const {
    ScopedMallocAllower ma;
    return hess_L(x, y, H);
}

template <Config Conf>
std::unique_ptr<ProblemBase<Conf>> FunctionalProblem<Conf>::clone() const & {
    return std::make_unique<FunctionalProblem>(*this);
}

template <Config Conf>
std::unique_ptr<ProblemBase<Conf>> FunctionalProblem<Conf>::clone() && {
    return std::make_unique<FunctionalProblem>(std::move(*this));
}

} // namespace alpaqa
