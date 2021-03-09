#pragma once

#include <panoc-alm/problem.hpp>

namespace pa {
namespace detail {

/// Calculate both ψ(x) and the vector ŷ that can later be used to compute ∇ψ.
/// @f[ \psi(x^k) = f(x^k) + \frac{1}{2}
/// \text{dist}_\Sigma^2\left(g(x^k) + \Sigma^{-1}y,\;D\right) @f]
/// @f[ \hat{y}  @f]
inline real_t calc_ψ_ŷ(const Problem &p, ///< [in]  Problem description
                       const vec &x,     ///< [in]  Decision variable @f$ x @f$
                       const vec &y, ///< [in]  Lagrange multipliers @f$ y @f$
                       const vec &Σ, ///< [in]  Penalty weights @f$ \Sigma @f$
                       vec &ŷ        ///< [out] @f$ \hat{y} @f$
) {
    // g(x)
    p.g(x, ŷ);
    // ζ = g(x) + Σ⁻¹y
    ŷ += (y.array() / Σ.array()).matrix();
    // d = ζ - Π(ζ, D)
    ŷ = projecting_difference(ŷ, p.D);
    // dᵀŷ, ŷ = Σ d
    real_t dᵀŷ = 0;
    for (unsigned i = 0; i < p.m; ++i) {
        dᵀŷ += ŷ(i) * Σ(i) * ŷ(i); // TODO: vectorize
        ŷ(i) = Σ(i) * ŷ(i);
    }
    // ψ(x) = f(x) + ½ dᵀŷ
    real_t ψ = p.f(x) + 0.5 * dᵀŷ;

    return ψ;
}

/// Calculate ∇ψ(x) using ŷ.
inline void
calc_grad_ψ_from_ŷ(const Problem &p, ///< [in]  Problem description
                   const vec &x,     ///< [in]  Decision variable @f$ x @f$
                   const vec &ŷ,     ///< [in]  @f$ \hat{y} @f$
                   vec &grad_ψ,      ///< [out] @f$ \nabla \psi(x) @f$
                   vec &work_n       ///<       Dimension n
) {
    // ∇ψ = ∇f(x) + ∇g(x) ŷ
    p.grad_f(x, grad_ψ);
    p.grad_g(x, ŷ, work_n);
    grad_ψ += work_n;
}

/// Calculate both ψ(x) and its gradient ∇ψ(x).
/// @f[ \psi(x^k) = f(x^k) + \frac{1}{2}
/// \text{dist}_\Sigma^2\left(g(x^k) + \Sigma^{-1}y,\;D\right) @f]
/// @f[ \nabla \psi(x) = \nabla f(x) + \nabla g(x)\ \hat{y}(x) @f]
inline real_t
calc_ψ_grad_ψ(const Problem &p, ///< [in]  Problem description
              const vec &x,     ///< [in]  Decision variable @f$ x @f$
              const vec &y,     ///< [in]  Lagrange multipliers @f$ y @f$
              const vec &Σ,     ///< [in]  Penalty weights @f$ \Sigma @f$
              vec &grad_ψ,      ///< [out] @f$ \nabla \psi(x) @f$
              vec &work_n,      ///<       Dimension n
              vec &work_m       ///<       Dimension m
) {
    // ψ(x) = f(x) + ½ dᵀŷ
    real_t ψ = calc_ψ_ŷ(p, x, y, Σ, work_m);
    // ∇ψ = ∇f(x) + ∇g(x) ŷ
    calc_grad_ψ_from_ŷ(p, x, work_m, grad_ψ, work_n);
    return ψ;
}

/// Calculate the gradient ∇ψ(x).
/// @f[ \nabla \psi(x) = \nabla f(x) + \nabla g(x)\ \hat{y}(x) @f]
inline void calc_grad_ψ(const Problem &p, ///< [in]  Problem description
                        const vec &x,     ///< [in]  Decision variable @f$ x @f$
                        const vec &y, ///< [in]  Lagrange multipliers @f$ y @f$
                        const vec &Σ, ///< [in]  Penalty weights @f$ \Sigma @f$
                        vec &grad_ψ,  ///< [out] @f$ \nabla \psi(x) @f$
                        vec &work_n,  ///<       Dimension n
                        vec &work_m   ///<       Dimension m
) {
    // g(x)
    p.g(x, work_m);
    // ζ = g(x) + Σ⁻¹y
    work_m += (y.array() / Σ.array()).matrix();
    // d = ζ - Π(ζ, D)
    work_m = projecting_difference(work_m, p.D);
    // ŷ = Σ d
    work_m = Σ.asDiagonal() * work_m;

    // ∇ψ = ∇f(x) + ∇g(x) ŷ
    p.grad_f(x, grad_ψ);
    p.grad_g(x, work_m, work_n);
    grad_ψ += work_n;
}

/// Calculate the error between ẑ and g(x).
/// @f[ \hat{z}^k = \Pi_D\left(g(x^k) + \Sigma^{-1}y\right) @f]
inline void
calc_err_z(const Problem &p, ///< [in]  Problem description
           const vec &x̂,     ///< [in]  Decision variable @f$ \hat{x} @f$
           const vec &y,     ///< [in]  Lagrange multipliers @f$ y @f$
           const vec &Σ,     ///< [in]  Penalty weights @f$ \Sigma @f$
           vec &err_z        ///< [out] @f$ g(\hat{x}) - \hat{z} @f$
) {
    // g(x̂)
    p.g(x̂, err_z);
    // ζ = g(x̂) + Σ⁻¹y
    // ẑ = Π(ζ, D)
    // g(x) - ẑ
    err_z = err_z - project(err_z + Σ.asDiagonal().inverse() * y, p.D); // TODO: catastrophic cancellation?
}

/**
 * Projected gradient step
 * @f[ \begin{aligned} 
 * \hat{x}^k &= T_{\gamma^k}\left(x^k\right) \\ 
 * &= \Pi_C\left(x^k - \gamma^k \nabla \psi(x^k)\right) \\ 
 * p^k &= \hat{x}^k - x^k \\ 
 * \end{aligned} @f]
 */
inline void calc_x̂(const Problem &prob, ///< [in]  Problem description
                   real_t γ,            ///< [in]  Step size
                   const vec &x,        ///< [in]  Decision variable @f$ x @f$
                   const vec &grad_ψ,   ///< [in]  @f$ \nabla \psi(x^k) @f$
                   vec &x̂, ///< [out] @f$ \hat{x}^k = T_{\gamma^k}(x^k) @f$
                   vec &p  ///< [out] @f$ \hat{x}^k - x^k @f$
) {
    using binary_real_f = real_t (*)(real_t, real_t);
    if (0) { // Naive
        x̂ = project(x - γ * grad_ψ, prob.C);
        p = x̂ - x; // catastrophic cancellation if step is small or x is large
    } else {
        p = (-γ * grad_ψ)
                .binaryExpr(prob.C.lowerbound - x, binary_real_f(std::fmax))
                .binaryExpr(prob.C.upperbound - x, binary_real_f(std::fmin));
        x̂ = x + p;
    }
}

/// @f[ \left\| \gamma_k^{-1} (x^k - \hat x^k) + \nabla \psi(\hat x^k) -
/// \nabla \psi(x^k) \right\|_\infty @f]
inline real_t calc_error_stop_crit(
    const vec &pₖ, ///< [in]  Projected gradient step @f$ \hat x^k - x^k @f$
    real_t γ,      ///< [in]  Step size
    const vec &grad_̂ψₖ, ///< [in]  Gradient in @f$ \hat x^k @f$
    const vec &grad_ψₖ  ///< [in]  Gradient in @f$ x^k @f$
) {
    auto err = (1 / γ) * pₖ + (grad_̂ψₖ - grad_ψₖ);
    // These parentheses     ^^^               ^^^
    // are important to prevent catastrophic cancellation when the step is small
    auto ε = vec_util::norm_inf(err);
    return ε;
}

} // namespace detail
} // namespace pa