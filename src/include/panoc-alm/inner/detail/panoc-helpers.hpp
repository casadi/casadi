#pragma once

#include <panoc-alm/util/atomic_stop_signal.hpp>
#include <panoc-alm/util/problem.hpp>
#include <panoc-alm/util/solverstatus.hpp>

namespace pa::detail {

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
    ŷ += Σ.asDiagonal().inverse() * y;
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
    p.grad_g_prod(x, ŷ, work_n);
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
    p.grad_g_prod(x, work_m, work_n);
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
    err_z = err_z - project(err_z + Σ.asDiagonal().inverse() * y,
                            p.D); // TODO: catastrophic cancellation?
}

/**
 * Projected gradient step
 * @f[ \begin{aligned} 
 * \hat{x}^k &= T_{\gamma^k}\left(x^k\right) \\ 
 * &= \Pi_C\left(x^k - \gamma^k \nabla \psi(x^k)\right) \\ 
 * p^k &= \hat{x}^k - x^k \\ 
 * \end{aligned} @f]
 */
inline auto
projected_gradient_step(const Box &C,     ///< [in]  Set to project onto
                        real_t γ,         ///< [in]  Step size
                        const vec &x,     ///< [in]  Decision variable @f$ x @f$
                        const vec &grad_ψ ///< [in]  @f$ \nabla \psi(x^k) @f$
) {
    using binary_real_f = real_t (*)(real_t, real_t);
    return (-γ * grad_ψ)
        .binaryExpr(C.lowerbound - x, binary_real_f(std::fmax))
        .binaryExpr(C.upperbound - x, binary_real_f(std::fmin));
}

inline void calc_x̂(const Problem &prob, ///< [in]  Problem description
                   real_t γ,            ///< [in]  Step size
                   const vec &x,        ///< [in]  Decision variable @f$ x @f$
                   const vec &grad_ψ,   ///< [in]  @f$ \nabla \psi(x^k) @f$
                   vec &x̂, ///< [out] @f$ \hat{x}^k = T_{\gamma^k}(x^k) @f$
                   vec &p  ///< [out] @f$ \hat{x}^k - x^k @f$
) {
    p = projected_gradient_step(prob.C, γ, x, grad_ψ);
    x̂ = x + p;
}

/// @f[ \left\| \gamma_k^{-1} (x^k - \hat x^k) + \nabla \psi(\hat x^k) -
/// \nabla \psi(x^k) \right\|_\infty @f]
inline real_t calc_error_stop_crit(
    const vec &pₖ, ///< [in]  Projected gradient step @f$ \hat x^k - x^k @f$
    real_t γ,      ///< [in]  Step size
    const vec &grad_̂ψₖ, ///< [in]  Gradient in @f$ \hat x^k @f$
    const vec &grad_ψₖ  ///< [in]  Gradient in @f$ x^k @f$
) {
    auto err = (1 / γ) * pₖ + (grad_ψₖ - grad_̂ψₖ);
    // These parentheses     ^^^               ^^^
    // are important to prevent catastrophic cancellation when the step is small
    auto ε = vec_util::norm_inf(err);
    return ε;
}

/// Increase the estimate of the Lipschitz constant of the objective gradient
/// and decrease the step size until quadratic upper bound or descent lemma is
/// satisfied:
/// @f[ \psi(x + p) \le \psi(x) + \nabla\psi(x)^\top p + \frac{L}{2} \|p\|^2 @f]
/// The projected gradient iterate @f$ \hat x^k @f$ and step @f$ p^k @f$ are
/// updated with the new step size @f$ \gamma^k @f$, and so are other quantities
/// that depend on them, such as @f$ \nabla\psi(x^k)^\top p^k @f$ and
/// @f$ \|p\|^2 @f$. The intermediate vector @f$ \hat y(x^k) @f$ can be used to
/// compute the gradient @f$ \nabla\psi(\hat x^k) @f$ or to update the Lagrange
/// multipliers.
///
/// @return The original step size, before it was reduced by this function.
inline real_t descent_lemma(
    /// [in]  Problem description
    const Problem &problem,
    /// [in]    Tolerance used to ignore rounding errors when the function
    ///         @f$ \psi(x) @f$ is relatively flat or the step size is very
    ///         small, which could cause @f$ \psi(x^k) < \psi(\hat x^k) @f$,
    ///         which is mathematically impossible but could occur in finite
    ///         precision floating point arithmetic.
    real_t rounding_tolerance,
    /// [in]    Current iterate @f$ x^k @f$
    const vec &xₖ,
    /// [in]    Objective function @f$ \psi(x^k) @f$
    real_t ψₖ,
    /// [in]    Gradient of objective @f$ \nabla\psi(x^k) @f$
    const vec &grad_ψₖ,
    /// [in]    Lagrange multipliers @f$ y @f$
    const vec &y,
    /// [in]    Penalty weights @f$ \Sigma @f$
    const vec &Σ,
    /// [out]   Projected gradient iterate @f$ \hat x^k @f$
    vec &x̂ₖ,
    /// [out]   Projected gradient step @f$ p^k @f$
    vec &pₖ,
    /// [out]   Intermediate vector @f$ \hat y(x^k) @f$
    vec &ŷx̂ₖ,
    /// [inout] Objective function @f$ \psi(\hat x^k) @f$
    real_t &ψx̂ₖ,
    /// [inout] Squared norm of the step @f$ \left\| p^k \right\|^2 @f$
    real_t &norm_sq_pₖ,
    /// [inout] Gradient of objective times step @f$ \nabla\psi(x^k)^\top p^k@f$
    real_t &grad_ψₖᵀpₖ,
    /// [inout] Lipschitz constant estimate @f$ L_{\nabla\psi}^k @f$
    real_t &Lₖ,
    /// [inout] Step size @f$ \gamma^k @f$
    real_t &γₖ) {

    real_t old_γₖ = γₖ;
    real_t margin = (1 + std::abs(ψₖ)) * rounding_tolerance;
    while (ψx̂ₖ - ψₖ > grad_ψₖᵀpₖ + 0.5 * Lₖ * norm_sq_pₖ + margin) {
        Lₖ *= 2;
        γₖ /= 2;

        // Calculate x̂ₖ and pₖ (with new step size)
        calc_x̂(problem, γₖ, xₖ, grad_ψₖ, /* in ⟹ out */ x̂ₖ, pₖ);
        // Calculate ∇ψ(xₖ)ᵀpₖ and ‖pₖ‖²
        grad_ψₖᵀpₖ = grad_ψₖ.dot(pₖ);
        norm_sq_pₖ = pₖ.squaredNorm();

        // Calculate ψ(x̂ₖ) and ŷ(x̂ₖ)
        ψx̂ₖ = calc_ψ_ŷ(problem, x̂ₖ, y, Σ, /* in ⟹ out */ ŷx̂ₖ);
    }
    return old_γₖ;
}

/// Check all stop conditions (required tolerance reached, out of time,
/// maximum number of iterations exceeded, interrupted by user,
/// infinite iterate, no progress made)
template <class ParamsT, class DurationT>
inline SolverStatus check_all_stop_conditions(
    /// [in]    Parameters including `max_iter`, `max_time` and `max_no_progress`
    const ParamsT &params,
    /// [in]    Time elapsed since the start of the algorithm
    DurationT time_elapsed,
    /// [in]    The current iteration number
    unsigned iteration,
    /// [in]    A stop signal for the user to interrupt the algorithm
    const AtomicStopSignal &stop_signal,
    /// [in]    Desired primal tolerance
    real_t ε,
    /// [in]    Tolerance of the current iterate
    real_t εₖ,
    /// [in]    The number of successive iterations no progress was made
    unsigned no_progress) {

    bool out_of_time     = time_elapsed > params.max_time;
    bool out_of_iter     = iteration == params.max_iter;
    bool interrupted     = stop_signal.stop_requested();
    bool not_finite      = not std::isfinite(εₖ);
    bool conv            = εₖ <= ε;
    bool max_no_progress = no_progress > params.max_no_progress;
    return conv              ? SolverStatus::Converged
           : out_of_time     ? SolverStatus::MaxTime
           : out_of_iter     ? SolverStatus::MaxIter
           : not_finite      ? SolverStatus::NotFinite
           : max_no_progress ? SolverStatus::NoProgress
           : interrupted     ? SolverStatus::Interrupted
                             : SolverStatus::Unknown;
}

/// Compute the Hessian matrix of the augmented Lagrangian function
/// @f[ \nabla^2_{xx} L_\Sigma(x, y) =
///     \Big. \nabla_{xx}^2 L(x, y) \Big|_{\big(x,\, \hat y(x, y)\big)}
///   + \sum_{i\in\mathcal{I}} \Sigma_i\,\nabla g_i(x) \nabla g_i(x)^\top @f]
inline void calc_augmented_lagrangian_hessian(
    /// [in]  Problem description
    const Problem &problem,
    /// [in]    Current iterate @f$ x^k @f$
    const vec &xₖ,
    /// [in]   Intermediate vector @f$ \hat y(x^k) @f$
    const vec &ŷxₖ,
    /// [in]    Lagrange multipliers @f$ y @f$
    const vec &y,
    /// [in]    Penalty weights @f$ \Sigma @f$
    const vec &Σ,
    /// [out]   The constraint values @f$ g(x^k) @f$
    vec &g,
    /// [out]   Hessian matrix @f$ H(x, y) @f$
    mat &H,
    ///         Dimension n
    vec &work_n) {

    // Compute the Hessian of the Lagrangian
    problem.hess_L(xₖ, ŷxₖ, H);
    // Compute the Hessian of the augmented Lagrangian
    problem.g(xₖ, g);
    for (vec::Index i = 0; i < problem.m; ++i) {
        real_t ζ = g(i) + y(i) / Σ(i);
        bool inactive =
            problem.D.lowerbound(i) < ζ && ζ < problem.D.upperbound(i);
        if (not inactive) {
            problem.grad_gi(xₖ, i, work_n);
            H += work_n * Σ(i) * work_n.transpose();
        }
    }
}

/// Compute the Hessian matrix of the augmented Lagrangian function multiplied
/// by the given vector, using finite differences.
/// @f[ \nabla^2_{xx} L_\Sigma(x, y)\, v \approx
///     \frac{\nabla_x L_\Sigma(x+hv, y) - \nabla_x L_\Sigma(x, y)}{h} @f]
inline void calc_augmented_lagrangian_hessian_prod_fd(
    /// [in]    Problem description
    const Problem &problem,
    /// [in]    Current iterate @f$ x^k @f$
    const vec &xₖ,
    /// [in]    Lagrange multipliers @f$ y @f$
    const vec &y,
    /// [in]    Penalty weights @f$ \Sigma @f$
    const vec &Σ,
    /// [in]    Gradient @f$ \nabla \psi(x^k) @f$
    const vec &grad_ψ,
    /// [in]    Vector to multiply with the Hessian
    const vec &v,
    /// [out]   Hessian-vector product
    vec &Hv,
    ///         Dimension n
    vec &work_n1,
    ///         Dimension n
    vec &work_n2,
    ///         Dimension m
    vec &work_m) {

    real_t h = std::sqrt(std::numeric_limits<real_t>::epsilon());
    vec &xₖh = work_n1;
    xₖh      = xₖ + h * v;
    calc_grad_ψ(problem, xₖh, y, Σ, Hv, work_n2, work_m);
    Hv -= grad_ψ;
    Hv /= h;
}

inline real_t initial_lipschitz_estimate(
    /// [in]    Problem description
    const Problem &problem,
    /// [in]    Current iterate @f$ x^k @f$
    const vec &xₖ,
    /// [in]    Lagrange multipliers @f$ y @f$
    const vec &y,
    /// [in]    Penalty weights @f$ \Sigma @f$
    const vec &Σ,
    /// [in]    Finite difference step size relative to xₖ
    real_t ε,
    /// [in]    Minimum absolute finite difference step size
    real_t δ,
    /// [out]   @f$ \psi(x^k) @f$
    real_t &ψ,
    /// [out]   Gradient @f$ \nabla \psi(x^k) @f$
    vec &grad_ψ,
    ///         Dimension n
    vec &work_n1,
    ///         Dimension n
    vec &work_n2,
    ///         Dimension n
    vec &work_n3,
    ///         Dimension m
    vec &work_m) {

    auto h        = (xₖ * ε).cwiseAbs().cwiseMax(δ);
    work_n1       = xₖ + h;
    real_t norm_h = h.norm();
    // Calculate ∇ψ(x₀ + h)
    calc_grad_ψ(problem, work_n1, y, Σ, /* in ⟹ out */ work_n2, work_n3,
                work_m);
    // Calculate ψ(xₖ), ∇ψ(x₀)
    ψ = calc_ψ_grad_ψ(problem, xₖ, y, Σ, /* in ⟹ out */ grad_ψ, work_n1,
                      work_m);

    // Estimate Lipschitz constant
    real_t L = (work_n2 - grad_ψ).norm() / norm_h;
    if (L < std::numeric_limits<real_t>::epsilon())
        L = std::numeric_limits<real_t>::epsilon();
    return L;
}

} // namespace pa::detail