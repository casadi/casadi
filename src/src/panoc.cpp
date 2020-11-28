#include <cassert>
#include <cmath>

#include <panoc-alm/lbfgs.hpp>
#include <panoc-alm/panoc.hpp>

#include <iostream>

namespace pa {

namespace detail {

void project_y(vec &y,          // inout
               const vec &z_lb, // in
               const vec &z_ub, // in
               real_t M         // in
) {
    constexpr real_t inf = std::numeric_limits<real_t>::infinity();
    // TODO: Handle NaN correctly
    auto max_lb = [M](real_t y, real_t z_lb) {
        real_t y_lb = z_lb == -inf ? 0 : -M;
        return std::max(y, y_lb);
    };
    y = y.binaryExpr(z_lb, max_lb);

    auto min_ub = [M](real_t y, real_t z_ub) {
        real_t y_ub = z_ub == inf ? 0 : M;
        return std::min(y, y_ub);
    };
    y = y.binaryExpr(z_ub, min_ub);
}

} // namespace detail

/**
 * ẑₖ ← Π(g(x̂ₖ) + Σ⁻¹y, D)
 * @f[
 * \hat{z}^k \leftarrow \Pi_D\left(g(\hat{x}^k) + \Sigma^{-1}y\right)
 * @f]
 */
void calc_ẑ(const Problem &p, ///< [in]  Problem description
            const vec &gₖ,    ///< [in]  Constraint @f$ g(\hat{x}^k) @f$
            const vec &Σ⁻¹y,  ///< [in]  @f$ \Sigma^{-1} y @f$
            vec &ẑₖ           ///< [out] Slack variable @f$ \hat{z}^k @f$
) {
    // ẑₖ ← Π(g(x̂ₖ) + Σ⁻¹y, D)
    ẑₖ = project(gₖ + Σ⁻¹y, p.D);
}

/**
 * ŷₖ ← Σ (g(xₖ) - ẑₖ)
 * @f[
 * \hat{y}^k \leftarrow \Sigma\left(g(\hat{x}^k) - \hat{z}\right) + y
 * @f]
 */
void calc_ŷ(const vec &ẑₖ, ///< [in]  Slack variable @f$ \hat{z}^k @f$
            const vec &gₖ, ///< [in]  Constraint @f$ g(\hat{x}^{k+1}) @f$
            const vec &Σ,  ///< [in]  Constraint weights @f$ \Sigma @f$
            vec &ŷₖ        ///< [out] @f$ \hat{y}^k @f$
) {
    ŷₖ = Σ.array() * (gₖ - ẑₖ).array();
    // conversion to Eigen array ensures element-wise multiplication
}

void calc_ẑŷ(const Problem &p, ///< [in]  Problem description
             const vec &gₖ,    ///< [in]  Constraint @f$ g(\hat{x}^k) @f$
             const vec &Σ⁻¹y,  ///< [in]  @f$ \Sigma^{-1} y @f$
             const vec &Σ,     ///< [in]  Constraint weights @f$ \Sigma @f$
             vec &ŷₖ           ///< [out] @f$ \hat{y}^k @f$
) {
    // TODO: does this allocate?
    auto ẑₖ = project(gₖ + Σ⁻¹y, p.D);
    ŷₖ      = Σ.array() * (gₖ - ẑₖ).array();
}

real_t calc_ψ(const Problem &p, const vec &x, const vec &ẑₖ, const vec &Σ) {
    return p.f(x) + 0.5 * dist_squared(ẑₖ, p.D, Σ);
}

void calc_grad_ψ(const Problem &p, ///< [in]  Problem description
                 const vec &xₖ,    ///< [in]  Previous solution @f$ x^k @f$
                 const vec &ŷₖ,    ///< [in]  @f$ \hat{y}^k @f$
                 vec &grad_g,      ///< [out] @f$ \nabla g(x^k) @f$
                 vec &grad_ψₖ      ///< [out] @f$ \nabla \psi(x^k) @f$
) {
    // ∇ψₖ ← ∇f(x)
    p.grad_f(xₖ, grad_ψₖ);
    // ∇gₖ   ← ∇g(x) ŷₖ
    p.grad_g(xₖ, ŷₖ, grad_g);
    // ∇ψₖ₊₁ ← ∇f(x) + ∇gₖ(x) ŷₖ
    grad_ψₖ += grad_g;
}

real_t calc_error_stop_crit(real_t γ, const vec &rₖ, const vec &grad_̂ψₖ,
                            const vec &grad_ψₖ) {
    // TODO: does this allocate?
    auto err = (1 / γ) * rₖ + grad_̂ψₖ - grad_ψₖ;
    return detail::norm_inf(err);
}

void PANOCParams::verify() {
    assert(L > 0); // TODO: better error handling
    assert(γ > 0);
    assert(γ < 1 / L);
    assert(σ > 0);
    assert(σ < γ * (1 - γ * L) / 2);
}

void PANOCSolver::operator()(const Problem &problem, // in
                             vec &x,                 // inout
                             vec &z,                 // out
                             vec &y,                 // inout
                             vec &err_z,             // out
                             const vec &Σ,           // in
                             real_t ε) {
    const auto n = x.size();
    const auto m = z.size();

    // TODO: allocates
    LBFGS lbfgs(n, params.lbfgs_mem);

    vec xₖ = x;       // Value of x at the beginning of the iteration
    vec x̂ₖ(n);        // Value of x after a projected gradient step
    vec xₖ₊₁(n);      // xₖ for next iteration
    vec x̂ₖ₊₁(n);      // x̂ₖ for next iteration
    vec ẑₖ(m);        // ẑ(xₖ) = Π(g(xₖ) + Σ⁻¹y, D)
    vec ẑₖ₊₁(m);      // ẑ(xₖ) for next iteration
    vec ŷₖ(m);        // Σ (g(xₖ) - ẑₖ)
    vec rₖ(n);        // xₖ - x̂ₖ
    vec rₖ_tmp(n);    // Workspace for LBFGS
    vec rₖ₊₁(n);      // xₖ₊₁ - x̂ₖ₊₁
    vec dₖ(n);        // Newton step Hₖ rₖ
    vec grad_ψₖ(n);   // ∇ψ(xₖ)
    vec grad_̂ψₖ(n);   // ∇ψ(x̂ₖ)
    vec grad_ψₖ₊₁(n); // ∇ψ(xₖ₊₁)
    vec g(m);         // g(x)
    vec grad_g(n);    // ∇g(x)

    real_t ψₖ;
    real_t norm_sq_grad_ψₖ;
    real_t norm_sq_rₖ;

    // params.verify();
    real_t L = params.L;
    real_t σ = params.σ;
    real_t γ = params.γ;

    // Σ and y are constant in PANOC, so calculate Σ⁻¹y once in advance
    vec Σ⁻¹y = y.array() / Σ.array();

    // Estimate Lipschitz constant using finite difference
    vec h(n);
    h = (x * params.Lipschitz.ε).cwiseMax(params.Lipschitz.δ);
    x += h;

    // Calculate ∇ψ(x₀ + h)
    problem.g(x, g);
    calc_ẑŷ(problem, g, Σ⁻¹y, Σ, ŷₖ);
    calc_grad_ψ(problem, x, ŷₖ, grad_g, grad_ψₖ₊₁);

    // Calculate ẑ(x₀), ∇ψ(x₀)
    problem.g(xₖ, g);
    calc_ẑ(problem, g, Σ⁻¹y, ẑₖ);
    calc_ŷ(ẑₖ, g, Σ, ŷₖ);
    calc_grad_ψ(problem, xₖ, ŷₖ, grad_g, grad_ψₖ);

    // Estimate Lipschitz constant
    L = (grad_ψₖ₊₁ - grad_ψₖ).norm() / h.norm();
    γ = 0.95 / L;
    σ = γ * (1 - γ * L) / 2;

    // Calculate x̂₀, r₀ (gradient step)
    x̂ₖ = project(xₖ - γ * grad_ψₖ, problem.C);
    rₖ = xₖ - x̂ₖ;

    // Calculate ψ(x₀), ‖∇ψ(x₀)‖², ‖r₀‖²
    ψₖ              = calc_ψ(problem, x, ẑₖ, Σ);
    norm_sq_grad_ψₖ = grad_ψₖ.squaredNorm();
    norm_sq_rₖ      = rₖ.squaredNorm();

    for (unsigned k = 0; k < params.max_iter; ++k) {
        std::cout << "[" << __PRETTY_FUNCTION__ << "] " << k << std::endl;
        // Calculate g(x̂ₖ), ∇ψ(x̂ₖ)
        problem.g(x̂ₖ, g);
        calc_ẑŷ(problem, g, Σ⁻¹y, Σ, ŷₖ);
        calc_grad_ψ(problem, x̂ₖ, ŷₖ, grad_g, grad_̂ψₖ);

        // Check stop condition
        real_t εₖ = calc_error_stop_crit(γ, rₖ, grad_̂ψₖ, grad_ψₖ);
        if (εₖ <= ε) {
            x     = std::move(x̂ₖ);
            z     = std::move(ẑₖ);
            y     = std::move(ŷₖ);
            err_z = g - z;
            return;
        }

        // Calculate ψ(x̂ₖ), ∇ψ(xₖ)ᵀrₖ
        calc_ẑ(problem, g, Σ⁻¹y, ẑₖ);
        real_t ψ̂xₖ        = calc_ψ(problem, x̂ₖ, ẑₖ, Σ);
        real_t grad_ψₖᵀrₖ = grad_ψₖ.dot(rₖ);
        while (ψ̂xₖ > ψₖ - grad_ψₖᵀrₖ + 0.5 * L * rₖ.squaredNorm()) {
            lbfgs.reset();
            L *= 2;
            σ /= 2;
            γ /= 2;

            // Calculate x̂ₖ and rₖ (with new step size)
            x̂ₖ = project(xₖ - γ * grad_ψₖ, problem.C);
            rₖ = xₖ - x̂ₖ;

            // Calculate ∇ψ(xₖ)ᵀrₖ
            grad_ψₖᵀrₖ = grad_ψₖ.dot(rₖ);

            // Calculate ψ(x̂ₖ)
            problem.g(x̂ₖ, g);
            calc_ẑ(problem, g, Σ⁻¹y, ẑₖ);
            ψ̂xₖ = calc_ψ(problem, x̂ₖ, ẑₖ, Σ);
        }

        // Calculate Newton step
        rₖ_tmp = rₖ;
        lbfgs.apply(1, rₖ_tmp, dₖ);

        // Line search
        real_t φₖ = ψₖ - 0.5 * γ * norm_sq_grad_ψₖ + 0.5 / γ * norm_sq_rₖ;
        real_t σ_norm_γ⁻¹rₖ = σ * norm_sq_rₖ / (γ * γ);
        real_t φₖ₊₁, ψₖ₊₁, norm_sq_grad_ψₖ₊₁, norm_sq_rₖ₊₁;
        real_t τ = 1;
        do {
            // Calculate xₖ₊₁
            xₖ₊₁ = xₖ - (1 - τ) * rₖ - τ * dₖ; // TODO: check sign
            // Calculate ẑ(xₖ₊₁), ∇ψ(xₖ₊₁)
            problem.g(xₖ₊₁, g);
            calc_ẑ(problem, g, Σ⁻¹y, ẑₖ₊₁);
            calc_ŷ(ẑₖ₊₁, g, Σ, ŷₖ);
            calc_grad_ψ(problem, xₖ₊₁, ŷₖ, grad_g, grad_ψₖ₊₁);
            // Calculate x̂ₖ₊₁, rₖ₊₁ (next gradient step)
            x̂ₖ₊₁ = project(xₖ₊₁ - γ * grad_ψₖ₊₁, problem.C);
            rₖ₊₁ = xₖ₊₁ - x̂ₖ₊₁;

            // Calculate ψ(xₖ₊₁), ‖∇ψ(xₖ₊₁)‖², ‖rₖ₊₁‖²
            ψₖ₊₁              = calc_ψ(problem, xₖ₊₁, ẑₖ₊₁, Σ);
            norm_sq_grad_ψₖ₊₁ = grad_ψₖ₊₁.squaredNorm();
            norm_sq_rₖ₊₁      = rₖ₊₁.squaredNorm();
            // Calculate φ(xₖ₊₁)
            φₖ₊₁ = ψₖ₊₁ - 0.5 * γ * norm_sq_grad_ψₖ₊₁ + 0.5 / γ * norm_sq_rₖ₊₁;

            τ /= 2;
        } while (φₖ₊₁ > φₖ - σ_norm_γ⁻¹rₖ);

        // Update LBFGS
        lbfgs.update(xₖ₊₁ - xₖ, rₖ₊₁ - rₖ);

        // Advance step
        ψₖ              = ψₖ₊₁;
        x̂ₖ              = std::move(x̂ₖ₊₁);
        ẑₖ              = std::move(ẑₖ₊₁);
        rₖ              = std::move(rₖ₊₁);
        grad_ψₖ         = std::move(grad_ψₖ₊₁);
        norm_sq_grad_ψₖ = norm_sq_grad_ψₖ₊₁;
        norm_sq_rₖ      = norm_sq_rₖ₊₁;
    }
}

} // namespace pa
