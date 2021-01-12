#include <cassert>
#include <cmath>

#include <panoc-alm/lbfgs.hpp>
#include <panoc-alm/panoc.hpp>

#include <iomanip>
#include <iostream>
#include <stdexcept>

namespace pa {

namespace detail {

/// ẑ ← Π(g(x) + Σ⁻¹y, D)
/// @f[ \hat{z}(x) = \Pi_D\left(g(x) + \Sigma^{-1}y\right) @f]
void calc_ẑ(const Problem &p, ///< [in]  Problem description
            const vec &g,     ///< [in]  Constraint value @f$ g(x) @f$
            const vec &Σ⁻¹y,  ///< [in]  @f$ \Sigma^{-1} y @f$
            vec &ẑ            ///< [out] Slack variable @f$ \hat{z}(x) @f$
) {
    ẑ = project(g + Σ⁻¹y, p.D);
}

/// ŷ ← Σ (g(x) - ẑ)
/// @f[ \hat{y}(x) = \Sigma\left(g(x) - \hat{z}(x)\right) + y @f]
void calc_ŷ(const vec &ẑ, ///< [in]  Slack variable @f$ \hat{z}(x) @f$
            const vec &g, ///< [in]  Constraint value @f$ g(x) @f$
            const vec &y, ///< [in]  Lagrange multipliers
            const vec &Σ, ///< [in]  Constraint weights @f$ \Sigma @f$
            vec &ŷ        ///< [out] @f$ \hat{y}(x) @f$
) {
    ŷ = Σ.asDiagonal() * (g - ẑ) + y;
}

/// @f[ \hat{y}(x) = \Sigma\left(g(x) - \hat{z}(x)\right) + y @f]
/// @f$ \hat{z}(x) @f$ is computed implicitly.
void calc_ẑŷ(const Problem &p, ///< [in]  Problem description
             const vec &g,     ///< [in]  Constraint @f$ g(\hat{x}) @f$
             const vec &Σ⁻¹y,  ///< [in]  @f$ \Sigma^{-1} y @f$
             const vec &y,     ///< [in]  Lagrange multipliers
             const vec &Σ,     ///< [in]  Constraint weights @f$ \Sigma @f$
             vec &ŷ            ///< [out] @f$ \hat{y}(x) @f$
) {
    auto ẑ = project(g + Σ⁻¹y, p.D);
    ŷ      = Σ.asDiagonal() * (g - ẑ) + y;
}

/// @f[ \psi(x) = f(x) + \textstyle\frac{1}{2}\displaystyle
/// \text{dist}_\Sigma^2\left(g(x) + \Sigma^{-1}y\right) @f]
real_t calc_ψ(const Problem &p, ///< [in]  Problem description
              const vec &x,     ///< [in]  Decision variable @f$ x @f$
              const vec &g,     ///< [in]  Constraint @f$ g(\hat{x}) @f$
              const vec &Σ⁻¹y,  ///< [in]  @f$ \Sigma^{-1} y @f$
              const vec &ẑ,     ///< [in]  Slack variable @f$ \hat{z}(x) @f$
              const vec &Σ      ///< [in]  Constraint weights @f$ \Sigma @f$
) {
    auto diff    = g + Σ⁻¹y - ẑ; // TODO: this could have been cached
    auto dist_sq = vec_util::norm_squared_weighted(diff, Σ);
    return p.f(x) + 0.5 * dist_sq;
}

/// @f[ \nabla \psi(x) = \nabla f(x) + \nabla g(x)\ \hat{y}(x) @f]
void calc_grad_ψ(
    const Problem &p, ///< [in]  Problem description
    const vec &x,     ///< [in]  Decision variable @f$ x @f$
    const vec &ŷ,     ///< [in]  @f$ \hat{y}(x) = \Sigma\left(g(x) - \Pi_D\left(
                      ///            g(x) + \Sigma^{-1}y\right) + y\right) @f$
    vec &grad_g,      ///< [out] @f$ \nabla g(x) @f$
    vec &grad_ψ       ///< [out] @f$ \nabla \psi(x) @f$
) {
    p.grad_f(x, grad_ψ);    // ∇ψ ← ∇f(x)
    p.grad_g(x, ŷ, grad_g); // ∇g ← ∇g(x) ŷ
    grad_ψ += grad_g;       // ∇ψ ← ∇f(x) + ∇g(x) ŷ
}

/// @f[ \left\| \gamma^{-1} (x - \hat{x}) + \nabla \psi(\hat{x}) -
/// \nabla \psi(x) \right\|_\infty @f]
real_t calc_error_stop_crit(real_t γ, const vec &rₖ, const vec &grad_̂ψₖ,
                            const vec &grad_ψₖ) {
    auto err = (1 / γ) * rₖ + grad_̂ψₖ - grad_ψₖ;
    return vec_util::norm_inf(err);
}

} // namespace detail

PANOCSolver::Stats PANOCSolver::operator()(
    const Problem &problem, ///< [in]    Problem description
    vec &x,                 ///< [inout] Decision variable @f$ x @f$
    vec &z,                 ///< [out]   Slack variable @f$ x @f$
    vec &y,                 ///< [inout] Lagrange multiplier @f$ x @f$
    vec &err_z,             ///< [out]   Slack variable error @f$ g(x) - z @f$
    const vec &Σ,           ///< [in]    Constraint weights @f$ \Sigma @f$
    real_t ε                ///< [in]    Tolerance @f$ \epsilon @f$
) {
    using namespace detail;

    const auto n = x.size();
    const auto m = z.size();

    // TODO: allocates
    LBFGS lbfgs(n, params.lbfgs_mem);

    vec xₖ = x;       // Value of x at the beginning of the iteration
    vec x̂ₖ(n);        // Value of x after a projected gradient step
    vec xₖ₊₁(n);      // xₖ for next iteration
    vec xₖ₋₁(n);      // xₖ of previous iteration (TODO)
    vec rₖ₋₁(n);      // rₖ of previous iteration (TODO)
    vec x̂ₖ₊₁(n);      // x̂ₖ for next iteration
    vec ẑₖ(m);        // ẑ(xₖ) = Π(g(xₖ) + Σ⁻¹y, D)
    vec ẑₖ₊₁(m);      // ẑ(xₖ) for next iteration
    vec ŷₖ(m);        // Σ (g(xₖ) - ẑₖ)
    vec rₖ(n);        // xₖ - x̂ₖ
    vec rₖ_tmp(n);    // Workspace for L-BFGS
    vec rₖ₊₁(n);      // xₖ₊₁ - x̂ₖ₊₁
    vec dₖ(n);        // Newton step Hₖ rₖ
    vec grad_ψₖ(n);   // ∇ψ(xₖ)
    vec grad_̂ψₖ(n);   // ∇ψ(x̂ₖ)
    vec grad_ψₖ₊₁(n); // ∇ψ(xₖ₊₁)
    vec g(m);         // g(x)
    vec grad_g(n);    // ∇g(x)

    // Σ and y are constant in PANOC, so calculate Σ⁻¹y once in advance
    vec Σ⁻¹y = y.array() / Σ.array();

    // Estimate Lipschitz constant using finite difference
    vec h(n);
    // TODO: what about negative x?
    h = (x * params.Lipschitz.ε).cwiseMax(params.Lipschitz.δ);
    x += h;

    // Calculate ∇ψ(x₀ + h)
    problem.g(x, g);
    calc_ẑŷ(problem, g, Σ⁻¹y, y, Σ, ŷₖ);
    calc_grad_ψ(problem, x, ŷₖ, grad_g, grad_ψₖ₊₁);

    // Calculate ẑ(x₀), ∇ψ(x₀)
    problem.g(xₖ, g);
    calc_ẑ(problem, g, Σ⁻¹y, ẑₖ);
    calc_ŷ(ẑₖ, g, y, Σ, ŷₖ);
    calc_grad_ψ(problem, xₖ, ŷₖ, grad_g, grad_ψₖ);

    // Estimate Lipschitz constant
    real_t L = (grad_ψₖ₊₁ - grad_ψₖ).norm() / h.norm();
    if (L < std::numeric_limits<real_t>::epsilon())
        L = std::numeric_limits<real_t>::epsilon();
    real_t γ = params.Lipschitz.Lγ_factor / L;
    real_t σ = γ * (1 - γ * L) / 2;

#ifdef PRINT_DEBUG_COUT
    std::cout << "[PANOC] "
                 "Initial step size: γ = "
              << γ << std::endl;
#endif

    // Calculate x̂₀, r₀ (gradient step)
    x̂ₖ = project(xₖ - γ * grad_ψₖ, problem.C);
    rₖ = xₖ - x̂ₖ;

    // Calculate ψ(x₀), ∇ψ(x₀)ᵀr₀, ‖r₀‖²
    real_t ψₖ         = calc_ψ(problem, x, g, Σ⁻¹y, ẑₖ, Σ);
    real_t grad_ψₖᵀrₖ = grad_ψₖ.dot(rₖ);
    real_t norm_sq_rₖ = rₖ.squaredNorm();

    std::cout << std::scientific;

    // std::cout << "Initial γ = " << γ << std::endl;

    for (unsigned k = 0; k <= params.max_iter; ++k) {
#if 0
        constexpr auto digits = std::numeric_limits<real_t>::digits10 + 1;
        std::cout << "[" << std::setw(3) << k << "] "
                  << "ψ(xₖ) = " << std::setprecision(digits) << ψₖ << std::endl
                  << "      "
                  << "γ     = " << γ << std::endl;
        std::cout << std::endl;
#endif
#ifdef PRINT_DEBUG_COUT
        std::cout << "[PANOC] "
                  << "Iteration #" << k << std::endl;
        std::cout << "[PANOC] "
                  << "xₖ     = " << xₖ.transpose() << std::endl;
        std::cout << "[PANOC] "
                  << "γ      = " << γ << std::endl;
        std::cout << "[PANOC] "
                  << "ψ(xₖ)  = " << ψₖ << std::endl;
        std::cout << "[PANOC] "
                  << "∇ψ(xₖ) = " << grad_ψₖ.transpose() << std::endl;
        problem.g(xₖ, g);
        std::cout << "[PANOC] "
                  << "g(xₖ)  = " << g.transpose() << std::endl;
#endif

        // Calculate g(x̂ₖ), ŷ, ∇ψ(x̂ₖ)
        problem.g(x̂ₖ, g);
        calc_ẑŷ(problem, g, Σ⁻¹y, y, Σ, ŷₖ);
        calc_grad_ψ(problem, x̂ₖ, ŷₖ, grad_g, grad_̂ψₖ);

        // std::cout << "\x1b[0;35m"
        //              "rₖ    = "
        //           << rₖ.transpose() << "\x1b[0m" << std::endl;
        // std::cout << "\x1b[0;35m"
        //              "‖rₖ‖    = "
        //           << rₖ.norm() << "\x1b[0m" << std::endl;
        // std::cout << "\x1b[0;35m"
        //              "grad_̂ψₖ = "
        //           << grad_̂ψₖ.transpose() << "\x1b[0m" << std::endl;
        // std::cout << "\x1b[0;35m"
        //              "grad_ψₖ = "
        //           << grad_ψₖ.transpose() << "\x1b[0m" << std::endl;
        // std::cout << "\x1b[0;35m"
        //              "grad_̂ψₖ - grad_ψₖ = "
        //           << (grad_̂ψₖ - grad_ψₖ).transpose() << "\x1b[0m" << std::endl;
        // std::cout << "\x1b[0;35m"
        //              "‖grad_̂ψₖ - grad_ψₖ‖ = "
        //           << (grad_̂ψₖ - grad_ψₖ).norm() << "\x1b[0m" << std::endl;
        // std::cout << "\x1b[0;35m"
        //              "ε = "
        //           << ε << "\x1b[0m" << std::endl;

        // Check stop condition
        real_t εₖ = calc_error_stop_crit(γ, rₖ, grad_̂ψₖ, grad_ψₖ);
        if (εₖ <= ε || k == params.max_iter) {
            x     = std::move(x̂ₖ);
            z     = std::move(ẑₖ);
            y     = std::move(ŷₖ);
            err_z = g - z;

            Stats s;
            s.iterations = k;
            s.converged  = εₖ <= ε;
            s.failed     = false;
            return s;
        }
        if (!std::isfinite(εₖ)) {
            std::cerr << "[PANOC] "
                         "\x1b[0;31m"
                         "inf/NaN"
                         "\x1b[0m"
                      << std::endl;

            Stats s;
            s.iterations = k;
            s.converged  = false;
            s.failed     = true;
            return s;
        }

        // Calculate ψ(x̂ₖ)
        calc_ẑ(problem, g, Σ⁻¹y, ẑₖ);
        real_t ψ̂xₖ    = calc_ψ(problem, x̂ₖ, g, Σ⁻¹y, ẑₖ, Σ);
        real_t margin = 0; // 1e-6 * std::abs(ψₖ); // TODO: OpEn did this. Why?
        while (ψ̂xₖ > ψₖ + margin - grad_ψₖᵀrₖ + 0.5 * L * norm_sq_rₖ) {
            lbfgs.reset();
            L *= 2;
            σ /= 2;
            γ /= 2;

            // Calculate x̂ₖ and rₖ (with new step size)
            x̂ₖ = project(xₖ - γ * grad_ψₖ, problem.C);
            rₖ = xₖ - x̂ₖ;

            // Calculate ∇ψ(xₖ)ᵀrₖ, ‖rₖ‖²
            grad_ψₖᵀrₖ = grad_ψₖ.dot(rₖ);
            norm_sq_rₖ = rₖ.squaredNorm();

            // Calculate ψ(x̂ₖ)
            problem.g(x̂ₖ, g);
            calc_ẑ(problem, g, Σ⁻¹y, ẑₖ);
            ψ̂xₖ = calc_ψ(problem, x̂ₖ, g, Σ⁻¹y, ẑₖ, Σ);

#ifdef PRINT_DEBUG_COUT
            std::cerr << "[PANOC] "
                      << "\x1b[0;34m"
                      << "(" << k << ") "
                      << "Update step size: γ = " << γ << "\x1b[0m"
                      << std::endl;
            std::cout << "[PANOC] "
                      << "\x1b[0;34m"
                         "Update L <<<<<<<<"
                         "\x1b[0m"
                      << std::endl;
            std::cout << "[PANOC] "
                      << "L = " << L << std::endl;
            std::cout << "[PANOC] "
                      << "γ = " << γ << std::endl;
            std::cout << "[PANOC] "
                      << "∇ψ(xₖ) = " << grad_ψₖ.transpose() << std::endl;
            std::cout << "[PANOC] "
                      << "x̂ₖ     = " << x̂ₖ.transpose() << std::endl;
            std::cout << "[PANOC] "
                      << "rₖ     = " << rₖ.transpose() << std::endl;
            std::cout << "[PANOC] "
                      << "∇ψₖᵀrₖ = " << grad_ψₖᵀrₖ << std::endl;
            std::cout << "[PANOC] "
                      << "‖rₖ‖²  = " << norm_sq_rₖ << std::endl;
            std::cout << "[PANOC] "
                      << "ψ(x̂ₖ)  = " << ψ̂xₖ << std::endl;
            std::cout << "[PANOC] "
                      << "g(x̂ₖ)  = " << g.transpose() << std::endl;
            std::cout << "[PANOC] "
                      << "-----------" << std::endl;
#endif
        }

        // Calculate Newton step
        rₖ_tmp = rₖ;
        lbfgs.apply(1, rₖ_tmp, dₖ);

        // Line search
        real_t φₖ =
            ψₖ - 0.5 * γ * grad_ψₖ.squaredNorm() +
            0.5 / γ * (xₖ - γ * grad_ψₖ - x̂ₖ).squaredNorm(); // TODO: optimize
        real_t σ_norm_γ⁻¹rₖ = σ * norm_sq_rₖ / (γ * γ);
        real_t φₖ₊₁, ψₖ₊₁, grad_ψₖ₊₁ᵀrₖ₊₁, norm_sq_rₖ₊₁;
        real_t τ = 1;
        do {
            if (τ >= params.τ_min) {
                // Calculate xₖ₊₁
                xₖ₊₁ = xₖ - (1 - τ) * rₖ - τ * dₖ;
            } else {
                xₖ₊₁ = x̂ₖ;
            }
            // xₖ₊₁ = project(xₖ₊₁, problem.C);
            // Calculate ẑ(xₖ₊₁), ∇ψ(xₖ₊₁)
            problem.g(xₖ₊₁, g);
            calc_ẑ(problem, g, Σ⁻¹y, ẑₖ₊₁);
            calc_ŷ(ẑₖ₊₁, g, y, Σ, ŷₖ);
            calc_grad_ψ(problem, xₖ₊₁, ŷₖ, grad_g, grad_ψₖ₊₁);
            // Calculate x̂ₖ₊₁, rₖ₊₁ (next gradient step)
            x̂ₖ₊₁ = project(xₖ₊₁ - γ * grad_ψₖ₊₁, problem.C);
            rₖ₊₁ = xₖ₊₁ - x̂ₖ₊₁;

            // Calculate ψ(xₖ₊₁), ‖∇ψ(xₖ₊₁)‖², ‖rₖ₊₁‖²
            ψₖ₊₁ = calc_ψ(problem, xₖ₊₁, g, Σ⁻¹y, ẑₖ₊₁, Σ);
            grad_ψₖ₊₁ᵀrₖ₊₁ = grad_ψₖ₊₁.dot(rₖ₊₁);
            norm_sq_rₖ₊₁   = rₖ₊₁.squaredNorm();
            // Calculate φ(xₖ₊₁)
            // TODO: optimize
            φₖ₊₁ = ψₖ₊₁ - 0.5 * γ * grad_ψₖ₊₁.squaredNorm() +
                   0.5 / γ * (xₖ₊₁ - γ * grad_ψₖ₊₁ - x̂ₖ₊₁).squaredNorm();

            τ /= 2;
        } while ((φₖ₊₁ > φₖ - σ_norm_γ⁻¹rₖ || ψₖ₊₁ > ψₖ) && τ >= params.τ_min);

#if 0
        std::cout << "γ = " << γ << std::endl;
        if (τ * 2 != 1) {
            std::cout << "\x1b[0;35m"
                         "τ = "
                      << τ * 2 << "\x1b[0m" << std::endl;
        } else {
            std::cout << "\x1b[0;33m"
                         "τ = "
                      << τ * 2 << "\x1b[0m" << std::endl;
        }
#endif
        // std::cout << "\x1b[0;33m"
        //              "rₖ = "
        //           << rₖ.transpose() << "\x1b[0m" << std::endl;
        // std::cout << "\x1b[0;33m"
        //              "dₖ = "
        //           << dₖ.transpose() << "\x1b[0m" << std::endl;

        if (τ < params.τ_min) {
            std::cerr << "[PANOC] "
                         "\x1b[0;31m"
                         "Warning: Line search failed"
                         "\x1b[0m"
                      << std::endl;
        }

        // Update L-BFGS
        lbfgs.update(xₖ₊₁ - xₖ, rₖ₊₁ - rₖ);

        // Advance step
        ψₖ         = ψₖ₊₁;
        xₖ₋₁       = std::move(xₖ);
        rₖ₋₁       = std::move(rₖ);
        xₖ         = std::move(xₖ₊₁);
        x̂ₖ         = std::move(x̂ₖ₊₁);
        ẑₖ         = std::move(ẑₖ₊₁);
        rₖ         = std::move(rₖ₊₁);
        grad_ψₖ    = std::move(grad_ψₖ₊₁);
        grad_ψₖᵀrₖ = grad_ψₖ₊₁ᵀrₖ₊₁;
        norm_sq_rₖ = norm_sq_rₖ₊₁;
    }
    throw std::logic_error("[PANOC] loop error");
}

} // namespace pa
