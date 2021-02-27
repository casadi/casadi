#include <cassert>
#include <chrono>
#include <cmath>

#include <panoc-alm/box.hpp>
#include <panoc-alm/lbfgs.hpp>
#include <panoc-alm/panoc.hpp>
#include <panoc-alm/vec.hpp>

#include <iomanip>
#include <iostream>
#include <stdexcept>

namespace pa {

using std::chrono::duration_cast;
using std::chrono::microseconds;

namespace detail {

/// Calculate both ψ(x) and the vector ŷ that can later be used to compute ∇ψ.
/// @f[ \psi(x^k) = f(x^k) + \frac{1}{2}
/// \text{dist}_\Sigma^2\left(g(x^k) + \Sigma^{-1}y,\;D\right) @f]
/// @f[ \hat{y}  @f]
real_t calc_ψ_ŷ(const Problem &p, ///< [in]  Problem description
                const vec &x,     ///< [in]  Decision variable @f$ x @f$
                const vec &y,     ///< [in]  Lagrange multipliers @f$ y @f$
                const vec &Σ,     ///< [in]  Penalty weights @f$ \Sigma @f$
                vec &ŷ            ///< [out] @f$ \hat{y} @f$
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
void calc_grad_ψ_from_ŷ(const Problem &p, ///< [in]  Problem description
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
real_t calc_ψ_grad_ψ(const Problem &p, ///< [in]  Problem description
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
void calc_grad_ψ(const Problem &p, ///< [in]  Problem description
                 const vec &x,     ///< [in]  Decision variable @f$ x @f$
                 const vec &y,     ///< [in]  Lagrange multipliers @f$ y @f$
                 const vec &Σ,     ///< [in]  Penalty weights @f$ \Sigma @f$
                 vec &grad_ψ,      ///< [out] @f$ \nabla \psi(x) @f$
                 vec &work_n,      ///<       Dimension n
                 vec &work_m       ///<       Dimension m
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

/// Calculate ẑ.
/// @f[ \hat{z}^k = \Pi_D\left(g(x^k) + \Sigma^{-1}y\right) @f]
void calc_ẑ(const Problem &p, ///< [in]  Problem description
            const vec &x̂,     ///< [in]  Decision variable @f$ \hat{x} @f$
            const vec &y,     ///< [in]  Lagrange multipliers @f$ y @f$
            const vec &Σ,     ///< [in]  Penalty weights @f$ \Sigma @f$
            vec &ẑ,           ///< [out] @f$ \hat{z} @f$
            vec &err_z        ///< [out] @f$ g(\hat{x}) - \hat{z} @f$
) {
    // g(x̂)
    p.g(x̂, err_z);
    // ζ = g(x̂) + Σ⁻¹y
    ẑ = err_z + (y.array() / Σ.array()).matrix();
    // ẑ = Π(ζ, D)
    ẑ = project(ẑ, p.D);
    // g(x̂) - ẑ
    err_z -= ẑ;
}

/**
 * Projected gradient step
 * @f[ \begin{aligned} 
 * \hat{x}^k &= T_{\gamma^k}\left(x^k\right) \\ 
 * &= \Pi_C\left(x^k - \gamma^k \nabla \psi(x^k)\right) \\ 
 * p^k &= \hat{x}^k - x^k \\ 
 * \end{aligned} @f]
 */
void calc_x̂(const Problem &prob, ///< [in]  Problem description
            real_t γ,            ///< [in]  Step size
            const vec &x,        ///< [in]  Decision variable @f$ x @f$
            const vec &grad_ψ,   ///< [in]  @f$ \nabla \psi(x^k) @f$
            vec &x̂, ///< [out] @f$ \hat{x}^k = T_{\gamma^k}(x^k) @f$
            vec &p  ///< [out] @f$ \hat{x}^k - x^k @f$
) {
    x̂ = project(x - γ * grad_ψ, prob.C);
    p = x̂ - x;
}

/// @f[ \left\| \gamma^{-1} (x - \hat{x}) + \nabla \psi(\hat{x}) -
/// \nabla \psi(x) \right\|_\infty @f]
real_t calc_error_stop_crit(const vec &pₖ, real_t γ, const vec &grad_̂ψₖ,
                            const vec &grad_ψₖ) {
    auto err = (1 / γ) * pₖ + grad_̂ψₖ - grad_ψₖ;
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
    auto start_time              = std::chrono::steady_clock::now();
    unsigned linesearch_failures = 0;
    unsigned lbfgs_failures      = 0;
    unsigned lbfgs_rejected      = 0;

    std::cout << std::scientific;

    using namespace detail;

    const auto n = x.size();
    const auto m = z.size();

    // TODO: allocates
    LBFGS lbfgs;
    SpecializedLBFGS slbfgs;
    if (params.experimental.specialized_lbfgs)
        slbfgs.resize(n, params.lbfgs_mem);
    else
        lbfgs.resize(n, params.lbfgs_mem);

    vec xₖ = x;       // Value of x at the beginning of the iteration
    vec x̂ₖ(n);        // Value of x after a projected gradient step
    vec xₖ₊₁(n);      // xₖ for next iteration
    vec x̂ₖ₊₁(n);      // x̂ₖ for next iteration
    vec ŷx̂ₖ(m);       // Σ (g(x̂ₖ) - ẑₖ)
    vec ŷx̂ₖ₊₁(m);     // ŷ(x̂ₖ) for next iteration
    vec pₖ(n);        // xₖ - x̂ₖ
    vec pₖ₊₁(n);      // xₖ₊₁ - x̂ₖ₊₁
    vec qₖ(n);        // Newton step Hₖ pₖ
    vec grad_ψₖ(n);   // ∇ψ(xₖ)
    vec grad_̂ψₖ(n);   // ∇ψ(x̂ₖ)
    vec grad_ψₖ₊₁(n); // ∇ψ(xₖ₊₁)

    vec work_n(n);
    vec work_m(m);

    // Estimate Lipschitz constant using finite difference
    vec h(n);
    // TODO: what about negative x?
    h = (x * params.Lipschitz.ε).cwiseMax(params.Lipschitz.δ);
    x += h;

    // Calculate ∇ψ(x₀ + h)
    calc_grad_ψ(problem, x, y, Σ, // in
                grad_ψₖ₊₁,        // out
                work_n, work_m);  // work
    // Calculate ψ(xₖ), ∇ψ(x₀)
    real_t ψₖ = calc_ψ_grad_ψ(problem, xₖ, y, Σ, // in
                              grad_ψₖ,           // out
                              work_n, work_m);   // work

    // Estimate Lipschitz constant
    real_t Lₖ = (grad_ψₖ₊₁ - grad_ψₖ).norm() / h.norm();
    // TODO: fail if L is nan?
    if (Lₖ < std::numeric_limits<real_t>::epsilon())
        Lₖ = std::numeric_limits<real_t>::epsilon();
    real_t γₖ                     = params.Lipschitz.Lγ_factor / Lₖ;
    real_t σₖ                     = γₖ * (1 - γₖ * Lₖ) / 2;
    [[maybe_unused]] real_t old_γ = γₖ;

    // Calculate x̂₀, p₀ (projected gradient step)
    calc_x̂(problem, γₖ, xₖ, grad_ψₖ, // in
           x̂ₖ, pₖ);                  // out
    // Calculate ψ(x̂ₖ) and ŷ(x̂ₖ)
    real_t ψ̂xₖ = calc_ψ_ŷ(problem, x̂ₖ, y, Σ, // in
                          ŷx̂ₖ);              // out

    real_t margin = 0; // 1e-6 * std::abs(ψₖ); // TODO: OpEn did this. Why?
    real_t grad_ψₖᵀpₖ = grad_ψₖ.dot(pₖ);
    real_t norm_sq_pₖ = pₖ.squaredNorm();
    // Decrease step size until quadratic upper bound is satisfied
    while (!params.experimental.update_lipschitz_in_linesearch &&
           (ψ̂xₖ > ψₖ + margin + grad_ψₖᵀpₖ + 0.5 * Lₖ * norm_sq_pₖ)) {
        lbfgs.reset();
        Lₖ *= 2;
        σₖ /= 2;
        γₖ /= 2;

        // Calculate x̂ₖ and pₖ (with new step size)
        calc_x̂(problem, γₖ, xₖ, grad_ψₖ, // in
               x̂ₖ, pₖ);                  // out
        // Calculate ∇ψ(xₖ)ᵀpₖ and ‖pₖ‖²
        grad_ψₖᵀpₖ = grad_ψₖ.dot(pₖ);
        norm_sq_pₖ = pₖ.squaredNorm();

        // Calculate ψ(x̂ₖ) and ŷ(x̂ₖ)
        ψ̂xₖ = calc_ψ_ŷ(problem, x̂ₖ, y, Σ, // in
                       ŷx̂ₖ);              // out
    }

    if (params.experimental.specialized_lbfgs)
        slbfgs.initialize(xₖ, grad_ψₖ, x̂ₖ, γₖ);

    for (unsigned k = 0; k <= params.max_iter; ++k) {
        // Calculate ∇ψ(x̂ₖ)
        calc_grad_ψ_from_ŷ(problem, x̂ₖ, ŷx̂ₖ, // in
                           grad_̂ψₖ,          // out
                           work_n);          // work

        // Check stop condition
        real_t εₖ = calc_error_stop_crit(pₖ, γₖ, grad_̂ψₖ, grad_ψₖ);

        // Print progress
        if (k % 100 == 0) {
            std::cout << "[PANOC] " << std::setw(6) << k
                      << ": ψ = " << std::setw(13) << ψₖ
                      << ", ‖∇ψ‖ = " << std::setw(13) << grad_ψₖ.norm()
                      << ", ‖p‖ = " << std::setw(13) << std::sqrt(norm_sq_pₖ)
                      << ", γ = " << std::setw(13) << γₖ
                      << ", εₖ = " << std::setw(13) << εₖ << "\r\n";
        }

        auto time_elapsed = std::chrono::steady_clock::now() - start_time;
        bool out_of_time  = time_elapsed > params.max_time;
        if (εₖ <= ε || k == params.max_iter || out_of_time) {
            // TODO: We could cache g(x) and ẑ, but would that faster?
            //       It saves 1 evaluation of g per ALM iteration, but requires
            //       many extra stores in the inner loops of PANOC.
            // TODO: move the computation of ẑ and g(x) to ALM?
            calc_ẑ(problem, x̂ₖ, y, Σ, z, err_z);
            x = std::move(x̂ₖ);
            y = std::move(ŷx̂ₖ);

            Stats s;
            s.iterations          = k;
            s.ε                   = εₖ;
            s.elapsed_time        = duration_cast<microseconds>(time_elapsed);
            s.status              = εₖ <= ε       ? SolverStatus::Converged
                                    : out_of_time ? SolverStatus::MaxTime
                                                  : SolverStatus::MaxIter;
            s.lbfgs_failures      = lbfgs_failures;
            s.lbfgs_rejected      = lbfgs_rejected;
            s.linesearch_failures = linesearch_failures;
            return s;
        } else if (!std::isfinite(εₖ)) {
            std::cerr << "[PANOC] "
                         "\x1b[0;31m"
                         "inf/NaN"
                         "\x1b[0m"
                      << std::endl;
            std::cout << "[k]   " << k << std::endl;
            std::cout << "qₖ₋₁: " << qₖ.transpose() << std::endl;
            std::cout << "xₖ:   " << xₖ.transpose() << std::endl;
            std::cout << "x̂ₖ:   " << x̂ₖ.transpose() << std::endl;
            std::cout << "ŷx̂ₖ:  " << ŷx̂ₖ.transpose() << std::endl;
            std::cout << "pₖ:   " << pₖ.transpose() << std::endl;
            std::cout << "γₖ:   " << γₖ << std::endl;
            std::cout << "∇_̂ψₖ:  " << grad_̂ψₖ.transpose() << std::endl;
            std::cout << "∇ψₖ:  " << grad_ψₖ.transpose() << std::endl;

            Stats s;
            s.iterations          = k;
            s.ε                   = εₖ;
            s.elapsed_time        = duration_cast<microseconds>(time_elapsed);
            s.status              = SolverStatus::NotFinite;
            s.lbfgs_failures      = lbfgs_failures;
            s.lbfgs_rejected      = lbfgs_rejected;
            s.linesearch_failures = linesearch_failures;
            return s;
        }

        // Calculate Newton step
        if (k > 0) {
            qₖ = pₖ;
            if (params.experimental.specialized_lbfgs)
                slbfgs.apply(qₖ);
            else
                lbfgs.apply(qₖ);
        }

#ifdef PRINT_DEBUG_COUT
        if (old_γ != γₖ) {
            std::cerr << "[PANOC] "
                      << "\x1b[0;34m"
                      << "(" << k << ") "
                      << "Update step size: γ = " << γₖ << "\x1b[0m"
                      << std::endl;
            std::cout << "[PANOC] "
                      << "\x1b[0;34m"
                         "Update L <<<<<<<<"
                         "\x1b[0m"
                      << std::endl;
        }
        std::cout << "[PANOC] "
                  << "L = " << Lₖ << std::endl;
        std::cout << "[PANOC] "
                  << "γ = " << γₖ << std::endl;
        std::cout << "[PANOC] "
                  << "∇ψ(xₖ) = " << grad_ψₖ.transpose() << std::endl;
        std::cout << "[PANOC] "
                  << "xₖ     = " << xₖ.transpose() << std::endl;
        std::cout << "[PANOC] "
                  << "x̂ₖ     = " << x̂ₖ.transpose() << std::endl;
        std::cout << "[PANOC] "
                  << "pₖ     = " << pₖ.transpose() << std::endl;
        std::cout << "[PANOC] "
                  << "qₖ     = " << qₖ.transpose() << std::endl;
        std::cout << "[PANOC] "
                  << "∇ψₖᵀrₖ = " << grad_ψₖᵀpₖ << std::endl;
        std::cout << "[PANOC] "
                  << "‖rₖ‖²  = " << norm_sq_pₖ << std::endl;
        std::cout << "[PANOC] "
                  << "ψ(x̂ₖ)  = " << ψ̂xₖ << std::endl;
        std::cout << "[PANOC] "
                  << "-----------" << std::endl;
#endif

        // Line search
        // TODO: don't we already have this from the previous iteration?
        const real_t φₖ = ψₖ + 1 / (2 * γₖ) * norm_sq_pₖ + grad_ψₖᵀpₖ;
        const real_t σ_norm_γ⁻¹pₖ = σₖ * norm_sq_pₖ / (γₖ * γₖ);
        real_t φₖ₊₁, ψₖ₊₁, grad_ψₖ₊₁ᵀpₖ₊₁, norm_sq_pₖ₊₁;
        real_t τ = 1;
        real_t Lₖ₊₁, σₖ₊₁, γₖ₊₁;
        real_t ls_cond;

        if (k == 0)
            τ = 0;
        if (qₖ.hasNaN()) {
            lbfgs.reset();
            slbfgs.reset();
            std::cerr << "[PANOC] "
                         "\x1b[0;34m"
                         "Newton step NaN"
                         "\x1b[0m"
                      << std::endl;
            τ = 0;
            ++lbfgs_failures;
        }
        do {
            Lₖ₊₁  = Lₖ;
            σₖ₊₁  = σₖ;
            γₖ₊₁  = γₖ;
            old_γ = γₖ;
            if (τ / 2 >= params.τ_min) {
                // Calculate xₖ₊₁
                xₖ₊₁ = xₖ + (1 - τ) * pₖ + τ * qₖ;
            } else {
                xₖ₊₁ = x̂ₖ;
            }
            // Calculate ψ(xₖ₊₁), ∇ψ(xₖ₊₁)
            ψₖ₊₁ = calc_ψ_grad_ψ(problem, xₖ₊₁, y, Σ, // in
                                 grad_ψₖ₊₁,           // out
                                 work_n, work_m);     // work
            // Calculate x̂ₖ₊₁, pₖ₊₁ (projected gradient step)
            calc_x̂(problem, γₖ₊₁, xₖ₊₁, grad_ψₖ₊₁, // in
                   x̂ₖ₊₁, pₖ₊₁);                    // out
            // Calculate ψ(x̂ₖ₊₁) and ŷ(x̂ₖ₊₁)
            real_t ψ̂xₖ₊₁ = calc_ψ_ŷ(problem, x̂ₖ₊₁, y, Σ, // in
                                    ŷx̂ₖ₊₁);              // out

            real_t margin =
                0; // 1e-6 * std::abs(ψₖ); // TODO: OpEn did this. Why?
            grad_ψₖ₊₁ᵀpₖ₊₁ = grad_ψₖ₊₁.dot(pₖ₊₁);
            norm_sq_pₖ₊₁   = pₖ₊₁.squaredNorm();
            // Decrease step size until quadratic upper bound is satisfied
            while (ψ̂xₖ₊₁ > ψₖ₊₁ + margin + grad_ψₖ₊₁ᵀpₖ₊₁ +
                               0.5 * Lₖ₊₁ * norm_sq_pₖ₊₁ &&
                   params.experimental.update_lipschitz_in_linesearch) {
                lbfgs.reset();
                Lₖ₊₁ *= 2;
                σₖ₊₁ /= 2;
                γₖ₊₁ /= 2;

                // Calculate x̂ₖ₊₁ and pₖ₊₁ (with new step size)
                calc_x̂(problem, γₖ₊₁, xₖ₊₁, grad_ψₖ₊₁, // in
                       x̂ₖ₊₁, pₖ₊₁);                    //out
                // Calculate ∇ψ(xₖ₊₁)ᵀpₖ₊₁ and ‖pₖ₊₁‖²
                grad_ψₖ₊₁ᵀpₖ₊₁ = grad_ψₖ₊₁.dot(pₖ₊₁);
                norm_sq_pₖ₊₁   = pₖ₊₁.squaredNorm();
                // Calculate ψ(x̂ₖ₊₁) and ŷ(x̂ₖ₊₁)
                ψ̂xₖ₊₁ = calc_ψ_ŷ(problem, x̂ₖ₊₁, y, Σ, ŷx̂ₖ₊₁);
            }

            φₖ₊₁ = ψₖ₊₁ + 1 / (2 * γₖ₊₁) * norm_sq_pₖ₊₁ + grad_ψₖ₊₁ᵀpₖ₊₁;

            τ /= 2;

            ls_cond = φₖ₊₁ - (φₖ - σ_norm_γ⁻¹pₖ);
        } while (ls_cond > 0 && τ >= params.τ_min);

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

        if (τ < params.τ_min && k != 0) {
            // std::cerr << "[PANOC] "
            //              "\x1b[0;31m"
            //              "Warning: Line search failed"
            //              "\x1b[0m"
            //           << std::endl;
            ++linesearch_failures;
        }

        // Update L-BFGS
        if (params.experimental.specialized_lbfgs)
            lbfgs_rejected +=
                !slbfgs.update(xₖ₊₁, grad_ψₖ₊₁, x̂ₖ₊₁, problem.C, γₖ₊₁);
        else
            lbfgs_rejected += !lbfgs.update(xₖ₊₁ - xₖ, pₖ - pₖ₊₁);

        // Advance step
        Lₖ = Lₖ₊₁;
        σₖ = σₖ₊₁;
        γₖ = γₖ₊₁;

        xₖ = std::move(xₖ₊₁);
        ψₖ = ψₖ₊₁;

        x̂ₖ  = std::move(x̂ₖ₊₁);
        ŷx̂ₖ = std::move(ŷx̂ₖ₊₁);
        pₖ  = std::move(pₖ₊₁);

        grad_ψₖ    = std::move(grad_ψₖ₊₁);
        grad_ψₖᵀpₖ = grad_ψₖ₊₁ᵀpₖ₊₁;
        norm_sq_pₖ = norm_sq_pₖ₊₁;
    }
    throw std::logic_error("[PANOC] loop error");
}

} // namespace pa
