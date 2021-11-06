#pragma once

#include <alpaqa/inner/detail/panoc-helpers.hpp>
#include <alpaqa/inner/directions/decl/specialized-lbfgs.hpp>
#include <alpaqa/inner/directions/lbfgs.hpp>

#include <cmath>

namespace alpaqa {

inline void SpecializedLBFGS::initialize(crvec x₀, crvec grad₀) {
    idx  = 0;
    full = false;
    x(0) = x₀;
    g(0) = grad₀;
}

/// Standard L-BFGS update without changing the step size γ.
inline bool SpecializedLBFGS::standard_update(crvec xₖ, crvec xₖ₊₁,
                                              crvec pₖ, crvec pₖ₊₁,
                                              crvec gradₖ₊₁) {
    const auto s = xₖ₊₁ - xₖ;
    const auto y = pₖ - pₖ₊₁;

    real_t yᵀs = y.dot(s);
    real_t sᵀs = s.squaredNorm();
    real_t pᵀp = pₖ₊₁.squaredNorm();
    real_t ρ   = 1 / yᵀs;

    if (not LBFGS::update_valid(params, yᵀs, sᵀs, pᵀp))
        return false;

    // Store the new s and y vectors
    this->s(idx) = s;
    this->y(idx) = y;
    this->ρ(idx) = ρ;

    // Store x and the gradient
    this->x(succ(idx)) = xₖ₊₁;
    this->g(succ(idx)) = gradₖ₊₁;

    // Increment the index in the circular buffer
    idx = succ(idx);
    full |= idx == 0;

    return true;
}

/// L-BFGS update when changing the step size γ, recomputing everything.
inline bool SpecializedLBFGS::full_update(crvec xₖ, crvec xₖ₊₁,
                                          crvec pₖ_old_γ, crvec pₖ₊₁,
                                          crvec gradₖ₊₁, const Box &C,
                                          real_t γ) {
    auto &&sₖ = xₖ₊₁ - xₖ;
    auto &&yₖ = this->w(); // temporary workspace
    // Old pₖ is no longer valid, recompute with new γ
    (void)pₖ_old_γ;
    auto &&pₖ = this->p();
    pₖ        = detail::projected_gradient_step(C, γ, x(idx), g(idx));
    yₖ        = pₖ - pₖ₊₁;

    assert(x(idx) == xₖ);

    real_t yᵀs = yₖ.dot(sₖ);
    real_t sᵀs = sₖ.squaredNorm();
    real_t pᵀp = pₖ₊₁.squaredNorm();
    real_t ρₖ  = 1 / yᵀs;

    if (not LBFGS::update_valid(params, yᵀs, sᵀs, pᵀp))
        return false;

    // Recompute all residuals with new γ
    // yₖ = pₖ - pₖ₊₁
    // pₖ = Π(-γ∇ψ(xₖ), C - xₖ)
    size_t endidx = full ? idx : pred(0);
    for (size_t i = pred(idx); i != endidx; i = pred(i)) {
        this->y(i) = -pₖ /* i+1 */;
        pₖ = detail::projected_gradient_step(C, γ, this->x(i), this->g(i));
        this->y(i) += pₖ /* i */;
    }
    // Store the new s and y vectors
    this->s(idx) = sₖ;
    this->y(idx) = yₖ;
    this->ρ(idx) = ρₖ;

    // Store x and the gradient
    this->x(succ(idx)) = xₖ₊₁;
    this->g(succ(idx)) = gradₖ₊₁;

    // Increment the index in the circular buffer
    idx = succ(idx);
    full |= idx == 0;

    return true;
}

inline bool SpecializedLBFGS::update(crvec xₖ, crvec xₖ₊₁,
                                     crvec pₖ, crvec pₖ₊₁,
                                     crvec gradₖ₊₁, const Box &C,
                                     real_t γ) {
    bool ret = (γ == old_γ || old_γ == 0)
                   ? standard_update(xₖ, xₖ₊₁, pₖ, pₖ₊₁, gradₖ₊₁)
                   : full_update(xₖ, xₖ₊₁, pₖ, pₖ₊₁, gradₖ₊₁, C, γ);
    old_γ    = γ;
    return ret;
}

template <class Vec>
void SpecializedLBFGS::apply(Vec &&q) {
    // TODO: dry, reuse standard LBFGS::apply
    auto update1 = [&](size_t i) {
        α(i) = ρ(i) * (s(i).dot(q));
        q -= α(i) * y(i);
    };
    if (idx)
        for (size_t i = idx; i-- > 0;)
            update1(i);
    if (full)
        for (size_t i = history(); i-- > idx;)
            update1(i);

    // q = H₀ * q; // TODO: diagonal matrix H₀?

    auto update2 = [&](size_t i) {
        real_t β = ρ(i) * (y(i).dot(q));
        q += (α(i) - β) * s(i);
    };
    if (full)
        for (size_t i = idx; i < history(); ++i)
            update2(i);
    for (size_t i = 0; i < idx; ++i)
        update2(i);
}

inline void SpecializedLBFGS::resize(size_t n, size_t history) {
    sto.resize(n + 1, history * 4 + 2);
    sto.fill(std::numeric_limits<real_t>::quiet_NaN());
    idx  = 0;
    full = false;
}

inline void SpecializedLBFGS::reset() {
    x(0) = x(idx);
    g(0) = x(idx);
    idx  = 0;
    full = false;
}

} // namespace alpaqa

#include <alpaqa/inner/directions/decl/panoc-direction-update.hpp>

namespace alpaqa {

template <>
struct PANOCDirection<SpecializedLBFGS> {

    static void initialize(SpecializedLBFGS &lbfgs, crvec x₀,
                           crvec x̂₀, crvec p₀, crvec grad₀) {
        lbfgs.initialize(x₀, grad₀);
        (void)x̂₀;
        (void)p₀;
    }

    static bool update(SpecializedLBFGS &lbfgs, crvec xₖ, crvec xₖ₊₁,
                       crvec pₖ, crvec pₖ₊₁, crvec gradₖ₊₁,
                       const Box &C, real_t γ) {
        return lbfgs.update(xₖ, xₖ₊₁, pₖ, pₖ₊₁, gradₖ₊₁, C, γ);
    }

    static bool apply(SpecializedLBFGS &lbfgs, crvec xₖ, crvec x̂ₖ,
                      crvec pₖ, real_t γ, rvec qₖ) {
        (void)xₖ;
        (void)x̂ₖ;
        (void)γ; // TODO: add this parameter to SLBFGS
        qₖ = pₖ;
        lbfgs.apply(qₖ);
        return true;
    }

    static void changed_γ(SpecializedLBFGS &lbfgs, real_t γₖ, real_t old_γₖ) {
        (void)lbfgs;
        (void)γₖ;
        (void)old_γₖ;
    }
};

} // namespace alpaqa