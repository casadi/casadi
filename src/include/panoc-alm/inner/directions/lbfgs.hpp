#pragma once

#include <panoc-alm/inner/directions/decl/lbfgs.hpp>

namespace pa {

inline bool LBFGS::update_valid(LBFGSParams params, real_t yᵀs, real_t sᵀs,
                                real_t pᵀp) {
    // Smallest number we want to divide by without overflow
    const real_t min_divisor = std::sqrt(std::numeric_limits<real_t>::min());

    // Check if this L-BFGS update is accepted
    if (not std::isfinite(yᵀs))
        return false;
    if (sᵀs < min_divisor)
        return false;
    if (yᵀs < min_divisor)
        return false;

    // CBFGS condition: https://epubs.siam.org/doi/10.1137/S1052623499354242
    real_t α = params.cbfgs.α;
    real_t ϵ = params.cbfgs.ϵ;
    // Condition: yᵀs / sᵀs >= ϵ ‖p‖^α
    bool cbfgs_cond = yᵀs / sᵀs >= ϵ * std::pow(pᵀp, α / 2);
    if (not cbfgs_cond)
        return false;

    return true;
}

inline bool LBFGS::update(const vec &xₖ, const vec &xₖ₊₁, const vec &pₖ,
                          const vec &pₖ₊₁, Sign sign) {
    const auto s = xₖ₊₁ - xₖ;
    const auto y = sign == Sign::Positive ? pₖ₊₁ - pₖ : pₖ - pₖ₊₁;
    real_t yᵀs   = y.dot(s);
    real_t sᵀs   = s.squaredNorm();
    real_t pᵀp   = pₖ₊₁.squaredNorm();
    real_t ρ     = 1 / yᵀs;

    if (not update_valid(params, yᵀs, sᵀs, pᵀp))
        return false;

    // Store the new s and y vectors
    this->s(idx) = s;
    this->y(idx) = y;
    this->ρ(idx) = ρ;

    // Increment the index in the circular buffer
    idx = succ(idx);
    full |= idx == 0;

    return true;
}

template <class Vec>
void LBFGS::apply(Vec &&q, real_t γ) {
    if (idx == 0 && not full)
        return;
    auto new_idx = idx > 0 ? idx - 1 : history() - 1;
    if (γ < 0)
        γ = 1. / (ρ(new_idx) * y(new_idx).squaredNorm());

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

    q *= γ;

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

template <class Vec, class IndexVec>
void LBFGS::apply(Vec &&q, real_t γ, const IndexVec &indices) {
    if (idx == 0 && not full)
        return;
    // Eigen 3.3.9 doesn't yet support indexing using a vector of indices
    // so we'll have to do it manually
    // TODO: abstract this away in an expression template?
    auto dot_ll = [&indices](const auto &a, const auto &b) {
        real_t acc = 0;
        for (auto j : indices)
            acc += a(j) * b(j);
        return acc;
    };
    auto update1 = [&](size_t i) {
        ρ(i) = 1. / dot_ll(s(i), y(i));
        if (ρ(i) <= 0)
            return;
        α(i) = ρ(i) * dot_ll(s(i), q);
        for (auto j : indices)
            q(j) -= α(i) * y(i)(j);
    };
    if (idx)
        for (size_t i = idx; i-- > 0;)
            update1(i);
    if (full)
        for (size_t i = history(); i-- > idx;)
            update1(i);

    for (auto j : indices)
        q(j) *= γ;

    auto update2 = [&](size_t i) {
        if (ρ(i) <= 0)
            return;
        real_t β = ρ(i) * dot_ll(y(i), q);
        for (auto j : indices)
            q(j) += (α(i) - β) * s(i)(j);
    };
    if (full)
        for (size_t i = idx; i < history(); ++i)
            update2(i);
    for (size_t i = 0; i < idx; ++i)
        update2(i);
}

inline void LBFGS::reset() {
    idx  = 0;
    full = false;
}

inline void LBFGS::resize(size_t n, size_t history) {
    sto.resize(n + 1, history * 2);
    reset();
}

inline void LBFGS::scale_y(real_t factor) {
    if (full) {
        for (size_t i = 0; i < history(); ++i) {
            y(i) *= factor;
            ρ(i) *= 1. / factor;
        }
    } else {
        for (size_t i = 0; i < idx; ++i) {
            y(i) *= factor;
            ρ(i) *= 1. / factor;
        }
    }
}

} // namespace pa

#include <panoc-alm/inner/directions/decl/panoc-direction-update.hpp>

namespace pa {

template <>
struct PANOCDirection<LBFGS> {

    static void initialize(LBFGS &lbfgs, const vec &x₀, const vec &x̂₀,
                           const vec &p₀, const vec &grad₀) {
        (void)lbfgs;
        (void)x₀;
        (void)x̂₀;
        (void)p₀;
        (void)grad₀;
    }

    static bool update(LBFGS &lbfgs, const vec &xₖ, const vec &xₖ₊₁,
                       const vec &pₖ, const vec &pₖ₊₁, const vec &grad_new,
                       const Box &C, real_t γ_new) {
        (void)grad_new;
        (void)C;
        (void)γ_new;
        return lbfgs.update(xₖ, xₖ₊₁, pₖ, pₖ₊₁, LBFGS::Sign::Negative);
    }

    static bool apply(LBFGS &lbfgs, const vec &xₖ, const vec &x̂ₖ, const vec &pₖ,
                      real_t γ, vec &qₖ) {
        (void)xₖ;
        (void)x̂ₖ;
        qₖ = pₖ;
        lbfgs.apply(qₖ, γ);
        return true;
    }

    static void changed_γ(LBFGS &lbfgs, real_t γₖ, real_t old_γₖ) {
        if (lbfgs.get_params().rescale_when_γ_changes)
            lbfgs.scale_y(γₖ / old_γₖ);
        else
            lbfgs.reset();
    }
};

} // namespace pa