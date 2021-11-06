#pragma once

#include <alpaqa/inner/directions/decl/lbfgs.hpp>
#include <stdexcept>
#include <type_traits>

namespace alpaqa {

inline bool LBFGS::update_valid(LBFGSParams params, real_t yᵀs, real_t sᵀs,
                                real_t pᵀp) {
    // Smallest number we want to divide by without overflow
    const real_t min_divisor = std::sqrt(std::numeric_limits<real_t>::min());

    // Check if this L-BFGS update is accepted
    if (not std::isfinite(yᵀs))
        return false;
    if (yᵀs < min_divisor)
        return false;
    if (sᵀs < min_divisor)
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

inline bool LBFGS::update(crvec xₖ, crvec xₖ₊₁, crvec pₖ, crvec pₖ₊₁, Sign sign,
                          bool forced) {
    const auto s = xₖ₊₁ - xₖ;
    const auto y = sign == Sign::Positive ? pₖ₊₁ - pₖ : pₖ - pₖ₊₁;
    real_t yᵀs   = y.dot(s);
    real_t ρ     = 1 / yᵀs;
    if (not forced) {
        real_t sᵀs = s.squaredNorm();
        real_t pᵀp = params.cbfgs.ϵ > 0 ? pₖ₊₁.squaredNorm() : 0;
        if (not update_valid(params, yᵀs, sᵀs, pᵀp))
            return false;
    }

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
bool LBFGS::apply(Vec &&q, real_t γ) {
    // Only apply if we have previous vectors s and y
    if (idx == 0 && not full)
        return false;

    // If the step size is negative, compute it as sᵀy/yᵀy
    if (γ < 0) {
        auto new_idx = idx > 0 ? idx - 1 : history() - 1;
        real_t yᵀy   = y(new_idx).squaredNorm();
        γ            = 1. / (ρ(new_idx) * yᵀy);
    }

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

    // r ← H₀ q
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

    return true;
}

template <class Vec, class IndexVec>
bool LBFGS::apply(Vec &&q, real_t γ, const IndexVec &J) {
    // Only apply if we have previous vectors s and y
    if (idx == 0 && not full)
        return false;
    using Index = typename std::remove_reference_t<Vec>::Index;
    bool fullJ  = q.size() == Index(J.size());

    // Eigen 3.3.9 doesn't yet support indexing using a vector of indices
    // so we'll have to do it manually
    // TODO: Abstract this away in an expression template / nullary expression?
    //       Or wait for Eigen update?

    // Dot product of two vectors, adding only the indices in set J
    auto dotJ = [&J, fullJ](const auto &a, const auto &b) {
        if (fullJ) {
            return a.dot(b);
        } else {
            real_t acc = 0;
            for (auto j : J)
                acc += a(j) * b(j);
            return acc;
        }
    };

    auto update1 = [&](size_t i) {
        // Recompute ρ, it depends on the index set J. Note that even if ρ was
        // positive for the full vectors s and y, that's not necessarily the
        // case for the smaller vectors s(J) and y(J).
        if (not fullJ)
            ρ(i) = 1. / dotJ(s(i), y(i));

        if (ρ(i) <= 0) // Reject negative ρ to ensure positive definiteness
            return;

        α(i) = ρ(i) * dotJ(s(i), q);
        if (fullJ)
            q -= α(i) * y(i);
        else
            for (auto j : J)
                q(j) -= α(i) * y(i)(j);

        if (γ < 0) {
            // Compute step size based on most recent yᵀs/yᵀy > 0
            real_t yᵀy = dotJ(y(i), y(i));
            γ          = 1. / (ρ(i) * yᵀy);
        }
    };
    if (idx)
        for (size_t i = idx; i-- > 0;)
            update1(i);
    if (full)
        for (size_t i = history(); i-- > idx;)
            update1(i);

    // If all ρ <= 0, fail
    if (γ < 0)
        return false;

    // r ← H₀ q
    if (fullJ)
        q *= γ;
    else
        for (auto j : J)
            q(j) *= γ;

    auto update2 = [&](size_t i) {
        if (ρ(i) <= 0)
            return;
        real_t β = ρ(i) * dotJ(y(i), q);
        if (fullJ)
            q += (α(i) - β) * s(i);
        else
            for (auto j : J)
                q(j) += (α(i) - β) * s(i)(j);
    };
    if (full)
        for (size_t i = idx; i < history(); ++i)
            update2(i);
    for (size_t i = 0; i < idx; ++i)
        update2(i);

    return true;
}

inline void LBFGS::reset() {
    idx  = 0;
    full = false;
}

inline void LBFGS::resize(size_t n) {
    if (params.memory < 1)
        throw std::invalid_argument("LBFGSParams::memory must be > 1");
    sto.resize(n + 1, params.memory * 2);
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

inline void PANOCDirection<LBFGS>::initialize(crvec x₀, crvec x̂₀, crvec p₀,
                                              crvec grad₀) {
    lbfgs.resize(x₀.size());
    (void)x̂₀;
    (void)p₀;
    (void)grad₀;
}

inline bool PANOCDirection<LBFGS>::update(crvec xₖ, crvec xₖ₊₁, crvec pₖ,
                                          crvec pₖ₊₁, crvec grad_new,
                                          const Box &C, real_t γ_new) {
    (void)grad_new;
    (void)C;
    (void)γ_new;
    return lbfgs.update(xₖ, xₖ₊₁, pₖ, pₖ₊₁, LBFGS::Sign::Negative);
}

inline bool PANOCDirection<LBFGS>::apply(crvec xₖ, crvec x̂ₖ, crvec pₖ, real_t γ,
                                         rvec qₖ) {
    (void)xₖ;
    (void)x̂ₖ;
    qₖ = pₖ;
    return lbfgs.apply(qₖ, γ);
}

inline void PANOCDirection<LBFGS>::changed_γ(real_t γₖ, real_t old_γₖ) {
    if (lbfgs.get_params().rescale_when_γ_changes)
        lbfgs.scale_y(γₖ / old_γₖ);
    else
        lbfgs.reset();
}

inline void PANOCDirection<LBFGS>::reset() { lbfgs.reset(); }

inline std::string PANOCDirection<LBFGS>::get_name() const {
    return lbfgs.get_name();
}

inline LBFGSParams PANOCDirection<LBFGS>::get_params() const {
    return lbfgs.get_params();
}

} // namespace alpaqa