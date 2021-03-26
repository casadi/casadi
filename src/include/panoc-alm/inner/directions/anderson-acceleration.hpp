#pragma once

#include "panoc-alm/util/vec.hpp"
#include <panoc-alm/inner/detail/limited-memory-qr.hpp>
#include <panoc-alm/inner/directions/decl/panoc-direction-update.hpp>

namespace pa {

inline void minimize_update_anderson(LimitedMemoryQR &qr, mat &G, const vec &rₖ,
                                     const vec &rₖ₋₁, const vec &gₖ, vec &γ_LS,
                                     vec &xₖ_aa) {
    // Update QR factorization for Anderson acceleration
    if (qr.num_columns() == qr.m())
        qr.remove_column();
    qr.add_column(rₖ - rₖ₋₁);

    // Solve least squares problem Anderson acceleration
    qr.solve(rₖ, γ_LS);

    // Compute Anderson acceleration next iterate yₑₓₜ = ∑ₙ₌₀ αₙ gₙ
    size_t g_idx = qr.ring_head();
    // α₀ = γ₀      for n = 0
    real_t α = γ_LS(0);
    xₖ_aa    = α * G.col(g_idx);

    for (size_t n = 1; n < qr.num_columns(); ++n) {
        g_idx = qr.ring_next(g_idx);
        // αₙ = γₙ - γₙ₋₁       for  0 < n < mₖ
        α = γ_LS(n) - γ_LS(n - 1);
        xₖ_aa += α * G.col(g_idx);
    }
    // αₘ = 1 - γₘ₋₁       for n = mₖ
    α = 1 - γ_LS(qr.num_columns() - 1);
    xₖ_aa += α * gₖ;

    // Add the new column to G
    G.col(qr.ring_tail()) = gₖ; // TODO: avoid copy, make G an array of vectors
}

/// Anderson Acceleration.
/// 
/// @ingroup    grp_PANOCDirectionProviders
class AndersonAccel {
  public:
    void resize(size_t n, size_t m) {
        size_t m_AA = std::min(n, m);
        qr.resize(n, m_AA); // TODO: support m > n?
        G.resize(n, m_AA);
        rₖ₋₁.resize(n);
        γ_LS.resize(m_AA);
    }

    void initialize(const vec &g₀, const vec &r₀) {
        G.col(0) = g₀;
        rₖ₋₁     = r₀;
    }

    void compute(const vec &gₖ, const vec &rₖ, vec &xₖ_aa) {
        minimize_update_anderson(qr, G, rₖ, rₖ₋₁, gₖ, γ_LS, xₖ_aa);
    }

    std::string get_name() const { return "AndersonAccel"; }

    void reinitialize(real_t γₖ, real_t old_γₖ) {
        reset();
        // TODO: show mathematically that this is a sensible thing to do:
        rₖ₋₁ *= γₖ / old_γₖ;
    }

    void reset() {
        size_t newest_g_idx = qr.ring_tail();
        if (newest_g_idx != 0)
            G.col(0) = G.col(newest_g_idx);
        qr.reset();
    }

  private:
    LimitedMemoryQR qr;
    mat G;
    vec rₖ₋₁;
    vec γ_LS;
};

template <>
struct PANOCDirection<AndersonAccel> {

    static void initialize(AndersonAccel &aa, const vec &x₀, const vec &x̂₀,
                           const vec &p₀, const vec &grad₀) {
        aa.initialize(x̂₀, p₀);
        (void)aa;
        (void)x₀;
        (void)grad₀;
    }

    static bool update(AndersonAccel &aa, const vec &xₖ, const vec &xₖ₊₁,
                       const vec &pₖ, const vec &pₖ₊₁, const vec &grad_new,
                       const Box &C, real_t γ_new) {
        (void)aa;
        (void)xₖ;
        (void)xₖ₊₁;
        (void)pₖ;
        (void)pₖ₊₁;
        (void)grad_new;
        (void)C;
        (void)γ_new;
        return true;
    }

    static bool apply(AndersonAccel &aa, const vec &xₖ, const vec &x̂ₖ,
                      const vec &pₖ, vec &qₖ) {
        (void)xₖ;
        (void)x̂ₖ;
        qₖ = pₖ;
        aa.compute(x̂ₖ, pₖ, qₖ);
        return true;
    }

    static void changed_γ(AndersonAccel &aa, real_t γₖ, real_t old_γₖ) {
        // TODO: Do we really want to flush every time γ changes?
        aa.reinitialize(γₖ, old_γₖ);
    }
};

} // namespace pa