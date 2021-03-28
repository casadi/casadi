#pragma once

#include <panoc-alm/inner/detail/anderson-helpers.hpp>
#include <panoc-alm/inner/directions/decl/panoc-direction-update.hpp>

namespace pa {

/// Anderson Acceleration.
///
/// @ingroup    grp_PANOCDirectionProviders
class AndersonAccel {
  public:
    void resize(size_t n, size_t m) {
        size_t m_AA = std::min(n, m); // TODO: support m > n?
        qr.resize(n, m_AA);
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