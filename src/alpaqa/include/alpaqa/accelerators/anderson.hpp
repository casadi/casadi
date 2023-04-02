#pragma once

#include <alpaqa/accelerators/internal/anderson-helpers.hpp>
#include <alpaqa/config/config.hpp>
#include <alpaqa/export.hpp>

#include <limits>
#include <stdexcept>

namespace alpaqa {

/// Parameters for the @ref AndersonAccel class.
template <Config Conf = DefaultConfig>
struct AndersonAccelParams {
    USING_ALPAQA_CONFIG(Conf);
    /// Length of the history to keep (the number of columns in the QR
    /// factorization).
    /// If this number is greater than the problem dimension, the memory is set
    /// to the problem dimension (otherwise the system is underdetermined).
    length_t memory = 10;
    /// Minimum divisor when solving close to singular systems,
    /// scaled by the maximum eigenvalue of R.
    real_t min_div_fac = real_t(1e2) * std::numeric_limits<real_t>::epsilon();
};

/**
 * Anderson Acceleration.
 *
 * Algorithm for accelerating fixed-point iterations for finding fixed points
 * of a function @f$ g @f$, i.e. @f$ g(x^\star) = x^\star @f$, or equivalently,
 * roots of the residual @f$ r(x) \triangleq g(x) - x @f$.
 *
 * @todo   Condition estimation of the QR factorization.
 *
 * @ingroup grp_Accelerators
 */
template <Config Conf = DefaultConfig>
class AndersonAccel {
  public:
    USING_ALPAQA_CONFIG(Conf);
    using Params = AndersonAccelParams<config_t>;

    AndersonAccel() = default;
    /// @param  params
    ///         Parameters.
    AndersonAccel(Params params) : params(params) {}
    /// @param  params
    ///         Parameters.
    /// @param  n
    ///         Problem dimension (size of the vectors).
    AndersonAccel(Params params, length_t n) : params(params) { resize(n); }

    /// Change the problem dimension. Flushes the history.
    /// @param  n
    ///         Problem dimension (size of the vectors).
    void resize(length_t n) {
        length_t m_AA = std::min(n, params.memory); // TODO: support m > n?
        qr.resize(n, m_AA);
        G.resize(n, m_AA);
        rₗₐₛₜ.resize(n);
        γ_LS.resize(m_AA);
        initialized = false;
    }

    /// Call this function on the first iteration to initialize the accelerator.
    void initialize(crvec g_0, crvec r_0) {
        assert(g_0.size() == n());
        assert(r_0.size() == n());
        G.col(0) = g_0;
        rₗₐₛₜ    = r_0;
        qr.reset();
        initialized = true;
    }

    /// Compute the accelerated iterate @f$ x^k_\text{AA} @f$, given the
    /// function value at the current iterate @f$ g^k = g(x^k) @f$ and the
    /// corresponding residual @f$ r^k = g^k - x^k @f$.
    void compute(crvec gₖ, crvec rₖ, rvec xₖ_aa) {
        if (!initialized)
            throw std::logic_error("AndersonAccel::compute() called before "
                                   "AndersonAccel::initialize()");
        minimize_update_anderson(qr, G,                             // inout
                                 rₖ, rₗₐₛₜ, gₖ, params.min_div_fac, // in
                                 γ_LS, xₖ_aa);                      // out
        rₗₐₛₜ = rₖ;
    }
    /// @copydoc compute(crvec, crvec, rvec)
    void compute(crvec gₖ, vec &&rₖ, rvec xₖ_aa) {
        if (!initialized)
            throw std::logic_error("AndersonAccel::compute() called before "
                                   "AndersonAccel::initialize()");
        minimize_update_anderson(qr, G,                             // inout
                                 rₖ, rₗₐₛₜ, gₖ, params.min_div_fac, // in
                                 γ_LS, xₖ_aa);                      // out
        rₗₐₛₜ = std::move(rₖ);
    }

    /// Reset the accelerator (but keep the last function value and residual, so
    /// calling @ref initialize is not necessary).
    void reset() {
        index_t newest_g_idx = qr.ring_tail();
        if (newest_g_idx != 0)
            G.col(0) = G.col(newest_g_idx);
        qr.reset();
    }

    /// Get the problem dimension.
    length_t n() const { return qr.n(); }
    /// Get the maximum number of stored columns.
    length_t history() const { return qr.m(); }
    /// Get the number of columns currently stored in the buffer.
    length_t current_history() const { return qr.current_history(); }

    /// Get the parameters.
    const Params &get_params() const { return params; }

    std::string get_name() const {
        return "AndersonAccel<" + std::string(config_t::get_name()) + '>';
    }

    /// Scale the factorization
    void scale_R(real_t scal) { qr.scale_R(scal); }

  private:
    Params params;
    LimitedMemoryQR<config_t> qr;
    mat G;
    vec rₗₐₛₜ;
    vec γ_LS;
    bool initialized = false;
};

} // namespace alpaqa