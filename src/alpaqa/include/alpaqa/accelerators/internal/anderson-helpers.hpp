#pragma once

#include <alpaqa/accelerators/internal/limited-memory-qr.hpp>
#include <alpaqa/config/config.hpp>
#include <alpaqa/export.hpp>

namespace alpaqa {

/**
 * @brief   Solve one step of Anderson acceleration to find a fixed point of a 
 *          function g(x):
 * 
 * @f$ g(x^\star) - x^\star = 0 @f$
 * 
 * Updates the QR factorization of @f$ \mathcal{R}_k = QR @f$, solves the least 
 * squares problem to find @f$ \gamma_\text{LS} @f$, computes the next 
 * iterate @f$ x_{k+1} @f$, and stores the current function value @f$ g_k @f$
 * in the matrix @f$ G @f$, which is used as a circular buffer.
 * @f[ \begin{aligned}
 * \def\gammaLS{\gamma_\text{LS}}
 * m_k &= \min \{k, m\} \\
 * g_i &= g(x_i) \\
 * r_i &= r(x_i) g_i - x_i \\
 * \Delta r_i &= r_i - r_{i-1} \\
 * \mathcal{R}_k &= \begin{pmatrix} \Delta r_{k-m_k+1} & \dots & \Delta r_k \end{pmatrix} \in \R^{n\times m_k}
 * \\
 * \gammaLS &= \argmin_{\gamma \in \R^{m_k}}\; \norm{\mathcal{R}_k \gamma - r_k}_2 \\
 * \alpha_i &= \begin{cases} \gammaLS[0] & i = 0 \\
 * \gammaLS[i] - \gammaLS[i-1] & 0 \lt i \lt m_k \\
 * 1 - \gammaLS[m_k - 1] & i = m_k \end{cases} \\
 * \tilde G_k &= \begin{pmatrix} g_{k - m_k} & \dots & g_{k-1} \end{pmatrix} \\
 * G_k &= \begin{pmatrix} g_{k - m_k} & \dots & g_{k} \end{pmatrix} \\
 * &= \begin{pmatrix} \tilde G_k & g_{k} \end{pmatrix} \\
 * x_{k+1} &= \sum_{i=0}^{m_k} \alpha_i\,g_{k - m_k + i} \\
 * &= G_k \alpha \\
 * \end{aligned} @f]
 */
template <Config Conf = DefaultConfig>
inline void minimize_update_anderson(
    /// [inout] QR factorization of @f$ \mathcal{R}_k @f$
    LimitedMemoryQR<Conf> &qr,
    /// [inout] Matrix of previous function values @f$ \tilde G_k @f$
    ///         (stored as ring buffer with the same indices as `qr`)
    rmat<Conf> G̃,
    /// [in]    Current residual @f$ r_k @f$
    crvec<Conf> rₖ,
    /// [in]    Previous residual @f$ r_{k-1} @f$
    crvec<Conf> rₗₐₛₜ,
    /// [in]    Current function value @f$ g_k @f$
    crvec<Conf> gₖ,
    /// [in]    Minimum divisor when solving close to singular systems,
    ///         scaled by the maximum eigenvalue of R
    real_t<Conf> min_div_fac,
    /// [out]   Solution to the least squares system
    rvec<Conf> γ_LS,
    /// [out]   Next Anderson iterate
    rvec<Conf> xₖ_aa) {

    // Update QR factorization for Anderson acceleration
    if (qr.num_columns() == qr.m()) // if the history buffer is full
        qr.remove_column();
    qr.add_column(rₖ - rₗₐₛₜ);

    // Solve least squares problem Anderson acceleration
    // γ = argmin ‖ ΔR γ - rₖ ‖²
    qr.solve_col(rₖ, γ_LS, qr.get_max_eig() * min_div_fac);

    // Iterate over columns of G, whose indices match the indices of the matrix
    // R in the QR factorization, stored as a circular buffer.
    auto g_it  = qr.ring_iter().begin();
    auto g_end = qr.ring_iter().end();
    assert(g_it != g_end);

    // Compute Anderson acceleration next iterate yₑₓₜ = ∑ₙ₌₀ αₙ gₙ
    // α₀ = γ₀             if n = 0
    // αₙ = γₙ - γₙ₋₁      if 0 < n < mₖ
    // αₘ = 1 - γₘ₋₁       if n = mₖ
    auto α = γ_LS(0);
    xₖ_aa  = α * G̃.col((*g_it).circular);
    while (++g_it != g_end) {
        auto [i, g_idx] = *g_it; // [zero based index, circular index]
        α               = γ_LS(i) - γ_LS(i - 1);
        xₖ_aa += α * G̃.col(g_idx);
    }
    α = 1 - γ_LS(qr.num_columns() - 1);
    xₖ_aa += α * gₖ;

    // Add the new column to G
    G̃.col(qr.ring_tail()) = gₖ; // TODO: avoid copy, make G an array of vectors
}

} // namespace alpaqa