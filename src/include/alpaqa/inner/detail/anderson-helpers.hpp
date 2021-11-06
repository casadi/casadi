#pragma once

#include <alpaqa/inner/detail/limited-memory-qr.hpp>

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
 * \mathcal{R}_k &= \begin{pmatrix} \Delta r_{k-m_k} & \dots & \Delta r_{k-1} \end{pmatrix} \in \mathbb{R}^{n\times m_k} \\
 * \Delta r_i &= r_{i+1} - r_i \\
 * r_i &= g_i - x_i \\
 * g_i &= g(x_i) \\
 * \DeclareMathOperator*{\argmin}{argmin}
 * \gamma_\text{LS} &= \argmin_\gamma \left\| \mathcal{R}_k \gamma - r_k \right\|_2 \\
 * \alpha_i &= \begin{cases} \gamma_\text{LS}[0] & i = 0 \\
 *                           \gamma_\text{LS}[i] - \gamma_\text{LS}[i-1] & 0 < i < m_k \\
 *                           1 - \gamma_\text{LS}[m_k - 1] & i = m_k \end{cases} \\
 * x_{k+1} &= \sum_{i=0}^{m_k} \alpha_i\,g_i
 * \end{aligned} @f]
 */
inline void minimize_update_anderson(LimitedMemoryQR &qr, rmat G, crvec rₖ,
                                     crvec rₖ₋₁, crvec gₖ, rvec γ_LS,
                                     rvec xₖ_aa) {
    // Update QR factorization for Anderson acceleration
    if (qr.num_columns() == qr.m()) // if the history buffer is full
        qr.remove_column();
    qr.add_column(rₖ - rₖ₋₁);

    // Solve least squares problem Anderson acceleration
    // γ = argmin ‖ ΔR γ - rₖ ‖²
    qr.solve_col(rₖ, γ_LS);

    // Iterate over columns of G, whose indices match the indices of the matrix
    // R in the QR factorization, stored as a circular buffer.
    auto g_it  = qr.ring_iter().begin();
    auto g_end = qr.ring_iter().end();
    assert(g_it != g_end);

    // Compute Anderson acceleration next iterate yₑₓₜ = ∑ₙ₌₀ αₙ gₙ
    // α₀ = γ₀             if n = 0
    // αₙ = γₙ - γₙ₋₁      if 0 < n < mₖ
    // αₘ = 1 - γₘ₋₁       if n = mₖ
    real_t α = γ_LS(0);
    xₖ_aa    = α * G.col((*g_it).circular);
    while (++g_it != g_end) {
        auto [i, g_idx] = *g_it; // [zero based index, circular index]
        α               = γ_LS(i) - γ_LS(i - 1);
        xₖ_aa += α * G.col(g_idx);
    }
    α = 1 - γ_LS(qr.num_columns() - 1);
    xₖ_aa += α * gₖ;

    // Add the new column to G
    G.col(qr.ring_tail()) = gₖ; // TODO: avoid copy, make G an array of vectors
}

} // namespace alpaqa