#include "eigen-matchers.hpp"
#include <gtest/gtest.h>

#include <Eigen/QR>
#include <alpaqa/inner/detail/anderson-helpers.hpp>

using alpaqa::crmat;
using alpaqa::crvec;
using alpaqa::mat;
using alpaqa::real_t;
using alpaqa::vec;

mat rotate_add(crmat m, crvec v) {
    mat result(m.rows(), m.cols());
    result.block(0, 0, m.rows(), m.cols() - 1) =
        m.block(0, 1, m.rows(), m.cols() - 1);
    result.rightCols(1) = v;
    return result;
}

TEST(Anderson, minimize) {
    size_t n = 4;
    size_t m = 3;
    size_t K = 7;

    mat R(n, K);
    mat ΔR(n, K - 1);
    mat X(n, K);
    mat G(n, K);
    alpaqa::LimitedMemoryQR qr(4, 3);

    X << 1, 2, 3, 5, 7, 9, -1, //
        2, 4, 3, 7, 1, 2, 4,   //
        1, 3, 2, 1, 2, 1, 3,   //
        5, 2, 1, 3, 3, 4, -1;

    G << 3, 3, 2, 1, 5, 3, 7,  //
        4, 4, 1, 2, 7, 1, -2,  //
        0, 6, 2, 2, -9, -4, 6, //
        0, 8, 1, 2, -1, 1, -1;

    // rₖ = gₖ - xₖ
    for (size_t i = 0; i < K; ++i)
        R.col(i) = G.col(i) - X.col(i);
    // Δrₖ = rₖ₊₁ - rₖ
    for (size_t i = 0; i < K - 1; ++i)
        ΔR.col(i) = R.col(i + 1) - R.col(i);
    // QR factorization of ΔR
    qr.add_column(ΔR.col(0));
    qr.add_column(ΔR.col(1));

    // First iteration ---------------------------------------------------------
    size_t k = 3;

    vec γ_LS(m);
    vec xₖ_aa(n);
    mat G₀ = G.block(0, 0, n, m);
    // Call function under test
    alpaqa::minimize_update_anderson(qr, G₀, R.col(k), R.col(k - 1), G.col(k), γ_LS,
                                 xₖ_aa);

    // Compute reference solution
    // Gref = (g₀ | g₁ | g₂ | g₃)
    mat Gref  = G.block(0, k - m, n, m + 1);
    mat ΔRref = ΔR.block(0, k - m, n, m);
    // γ = argmin ‖ ΔR γ - r₃ ‖²
    vec γ_exp = ΔRref.colPivHouseholderQr().solve(R.col(k));
    // α₀ = γ₀
    // αₙ = γₙ - γₙ₋₁
    // αₘ = 1  - γₘ₋₁
    vec α(m + 1);
    α(0) = γ_exp(0);
    for (size_t i = 1; i < m; ++i)
        α(i) = γ_exp(i) - γ_exp(i - 1);
    α(m) = 1 - γ_exp(m - 1);
    // x = ∑ₙ₌₀ αₙ gₙ
    vec x_exp = Gref * α;

    constexpr real_t ε = 5e-14;
    std::cout << std::setprecision(16) << std::scientific;
    std::cout << γ_LS.transpose() << std::endl;
    std::cout << γ_exp.transpose() << std::endl;
    EXPECT_THAT(print_wrap(γ_LS), EigenAlmostEqual(print_wrap(γ_exp), ε));

    std::cout << xₖ_aa.transpose() << std::endl;
    std::cout << x_exp.transpose() << std::endl;
    EXPECT_THAT(print_wrap(xₖ_aa), EigenAlmostEqual(print_wrap(x_exp), ε));

    // Oldest column of G should have been overwritten by g₃
    mat G₁(4, 3);
    G₁.col(1) = G.col(1);
    G₁.col(2) = G.col(2);
    G₁.col(0) = G.col(3);
    EXPECT_THAT(print_wrap(G₀), EigenAlmostEqual(print_wrap(G₁), ε));

    mat QᵀQ = qr.get_Q().transpose() * qr.get_Q();
    EXPECT_THAT(print_wrap(qr.get_Q() * qr.get_R()),
                EigenAlmostEqual(print_wrap(ΔRref), ε));
    EXPECT_THAT(print_wrap(QᵀQ),
                EigenAlmostEqual(print_wrap(mat::Identity(m, m)), ε));
    std::cout << "\nR:\n" << qr.get_R() << std::endl;
    std::cout << "\nQᵀQ:\n" << QᵀQ << std::endl;
    std::cout << std::endl;

    // Next iteration
    // -------------------------------------------------------------------------

    ++k;

    // Call function under test
    alpaqa::minimize_update_anderson(qr, G₁, R.col(k), R.col(k - 1), G.col(k), γ_LS,
                                 xₖ_aa);

    // Compute reference solution
    // Gref = (g₁ | g₂ | g₃ | g₄)
    Gref  = G.block(0, k - m, n, m + 1);
    ΔRref = ΔR.block(0, k - m, n, m);
    // γ = argmin ‖ ΔR γ - r₃ ‖²
    γ_exp = ΔRref.colPivHouseholderQr().solve(R.col(k));
    // α₀ = γ₀
    // αₙ = γₙ - γₙ₋₁
    // αₘ = 1  - γₘ₋₁
    α(0) = γ_exp(0);
    for (size_t i = 1; i < m; ++i)
        α(i) = γ_exp(i) - γ_exp(i - 1);
    α(m) = 1 - γ_exp(m - 1);
    // x = ∑ₙ₌₀ αₙ gₙ
    x_exp = Gref * α;

    std::cout << γ_LS.transpose() << std::endl;
    std::cout << γ_exp.transpose() << std::endl;
    EXPECT_THAT(print_wrap(γ_LS), EigenAlmostEqual(print_wrap(γ_exp), ε));

    std::cout << xₖ_aa.transpose() << std::endl;
    std::cout << x_exp.transpose() << std::endl;
    EXPECT_THAT(print_wrap(xₖ_aa), EigenAlmostEqual(print_wrap(x_exp), ε));

    // Oldest column of G should have been overwritten by g₄
    mat G₂(4, 3);
    G₂.col(2) = G.col(2);
    G₂.col(0) = G.col(3);
    G₂.col(1) = G.col(4);
    EXPECT_THAT(print_wrap(G₁), EigenAlmostEqual(print_wrap(G₂), ε));

    QᵀQ = qr.get_Q().transpose() * qr.get_Q();
    EXPECT_THAT(print_wrap(qr.get_Q() * qr.get_R()),
                EigenAlmostEqual(print_wrap(ΔRref), ε));
    EXPECT_THAT(print_wrap(QᵀQ),
                EigenAlmostEqual(print_wrap(mat::Identity(m, m)), ε));
    std::cout << "\nR:\n" << qr.get_R() << std::endl;
    std::cout << "\nQᵀQ:\n" << QᵀQ << std::endl;
    std::cout << std::endl;

    // Next iteration
    // -------------------------------------------------------------------------

    ++k;

    // Call function under test
    alpaqa::minimize_update_anderson(qr, G₂, R.col(k), R.col(k - 1), G.col(k), γ_LS,
                                 xₖ_aa);

    // Compute reference solution
    // Gref = (g₂ | g₃ | g₄ | g₅)
    Gref  = G.block(0, k - m, n, m + 1);
    ΔRref = ΔR.block(0, k - m, n, m);
    // γ = argmin ‖ ΔR γ - r₃ ‖²
    γ_exp = ΔRref.colPivHouseholderQr().solve(R.col(k));
    // α₀ = γ₀
    // αₙ = γₙ - γₙ₋₁
    // αₘ = 1  - γₘ₋₁
    α(0) = γ_exp(0);
    for (size_t i = 1; i < m; ++i)
        α(i) = γ_exp(i) - γ_exp(i - 1);
    α(m) = 1 - γ_exp(m - 1);
    // x = ∑ₙ₌₀ αₙ gₙ
    x_exp = Gref * α;

    std::cout << γ_LS.transpose() << std::endl;
    std::cout << γ_exp.transpose() << std::endl;
    EXPECT_THAT(print_wrap(γ_LS), EigenAlmostEqual(print_wrap(γ_exp), ε));

    std::cout << xₖ_aa.transpose() << std::endl;
    std::cout << x_exp.transpose() << std::endl;
    EXPECT_THAT(print_wrap(xₖ_aa), EigenAlmostEqual(print_wrap(x_exp), ε));

    // Oldest column of G should have been overwritten by g₄
    mat G₃(4, 3);
    G₃.col(0) = G.col(3);
    G₃.col(1) = G.col(4);
    G₃.col(2) = G.col(5);
    EXPECT_THAT(print_wrap(G₂), EigenAlmostEqual(print_wrap(G₃), ε));

    QᵀQ = qr.get_Q().transpose() * qr.get_Q();
    EXPECT_THAT(print_wrap(qr.get_Q() * qr.get_R()),
                EigenAlmostEqual(print_wrap(ΔRref), ε));
    EXPECT_THAT(print_wrap(QᵀQ),
                EigenAlmostEqual(print_wrap(mat::Identity(m, m)), ε));
    std::cout << "\nR:\n" << qr.get_R() << std::endl;
    std::cout << "\nQᵀQ:\n" << QᵀQ << std::endl;
    std::cout << std::endl;

    // Next iteration
    // -------------------------------------------------------------------------

    ++k;

    // Call function under test
    alpaqa::minimize_update_anderson(qr, G₃, R.col(k), R.col(k - 1), G.col(k), γ_LS,
                                 xₖ_aa);

    // Compute reference solution
    // Gref = (g₂ | g₃ | g₄ | g₅)
    Gref  = G.block(0, k - m, n, m + 1);
    ΔRref = ΔR.block(0, k - m, n, m);
    // γ = argmin ‖ ΔR γ - r₃ ‖²
    γ_exp = ΔRref.colPivHouseholderQr().solve(R.col(k));
    // α₀ = γ₀
    // αₙ = γₙ - γₙ₋₁
    // αₘ = 1  - γₘ₋₁
    α(0) = γ_exp(0);
    for (size_t i = 1; i < m; ++i)
        α(i) = γ_exp(i) - γ_exp(i - 1);
    α(m) = 1 - γ_exp(m - 1);
    // x = ∑ₙ₌₀ αₙ gₙ
    x_exp = Gref * α;

    std::cout << γ_LS.transpose() << std::endl;
    std::cout << γ_exp.transpose() << std::endl;
    EXPECT_THAT(print_wrap(γ_LS), EigenAlmostEqual(print_wrap(γ_exp), ε));

    std::cout << xₖ_aa.transpose() << std::endl;
    std::cout << x_exp.transpose() << std::endl;
    EXPECT_THAT(print_wrap(xₖ_aa), EigenAlmostEqual(print_wrap(x_exp), ε));

    // Oldest column of G should have been overwritten by g₄
    mat G₄(4, 3);
    G₄.col(1) = G.col(4);
    G₄.col(2) = G.col(5);
    G₄.col(0) = G.col(6);
    EXPECT_THAT(print_wrap(G₃), EigenAlmostEqual(print_wrap(G₄), ε));

    QᵀQ = qr.get_Q().transpose() * qr.get_Q();
    EXPECT_THAT(print_wrap(qr.get_Q() * qr.get_R()),
                EigenAlmostEqual(print_wrap(ΔRref), ε));
    EXPECT_THAT(print_wrap(QᵀQ),
                EigenAlmostEqual(print_wrap(mat::Identity(m, m)), ε));
    std::cout << "\nR:\n" << qr.get_R() << std::endl;
    std::cout << "\nQᵀQ:\n" << QᵀQ << std::endl;
    std::cout << std::endl;
}