#include <alpaqa/decl/alm.hpp>
#include <alpaqa/inner/decl/panoc.hpp>
#include <alpaqa/inner/directions/decl/lbfgs.hpp>

#include "eigen-matchers.hpp"

TEST(ALM, singleshooting1D) {
    using namespace alpaqa;

    Box C{vec(1), vec(1)};
    C.lowerbound << -1;
    C.upperbound << 1;
    Box D{vec(0), vec(0)};

    using Diag = Eigen::DiagonalMatrix<real_t, Eigen::Dynamic, Eigen::Dynamic>;

    auto A = alpaqa::mat::Constant(1, 1, 0.5);
    auto B = alpaqa::mat::Constant(1, 1, 1);

    auto Q = Diag(1);
    Q.diagonal() << 10;
    auto R = Diag(1);
    R.diagonal() << 1;

    auto x0 = vec(1);
    x0 << 1;

    auto f = [&](crvec u) {
        auto x1 = A * x0 + B * u;
        return x1.dot(Q * x1) + u.dot(R * u);
    };
    auto grad_f = [&](crvec u, rvec grad_f) {
        grad_f = 2 * B * Q * (A * x0 + B * u) + 2 * R * u;
    };
    auto g      = [&](crvec u, rvec g_u) { (void)u, void(g_u); };
    auto grad_g = [&](crvec u, crvec v, rvec grad_u_v) {
        (void)u, (void)v, grad_u_v(0) = 0;
    };

    Problem p{1, 0, C, D, f, grad_f, g, grad_g, {}, {}, {}};

    ALMParams almparam;
    almparam.ε        = 1e-4;
    almparam.δ        = 1e-4;
    almparam.Δ        = 100;
    almparam.Σ₀       = 20;
    almparam.ε₀       = 1;
    almparam.θ        = 0.25;
    almparam.ρ        = 1e-1;
    almparam.M        = 1e9;
    almparam.Σ_max    = 1e9;
    almparam.max_iter = 10;

    PANOCParams panocparam;
    panocparam.Lipschitz.ε = 1e-6;
    panocparam.Lipschitz.δ = 1e-12;
    panocparam.max_iter    = 100;

    LBFGSParams lbfgsparam;
    lbfgsparam.memory = 10;

    ALMSolver<> solver{almparam, {panocparam, lbfgsparam}};

    vec x(1);
    x << 1;
    vec y(0);
    solver(p, y, x);

    EXPECT_NEAR(x(0), -0.454545, 1e-4);
    std::cout << "u = " << x.transpose() << std::endl;
    std::cout << "x = " << (A * x0 + B * x).transpose() << std::endl;
}

TEST(ALM, multipleshooting1D) {
    using namespace alpaqa;

    Box C{vec(2), vec(2)};
    C.lowerbound << -1, -inf;
    C.upperbound << 1, inf;
    Box D{vec(1), vec(1)};
    D.lowerbound << 0;
    D.upperbound << 0;

    using Diag = Eigen::DiagonalMatrix<real_t, Eigen::Dynamic, Eigen::Dynamic>;

    auto A = alpaqa::mat::Constant(1, 1, 0.5);
    auto B = alpaqa::mat::Constant(1, 1, 1);

    auto Q = Diag(1);
    Q.diagonal() << 10;
    auto R = Diag(1);
    R.diagonal() << 1;

    auto x0 = vec(1);
    x0 << 1;

    auto f = [&](crvec ux) {
        auto u = ux.topRows(1);
        auto x = ux.bottomRows(1);
        return x.dot(Q * x) + u.dot(R * u);
    };
    auto grad_f = [&](crvec ux, rvec grad_f) {
        auto u               = ux.topRows(1);
        auto x               = ux.bottomRows(1);
        grad_f.topRows(1)    = 2 * R * u;
        grad_f.bottomRows(1) = 2 * Q * x;
    };
    auto g = [&](crvec ux, rvec g_u) {
        auto u         = ux.topRows(1);
        auto x         = ux.bottomRows(1);
        g_u.topRows(1) = A * x0 + B * u - x;
    };
    auto grad_g = [&](crvec ux, crvec v, rvec grad_u_v) {
        (void)ux;
        grad_u_v.topRows(1)    = B * v;
        grad_u_v.bottomRows(1) = -alpaqa::mat::Identity(1, 1) * v;
    };

    Problem p{2, 1, C, D, f, grad_f, g, grad_g, {}, {}, {}};

    ALMParams almparam;
    almparam.ε        = 1e-4;
    almparam.δ        = 1e-4;
    almparam.Δ        = 5; ///< Factor used in updating the penalty parameters
    almparam.Σ₀       = 1; ///< Initial penalty parameter
    almparam.ε₀       = 1e-4; ///< Initial tolerance on x
    almparam.θ        = 0.25;
    almparam.ρ        = 1e-1;
    almparam.M        = 1e9;
    almparam.Σ_max    = 1e9;
    almparam.max_iter = 10;

    PANOCParams panocparam;
    panocparam.Lipschitz.ε = 1e-6;
    panocparam.Lipschitz.δ = 1e-12;
    panocparam.max_iter    = 100;

    LBFGSParams lbfgsparam;
    lbfgsparam.memory = 10;

    ALMSolver<> solver{almparam, {panocparam, lbfgsparam}};

    vec x(2);
    x << 0.5, 0.5;
    vec y(1);
    y << 1;

    auto begin = std::chrono::high_resolution_clock::now();
    auto stats = solver(p, y, x);
    auto end   = std::chrono::high_resolution_clock::now();

    std::cout << "u = " << x.topRows(1).transpose() << std::endl;
    std::cout << "x = " << x.bottomRows(1).transpose() << std::endl;
    std::cout << "Inner: " << stats.inner.iterations
              << ", Outer: " << stats.outer_iterations << std::endl;

    auto duration =
        std::chrono::duration_cast<std::chrono::microseconds>(end - begin);
    std::cout << duration.count() << "µs" << std::endl;

    EXPECT_NEAR(x(0), -0.454545, 1e-4);
}

TEST(ALM, multipleshooting8D) {
    using namespace alpaqa;

    unsigned nx = 8, nu = 8;
    unsigned n = nu + nx;
    unsigned m = nx;

    Box C{vec(n), vec(n)};
    C.lowerbound.topRows(nu).fill(-1);
    C.lowerbound.bottomRows(nx).fill(-inf);
    C.upperbound.topRows(nu).fill(1);
    C.upperbound.bottomRows(nx).fill(inf);
    Box D{vec(m), vec(m)};
    D.lowerbound.fill(0);
    D.upperbound.fill(0);

#include "matrices/test-alm.mat.ipp"

    using Diag = Eigen::DiagonalMatrix<real_t, Eigen::Dynamic, Eigen::Dynamic>;

    auto Q = Diag(nx);
    Q.diagonal().fill(10);
    auto R = Diag(nu);
    R.diagonal().fill(1);

    auto x0 = vec(nx);
    x0.fill(1);

    auto f = [&](crvec ux) {
        auto u = ux.topRows(nu);
        auto x = ux.bottomRows(nx);
        return x.dot(Q * x) + u.dot(R * u);
    };
    auto grad_f = [&](crvec ux, rvec grad_f) {
        auto u                = ux.topRows(nu);
        auto x                = ux.bottomRows(nx);
        grad_f.topRows(nu)    = 2 * (R * u);
        grad_f.bottomRows(nx) = 2 * (Q * x);
    };
    auto g = [&](crvec ux, rvec g_u) {
        auto u         = ux.topRows(nu);
        auto x         = ux.bottomRows(nx);
        g_u.topRows(m) = A * x0 + B * u - x;
    };
    auto grad_g = [&](crvec ux, crvec v, rvec grad_u_v) {
        (void)ux;
        grad_u_v.topRows(nu)    = v.transpose() * B;
        grad_u_v.bottomRows(nx) = -alpaqa::mat::Identity(nx, nx) * v;
    };

    Problem p{n, m, C, D, f, grad_f, g, grad_g, {}, {}, {}};

    ALMParams almparam;
    almparam.ε        = 1e-4;
    almparam.δ        = 1e-4;
    almparam.Δ        = 5; ///< Factor used in updating the penalty parameters
    almparam.Σ₀       = 1; ///< Initial penalty parameter
    almparam.ε₀       = 1e-4; ///< Initial tolerance on x
    almparam.θ        = 0.25;
    almparam.ρ        = 1e-1;
    almparam.M        = 1e9;
    almparam.Σ_max    = 1e9;
    almparam.max_iter = 20;

    PANOCParams panocparam;
    panocparam.Lipschitz.ε = 1e-6;
    panocparam.Lipschitz.δ = 1e-12;
    panocparam.max_iter    = 200;

    LBFGSParams lbfgsparam;
    lbfgsparam.memory = 10;

    ALMSolver<> solver{almparam, {panocparam, lbfgsparam}};

    vec x(n);
    x.fill(5);
    vec y(m);
    y.fill(1);

    auto begin = std::chrono::high_resolution_clock::now();
    auto stats = solver(p, y, x);
    auto end   = std::chrono::high_resolution_clock::now();

    std::cout << "u = " << x.topRows(nu).transpose() << std::endl;
    std::cout << "x = " << x.bottomRows(nx).transpose() << std::endl;
    std::cout << "y = " << y.transpose() << std::endl;
    std::cout << "f(x) = " << f(x) << std::endl;
    auto gx = vec(m);
    g(x, gx);
    std::cout << "g(x) = " << gx.transpose() << std::endl;
    std::cout << "Inner: " << stats.inner.iterations
              << ", Outer: " << stats.outer_iterations << std::endl;

    auto duration =
        std::chrono::duration_cast<std::chrono::microseconds>(end - begin);
    std::cout << duration.count() << "µs" << std::endl;

    auto u  = x.topRows(nu);
    auto xs = x.bottomRows(nx);
    auto ε  = 6e-4;

    // TODO: write matcher for Eigen vectors
    EXPECT_NEAR(u(0), -1, ε);
    EXPECT_NEAR(u(1), -1, ε);
    EXPECT_NEAR(u(2), -1, ε);
    EXPECT_NEAR(u(3), -1, ε);
    EXPECT_NEAR(u(4), -1, ε);
    EXPECT_NEAR(u(5), -0.878565, ε);
    EXPECT_NEAR(u(6), -0.719274, ε);
    EXPECT_NEAR(u(7), -0.999999, ε);

    EXPECT_NEAR(xs(0), 0.312975, ε);
    EXPECT_NEAR(xs(1), -1.14893, ε);
    EXPECT_NEAR(xs(2), -0.03111, ε);
    EXPECT_NEAR(xs(3), 1.69392, ε);
    EXPECT_NEAR(xs(4), -0.11316, ε);
    EXPECT_NEAR(xs(5), 1.60049, ε);
    EXPECT_NEAR(xs(6), 0.177012, ε);
    EXPECT_NEAR(xs(7), 1.14043, ε);

    // TODO: scale tolerances when preconditioning is enabled
    // EXPECT_NEAR(y(0), 6.2595, ε);
    // EXPECT_NEAR(y(1), -22.9787, ε);
    // EXPECT_NEAR(y(2), -0.6222, ε);
    // EXPECT_NEAR(y(3), 33.8783, ε);
    // EXPECT_NEAR(y(4), -2.2632, ε);
    // EXPECT_NEAR(y(5), 32.0098, ε);
    // EXPECT_NEAR(y(6), 3.54023, ε);
    // EXPECT_NEAR(y(7), 22.8086, ε);
}
