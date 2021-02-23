#include <panoc-alm/alm.hpp>

#include "eigen-matchers.hpp"

#include <chrono>

constexpr static auto inf = std::numeric_limits<pa::real_t>::infinity();

TEST(ALM, riskaverse) {
    using namespace pa;

    unsigned nu = 2;
    unsigned nx = 4;
    unsigned ns = 2;
    unsigned ny = 3;
    unsigned n  = nu + ns + ny;
    unsigned m  = 3;

    auto u = [nu](auto &&x) { return x.segment(0, nu); };
    auto s = [nu, ns](auto &&x) { return x.segment(nu, ns); };
    auto y = [nu, ns, ny](auto &&x) { return x.segment(nu + ns, ny); };

    Box C{vec(n), vec(n)};
    u(C.lowerbound).fill(-10);
    u(C.upperbound).fill(10);
    s(C.lowerbound).fill(-inf);
    s(C.upperbound).fill(inf);
    y(C.lowerbound).fill(0);
    y(C.upperbound).fill(inf);
    Box D{vec(m), vec(m)};
    D.lowerbound.fill(-inf);
    D.upperbound.fill(0);

    real_t Ts = 0.05;

    pa::mat A = pa::mat::Identity(nx, nx);
    A(0, 2)   = Ts;
    A(1, 3)   = Ts;
    pa::mat B = pa::mat::Zero(nx, nu);
    B(2, 0)   = Ts;
    B(3, 1)   = Ts;

    auto f = [&](const vec &x, const vec &u) { return A * x + B * u; };

    using Diag = Eigen::DiagonalMatrix<real_t, Eigen::Dynamic, Eigen::Dynamic>;

    auto Q = Diag(nx);
    Q.diagonal().fill(10);
    auto R = Diag(nu);
    R.diagonal().fill(1);

    auto x0 = vec(nx);
    x0.fill(10);

    auto obj_f  = [&](const vec &ux) { return s(ux)(0); };
    auto grad_f = [&](const vec &ux, vec &grad_f) {
        (void)ux;
        grad_f.fill(0);
        s(grad_f)(0) = 1;
        // std::cout << "∇f = " << grad_f.transpose() << std::endl;
        // std::cout << "u  = " << u(ux).transpose() << std::endl;
        // std::cout << "s  = " << s(ux).transpose() << std::endl;
        // std::cout << "y  = " << y(ux).transpose() << std::endl;
    };
    auto g = [&](const vec &ux, vec &g_u) {
        g_u(0) = y(ux)(0) - y(ux)(1) - s(ux)(0);
        g_u(1) = f(x0, u(ux)).dot(Q * f(x0, u(ux))) - s(ux)(1);
        g_u(2) = x0.dot(Q * x0) + u(ux).dot(R * u(ux)) -
                 (y(ux)(0) - y(ux)(1) - y(ux)(2) - s(ux)(1));
    };
    auto grad_g = [&](const vec &ux, const vec &v, vec &grad_u_v) {
        pa::mat grad      = pa::mat::Zero(n, m);
        s(grad.col(0))(0) = -1;
        y(grad.col(0))(0) = 1;
        y(grad.col(0))(1) = -1;

        u(grad.col(1))    = 2 * B.transpose() * (Q * (A * x0 + B * u(ux)));
        s(grad.col(1))(1) = -1;

        u(grad.col(2))    = 2 * (R * u(ux));
        s(grad.col(2))(1) = 1;
        y(grad.col(2))(0) = -1;
        y(grad.col(2))(1) = 1;
        y(grad.col(2))(2) = 1;

        grad_u_v = grad * v;
    };

    Problem p{n, m, C, D, obj_f, grad_f, g, grad_g};

    ALMParams almparam;
    almparam.ε        = 1e-5;
    almparam.δ        = 1e-5;
    almparam.Δ        = 5; ///< Factor used in updating the penalty parameters
    almparam.Σ₀       = 1e-2; ///< Initial penalty parameter
    almparam.ε₀       = 1e-4; ///< Initial tolerance on x
    almparam.θ        = 0.25;
    almparam.ρ        = 1e-1;
    almparam.M        = 1e9;
    almparam.σₘₐₓ     = 1e9;
    almparam.max_iter = 100;

    PANOCParams panocparam;
    panocparam.Lipschitz.ε = 1e-6;
    panocparam.Lipschitz.δ = 1e-12;
    panocparam.lbfgs_mem   = 10;
    panocparam.max_iter    = 200;

    ALMSolver solver{almparam, panocparam};

    vec x(n);
    x.fill(0);
    vec λ(m);
    λ.fill(0);

    auto stats = solver(p, λ, x);

    constexpr unsigned N = 10;
    auto begin           = std::chrono::high_resolution_clock::now();
    for (volatile unsigned i = 0; i < N; ++i) {
        x.fill(0);
        λ.fill(0);
        stats = solver(p, λ, x);
    }
    auto end = std::chrono::high_resolution_clock::now();

    std::cout << "u = " << x.topRows(nu).transpose() << std::endl;
    std::cout << "x = " << x.bottomRows(nx).transpose() << std::endl;
    std::cout << "λ = " << λ.transpose() << std::endl;
    std::cout << "f(x) = " << obj_f(x) << std::endl;
    auto gx = vec(m);
    g(x, gx);
    std::cout << "g(x) = " << gx.transpose() << std::endl;
    std::cout << "Inner: " << stats.inner_iterations
              << ", Outer: " << stats.outer_iterations << std::endl;

    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(
        (end - begin) / N);
    std::cout << duration.count() << "µs" << std::endl;

    // EXPECT_NEAR(x(0), -0.454545, 1e-4);
}
