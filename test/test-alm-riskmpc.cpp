#include <panoc-alm/decl/alm.hpp>
#include <panoc-alm/inner/decl/panoc.hpp>
#include <panoc-alm/inner/directions/decl/lbfgs.hpp>

#include "eigen-matchers.hpp"

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
    ProblemWithCounters pc(p);

    ALMParams almparam;
    almparam.ε        = 1e-8;
    almparam.δ        = 1e-8;
    almparam.Δ        = 20; ///< Factor used in updating the penalty parameters
    almparam.Σ₀       = 0;   ///< Initial penalty parameter
    almparam.σ₀       = 1e-2; ///< Initial penalty parameter factor
    almparam.ε₀       = 1e-1; ///< Initial tolerance on x
    almparam.θ        = 0.25;
    almparam.ρ        = 1e-1;
    almparam.M        = 1e9;
    almparam.Σₘₐₓ     = 1e9;
    almparam.max_iter = 100;
    almparam.preconditioning = false;

    PANOCParams panocparam;
    panocparam.Lipschitz.ε                    = 1e-11;
    panocparam.Lipschitz.δ                    = 1e-11;
    panocparam.lbfgs_mem                      = 20;
    panocparam.max_iter                       = 1000;
    panocparam.update_lipschitz_in_linesearch = true;

    panocparam.print_interval = 0;
    almparam.print_interval   = 1;

    LBFGSParams lbfgsparam;

    ALMSolver<> solver{almparam, {panocparam, lbfgsparam}};

    vec x(n);
    x.fill(0);
    vec λ(m);
    λ.fill(0);

    ALMSolver<>::Stats stats;

    constexpr unsigned N = 1;
    auto begin           = std::chrono::high_resolution_clock::now();
    for (volatile unsigned i = 0; i < N; ++i) {
        x.fill(0);
        λ.fill(0);
        stats = solver(pc, λ, x);
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

    std::cout << "# eval f:  " << pc.evaluations.f << std::endl;
    std::cout << "# eval ∇f: " << pc.evaluations.grad_f << std::endl;
    std::cout << "# eval g:  " << pc.evaluations.g << std::endl;
    std::cout << "# eval ∇g: " << pc.evaluations.grad_g << std::endl;

    std::cout << "Status: " << stats.status << std::endl;

    // EXPECT_NEAR(x(0), -0.454545, 1e-4);
}
