#include <panoc-alm/alm.hpp>

#include "eigen-matchers.hpp"

auto inf = std::numeric_limits<pa::real_t>::infinity();

TEST(ALM, singleshooting1D) {
    using namespace pa;

    Box C{vec(1), vec(1)};
    C.lowerbound << -1;
    C.upperbound << 1;
    Box D{vec(0), vec(0)};

    using Diag = Eigen::DiagonalMatrix<real_t, Eigen::Dynamic, Eigen::Dynamic>;

    auto A = Eigen::MatrixXd::Constant(1, 1, 0.5);
    auto B = Eigen::MatrixXd::Constant(1, 1, 1);

    auto Q = Diag(1);
    Q.diagonal() << 10;
    auto R = Diag(1);
    R.diagonal() << 1;

    auto x0 = vec(1);
    x0 << 1;

    auto f = [&](const vec &u) {
        auto x1 = A * x0 + B * u;
        return x1.dot(Q * x1) + u.dot(R * u);
    };
    auto grad_f = [&](const vec &u, vec &grad_f) {
        grad_f = 2 * B * Q * (A * x0 + B * u) + 2 * R * u;
    };
    auto g      = [&](const vec &u, vec &g_u) { (void)u, void(g_u); };
    auto grad_g = [&](const vec &u, const vec &v, vec &grad_u_v) {
        (void)u, (void)v, grad_u_v(0) = 0;
    };

    Problem p{1, 0, C, D, f, grad_f, g, grad_g};

    ALMParams almparam;
    almparam.ε        = 1e-4;
    almparam.δ        = 1e-4;
    almparam.Δ        = 100;
    almparam.Σ₀       = 20;
    almparam.ε₀       = 1;
    almparam.θ        = 0.25;
    almparam.ρ        = 1e-1;
    almparam.M        = 1e9;
    almparam.σₘₐₓ     = 1e9;
    almparam.max_iter = 100;

    PANOCParams panocparam;
    panocparam.L           = 0;
    panocparam.γ           = 0;
    panocparam.σ           = 0;
    panocparam.Lipschitz.ε = 1e-6;
    panocparam.Lipschitz.δ = 1e-12;
    panocparam.lbfgs_mem   = 10;
    panocparam.max_iter    = 10000;

    ALMSolver solver{almparam, panocparam};

    vec x(1);
    x << 1;
    vec y(0);
    solver(p, y, x);

    EXPECT_NEAR(x(0), -0.454545, 1e-4);
}

TEST(ALM, multipleshooting1D) {
    using namespace pa;

    Box C{vec(2), vec(2)};
    C.lowerbound << -1, -inf;
    C.upperbound << 1, inf;
    Box D{vec(1), vec(1)};
    D.lowerbound << 0;
    D.upperbound << 0;

    using Diag = Eigen::DiagonalMatrix<real_t, Eigen::Dynamic, Eigen::Dynamic>;

    auto A = Eigen::MatrixXd::Constant(1, 1, 0.5);
    auto B = Eigen::MatrixXd::Constant(1, 1, 1);

    auto Q = Diag(1);
    Q.diagonal() << 10;
    auto R = Diag(1);
    R.diagonal() << 1;

    auto x0 = vec(1);
    x0 << 1;

    auto f = [&](const vec &ux) {
        auto u = ux.topRows(1);
        auto x = ux.bottomRows(1);
        return x.dot(Q * x) + u.dot(R * u);
    };
    auto grad_f = [&](const vec &ux, vec &grad_f) {
        auto u               = ux.topRows(1);
        auto x               = ux.bottomRows(1);
        grad_f.topRows(1)    = 2 * R * u;
        grad_f.bottomRows(1) = 2 * Q * x;
    };
    auto g = [&](const vec &ux, vec &g_u) {
        auto u         = ux.topRows(1);
        auto x         = ux.bottomRows(1);
        g_u.topRows(1) = A * x0 + B * u - x;
    };
    auto grad_g = [&](const vec &u, const vec &v, vec &grad_u_v) {
        grad_u_v.topRows(1)    = B * v;
        grad_u_v.bottomRows(1) = -Eigen::MatrixXd::Identity(1, 1) * v;
    };

    Problem p{2, 1, C, D, f, grad_f, g, grad_g};

    ALMParams almparam;
    almparam.ε        = 1e-4;
    almparam.δ        = 1e-4;
    almparam.Δ        = 100; ///< Factor used in updating the penalty parameters
    almparam.Σ₀       = 20;  ///< Initial penalty parameter
    almparam.ε₀       = 1;   ///< Initial tolerance on x
    almparam.θ        = 0.25;
    almparam.ρ        = 1e-1;
    almparam.M        = 1e9;
    almparam.σₘₐₓ     = 1e9;
    almparam.max_iter = 100;

    PANOCParams panocparam;
    panocparam.L           = 0;
    panocparam.γ           = 0;
    panocparam.σ           = 0;
    panocparam.Lipschitz.ε = 1e-6;
    panocparam.Lipschitz.δ = 1e-12;
    panocparam.lbfgs_mem   = 10;
    panocparam.max_iter    = 10000;

    ALMSolver solver{almparam, panocparam};

    vec x(2);
    x << 0.5, 0.5;
    vec y(1);
    y << 1;
    solver(p, y, x);

    std::cout << "u = " << x.topRows(1).transpose() << std::endl;
    std::cout << "x = " << x.bottomRows(1).transpose() << std::endl;

    EXPECT_NEAR(x(0), -0.454545, 1e-4);
}
