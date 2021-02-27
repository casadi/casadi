#include <gtest/gtest.h>
#include <panoc-alm-ref/fd.hpp>
#include <panoc-alm-ref/panoc-ref.hpp>
#include <panoc-alm/panoc.hpp>

constexpr static auto inf = std::numeric_limits<pa::real_t>::infinity();
constexpr static auto NaN = std::numeric_limits<pa::real_t>::quiet_NaN();

TEST(PANOC, DISABLED_ref) {
    using pa::Box;
    using pa::Problem;
    using pa::real_t;
    using pa::vec;

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

    auto f = [=](const vec &x, const vec &u) { return A * x + B * u; };

    using Diag = Eigen::DiagonalMatrix<real_t, Eigen::Dynamic, Eigen::Dynamic>;

    auto Q = Diag(nx);
    Q.diagonal().fill(10);
    auto R = Diag(nu);
    R.diagonal().fill(1);

    auto x0 = vec(nx);
    // x0 << 10, 20, 2, 3;
    x0.fill(10);

    real_t scale = 1;
    auto obj_f   = [=](const vec &ux) { return scale * s(ux)(0); };
    auto grad_f  = [=](const vec &ux, vec &grad_f) {
        (void)ux;
        grad_f.fill(0);
        s(grad_f)(0) = scale;
        // std::cout << "∇f = " << grad_f.transpose() << std::endl;
        // std::cout << "u  = " << u(ux).transpose() << std::endl;
        // std::cout << "s  = " << s(ux).transpose() << std::endl;
        // std::cout << "y  = " << y(ux).transpose() << std::endl;
    };
    auto g = [=](const vec &ux, vec &g_u) {
        g_u(0) = y(ux)(0) - y(ux)(1) - s(ux)(0);
        g_u(1) = f(x0, u(ux)).dot(Q * f(x0, u(ux))) - s(ux)(1);
        g_u(2) = x0.dot(Q * x0) + u(ux).dot(R * u(ux)) -
                 (y(ux)(0) - y(ux)(1) - y(ux)(2) - s(ux)(1));
    };
    auto grad_g_mat = [=](const vec &ux) {
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

        return grad;
    };
    auto grad_g = [=](const vec &ux, const vec &v, vec &grad_u_v) {
        grad_u_v = grad_g_mat(ux) * v;
    };

    Problem p{n, m, C, D, obj_f, grad_f, g, grad_g};

    pa::PANOCParams params;
    params.lbfgs_mem = 20;
    params.max_iter  = 1000;
    params.τ_min     = 1. / 16;

    pa::PANOCParams params_ref = params;
    params_ref.max_iter        = 1000;

    pa::PANOCSolver solver{params};
    pa_ref::PANOCSolver solver_ref{params_ref};

    vec λ     = vec::Ones(m);
    vec λ_ref = vec::Ones(m);

    vec x     = vec::Zero(n);
    vec x_ref = vec::Zero(n);

    vec z     = vec::Constant(m, NaN);
    vec z_ref = vec::Constant(m, NaN);

    vec err_z     = vec::Constant(m, NaN);
    vec err_z_ref = vec::Constant(m, NaN);

    vec Σ(m);
    // Σ << 1 / 5.4975667769191568e+02, 1 / 2.0648460636997861e+04,
    //     1 / 3.7043582593175870e+02;
    Σ.fill(20);

    real_t ε = 1e-4;

    auto stats     = solver(p, x, z, λ, err_z, Σ, ε);
    auto stats_ref = solver_ref(p, x_ref, z_ref, λ_ref, err_z, Σ, ε);

    EXPECT_EQ(stats.iterations, stats_ref.iterations);
    std::cout << "Normal" << std::endl;
    std::cout << "u = " << u(x).transpose() << "\ts = " << s(x).transpose()
              << "\ty = " << y(x).transpose() << std::endl;
    std::cout << stats.iterations << std::endl;
    std::cout << stats.status << std::endl << std::endl;

    std::cout << "Ref" << std::endl;
    std::cout << "u = " << u(x_ref).transpose()
              << "\ts = " << s(x_ref).transpose()
              << "\ty = " << y(x_ref).transpose() << std::endl;
    std::cout << stats_ref.iterations << std::endl;
    std::cout << stats_ref.status << std::endl << std::endl;

    {
        vec x = vec::Ones(n);
        pa::mat grad(n, m);
        for (unsigned i = 0; i < m; ++i) {
            auto fun = [=](const vec &x) {
                vec g(m), gh(m);
                p.g(x, g);
                return g(i);
            };
            grad.col(i) = pa_ref::finite_diff(fun, x);
        }
        std::cout << grad_g_mat(x) << std::endl << std::endl;
        std::cout << grad << std::endl << std::endl;
    }
}

// Check that evaluations of ψ and ∇ψ in PANOC are correct
TEST(PANOC, finiteDiffPsi) {
    using pa::Box;
    using pa::Problem;
    using pa::real_t;
    using pa::vec;

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

    auto f = [=](const vec &x, const vec &u) { return A * x + B * u; };

    using Diag = Eigen::DiagonalMatrix<real_t, Eigen::Dynamic, Eigen::Dynamic>;

    auto Q = Diag(nx);
    Q.diagonal().fill(10);
    auto R = Diag(nu);
    R.diagonal().fill(1);

    auto x0 = vec(nx);
    x0 << 10, 20, 2, 3;

    auto obj_f  = [=](const vec &ux) { return s(ux)(0); };
    auto grad_f = [=](const vec &ux, vec &grad_f) {
        (void)ux;
        grad_f.fill(0);
        s(grad_f)(0) = 1;
    };
    auto g = [=](const vec &ux, vec &g_u) {
        g_u(0) = y(ux)(0) - y(ux)(1) - s(ux)(0);
        g_u(1) = f(x0, u(ux)).dot(Q * f(x0, u(ux))) - s(ux)(1);
        g_u(2) = x0.dot(Q * x0) + u(ux).dot(R * u(ux)) -
                 (y(ux)(0) - y(ux)(1) - y(ux)(2) - s(ux)(1));
    };
    auto grad_g_mat = [=](const vec &ux) {
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

        return grad;
    };
    auto grad_g = [=](const vec &ux, const vec &v, vec &grad_u_v) {
        grad_u_v = grad_g_mat(ux) * v;
    };

    Problem p{n, m, C, D, obj_f, grad_f, g, grad_g};

    {
        vec λ = vec::Ones(m);
        vec Σ = vec::Ones(m);

        vec x = vec::Ones(n);
        vec g = pa_ref::detail::eval_g(p, x);
        vec grad_f(n);
        p.grad_f(x, grad_f);

        auto ψ = [=](const vec &x) {
            return pa_ref::detail::eval_ψ(p, x, λ, Σ);
        };

        std::cout << "x     = " << x.transpose() << std::endl;
        std::cout << "λ     = " << λ.transpose() << std::endl;
        std::cout << "Σ     = " << Σ.transpose() << std::endl;
        std::cout << "f(x)  = " << p.f(x) << std::endl;
        std::cout << "g(x)  = " << g.transpose() << std::endl;
        std::cout << "ẑ(x)  = "
                  << pa_ref::detail::eval_ẑ(p, x, λ, Σ).transpose()
                  << std::endl;
        std::cout << "ŷ(x)  = "
                  << pa_ref::detail::eval_ŷ(p, x, λ, Σ).transpose()
                  << std::endl;
        std::cout << "ψ(x)  = " << ψ(x) << std::endl;

        std::cout << "∇f(x) = " << grad_f.transpose() << std::endl;
        std::cout << "∇g(x) = \n";
        std::cout << grad_g_mat(x) << std::endl;
        std::cout << "∇ψ(x) = "
                  << pa_ref::detail::eval_grad_ψ(p, x, λ, Σ).transpose()
                  << std::endl;
        std::cout << "∇ψ(x) = " << pa_ref::finite_diff(ψ, x).transpose()
                  << " (finite differences)" << std::endl;
    }
}