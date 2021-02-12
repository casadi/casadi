#include <casadi/casadi.hpp>
#include <panoc-alm/alm.hpp>

namespace cs = casadi;

int main(int argc, char *argv[]) {
    const char *so_name =
        "./examples/CasADi/Rosenbrock/librosenbrock_functions.so";
    if (argc > 1)
        so_name = argv[1];
    cs::Function f      = cs::external("f", so_name);
    cs::Function grad_f = cs::external("grad_f", so_name);
    cs::Function g      = cs::external("g", so_name);
    cs::Function grad_g = cs::external("grad_g_v", so_name);

    constexpr auto inf = std::numeric_limits<double>::infinity();
    using pa::real_t;
    using pa::vec;

    pa::Problem p;
    p.n = 3;
    p.m = 1;
    p.C = pa::Box{vec::Constant(3, inf), vec::Constant(3, -inf)};
    p.D = pa::Box{vec::Constant(1, 0.), vec::Constant(1, 0.)};
    p.f = [f{std::move(f)}](const vec &x) -> real_t {
        double f_x;
        std::vector<const double *> in{x.data()};
        std::vector<double *> out{&f_x};
        // TODO: optimize allocations (including work allocations)
        f(in, out);
        return f_x;
    };
    p.grad_f = [grad_f{std::move(grad_f)}](const vec &x, vec &grad_fx) {
        std::vector<const double *> in{x.data()};
        std::vector<double *> out{grad_fx.data()};
        // TODO: optimize allocations (including work allocations)
        grad_f(in, out);
    };
    p.g = [g{std::move(g)}](const vec &x, vec &gx) {
        std::vector<const double *> in{x.data()};
        std::vector<double *> out{gx.data()};
        // TODO: optimize allocations (including work allocations)
        g(in, out);
    };
    p.grad_g = [grad_g{std::move(grad_g)}](const vec &x, const vec &y,
                                           vec &grad_gxy) {
        std::vector<const double *> in{x.data(), y.data()};
        std::vector<double *> out{grad_gxy.data()};
        // TODO: optimize allocations (including work allocations)
        grad_g(in, out);
    };

    pa::ALMParams almparam;
    almparam.ε        = 1e-8;
    almparam.δ        = 1e-8;
    almparam.Δ        = 100;
    almparam.Σ₀       = 1e-2;
    almparam.ε₀       = 1;
    almparam.θ        = 0.25;
    almparam.ρ        = 1e-1;
    almparam.M        = 1e9;
    almparam.σₘₐₓ     = 1e9;
    almparam.max_iter = 10;

    pa::PANOCParams panocparam;
    panocparam.Lipschitz.ε = 1e-6;
    panocparam.Lipschitz.δ = 1e-12;
    panocparam.lbfgs_mem   = 10;
    panocparam.max_iter    = 50;

    pa::ALMSolver solver{almparam, panocparam};

    vec x(3);
    x << 2.5, 3.0, 0.75;
    vec y(1);
    y << 1;

    auto stats = solver(p, y, x);

    std::cout << "x = " << x.transpose() << std::endl;
    std::cout << "y = " << y.transpose() << std::endl;

    std::cout << "inner: " << stats.inner_iterations << std::endl;
    std::cout << "outer: " << stats.outer_iterations << std::endl;
}