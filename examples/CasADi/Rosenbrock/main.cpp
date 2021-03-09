#include <panoc-alm/decl/alm.hpp>
#include <panoc-alm/inner/decl/panoc.hpp>
#include <panoc-alm/inner/lbfgs.hpp>

#include <panoc-alm/interop/casadi/CasADiLoader.hpp>

int main(int argc, char *argv[]) {
    const char *so_name =
        "examples/CasADi/Rosenbrock/librosenbrock_functions.so";
    if (argc > 1)
        so_name = argv[1];

    constexpr auto inf = std::numeric_limits<double>::infinity();
    using pa::real_t;
    using pa::vec;

    pa::Problem p;
    p.n = 3;
    p.m = 1;
    p.C = pa::Box{
        vec::Constant(3, inf),
        vec::Constant(3, -inf),
    };
    p.D = pa::Box{
        vec::Constant(1, 0.),
        vec::Constant(1, 0.),
    };
    p.f      = load_CasADi_objective(so_name);
    p.grad_f = load_CasADi_gradient_objective(so_name);
    p.g      = load_CasADi_constraints(so_name);
    p.grad_g = load_CasADi_gradient_constraints(so_name);

    pa::ALMParams almparam;
    almparam.ε        = 1e-8;
    almparam.δ        = 1e-8;
    almparam.Δ        = 10;
    almparam.Σ₀       = 2;
    almparam.ε₀       = 1;
    almparam.θ        = 0.25;
    almparam.ρ        = 1e-1;
    almparam.M        = 1e9;
    almparam.σₘₐₓ     = 1e9;
    almparam.max_iter = 100;

    pa::PANOCParams panocparam;
    panocparam.Lipschitz.ε = 1e-6;
    panocparam.Lipschitz.δ = 1e-12;
    panocparam.lbfgs_mem   = 10;
    panocparam.max_iter    = 500;

    pa::ALMSolver<> solver{almparam, panocparam};

    vec x(3);
    x << 2.5, 3.0, 0.75;
    vec y(1);
    y << 1;

    auto stats = solver(p, y, x);

    std::cout << "x = " << x.transpose() << std::endl;
    std::cout << "y = " << y.transpose() << std::endl;
    vec g(p.m);
    p.g(x, g);
    std::cout << "g = " << g.transpose() << std::endl;
    std::cout << "f = " << p.f(x) << std::endl;

    std::cout << "inner: " << stats.inner_iterations << std::endl;
    std::cout << "outer: " << stats.outer_iterations << std::endl;
}