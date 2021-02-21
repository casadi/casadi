#include <panoc-alm/interop/cutest/CUTEstLoader.hpp>
#include <panoc-alm/alm.hpp>

int main() {
    const char *so_fname =
        "/home/pieter/GitHub/PANOC-ALM/build/examples/CUTEst/"
        "Rosenbrock/CUTEst/ROSENBR/libcutest-ROSENBR.so";
    const char *outsdif_fname = "/home/pieter/GitHub/PANOC-ALM/build/examples/"
                                "CUTEst/Rosenbrock/CUTEst/ROSENBR/OUTSDIF.d";

    auto p = load_CUTEst_problem(so_fname, outsdif_fname);
    
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
    almparam.max_iter = 20;

    pa::PANOCParams panocparam;
    panocparam.Lipschitz.ε = 1e-6;
    panocparam.Lipschitz.δ = 1e-12;
    panocparam.lbfgs_mem   = 10;
    panocparam.max_iter    = 50;

    pa::ALMSolver solver{almparam, panocparam};
    auto stats = solver(p.problem, p.y0, p.x0);

    std::cout << "x = " << p.x0.transpose() << std::endl;
    std::cout << "y = " << p.y0.transpose() << std::endl;

    std::cout << "x_l = " << p.problem.C.lowerbound.transpose() << std::endl;
    std::cout << "x_u = " << p.problem.C.upperbound.transpose() << std::endl;

    std::cout << "inner: " << stats.inner_iterations << std::endl;
    std::cout << "outer: " << stats.outer_iterations << std::endl;
}