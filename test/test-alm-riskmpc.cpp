#include <iomanip>
#include <alpaqa/decl/alm.hpp>
#include <alpaqa/inner/decl/panoc.hpp>
#include <alpaqa/inner/directions/decl/lbfgs.hpp>
#include <alpaqa/reference-problems/riskaverse-mpc.hpp>

#include "eigen-matchers.hpp"

TEST(ALM, riskaverse) {
    using namespace alpaqa;

    Problem p = problems::riskaverse_mpc_problem();
    ProblemWithCounters pc(p);

    ALMParams almparam;
    almparam.ε        = 1e-8;
    almparam.δ        = 1e-8;
    almparam.Δ        = 20;  ///< Factor used in updating the penalty parameters
    almparam.Σ₀       = 0;   ///< Initial penalty parameter
    almparam.σ₀       = 1e3; ///< Initial penalty parameter factor
    almparam.ε₀       = 1e-1; ///< Initial tolerance on x
    almparam.θ        = 0.25;
    almparam.ρ        = 1e-1;
    almparam.M        = 1e9;
    almparam.Σ_max    = 1e9;
    almparam.max_iter = 100;
    almparam.preconditioning = false;
    // almparam.single_penalty_factor = true;

    PANOCParams panocparam;
    panocparam.max_iter                       = 1000;
    panocparam.update_lipschitz_in_linesearch = true;

    panocparam.print_interval = 0;
    almparam.print_interval   = 1;

    LBFGSParams lbfgsparam;
    lbfgsparam.memory = 20;

    ALMSolver<> solver{almparam, {panocparam, lbfgsparam}};

    vec x(p.n);
    x.fill(0);
    vec λ(p.m);
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

    std::cout << std::setprecision(17);
    std::cout << "u = " << x.topRows(2).transpose() << std::endl;
    std::cout << "x = " << x.bottomRows(5).transpose() << std::endl;
    std::cout << "λ = " << λ.transpose() << std::endl;
    std::cout << "f(x) = " << p.f(x) << std::endl;
    auto gx = vec(p.m);
    p.g(x, gx);
    std::cout << "g(x) = " << gx.transpose() << std::endl;
    std::cout << "Inner: " << stats.inner.iterations
              << ", Outer: " << stats.outer_iterations << std::endl;

    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(
        (end - begin) / N);
    std::cout << duration.count() << "µs" << std::endl;

    std::cout << "# eval f:  " << pc.evaluations->f << std::endl;
    std::cout << "# eval ∇f: " << pc.evaluations->grad_f << std::endl;
    std::cout << "# eval g:  " << pc.evaluations->g << std::endl;
    std::cout << "# eval ∇g: " << pc.evaluations->grad_g_prod << std::endl;

    std::cout << "Status: " << stats.status << std::endl;

    EXPECT_NEAR(x(0), -4.87804878, 5e-5);
    EXPECT_NEAR(x(1), -4.87804878, 5e-5);
}
