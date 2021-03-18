#include <panoc-alm/interop/cutest/CUTEstLoader.hpp>

#include <drivers/YAMLEncoder.hpp>

#include <atomic>
#include <chrono>
#include <csignal>
#include <cstring>
#include <fstream>
#include <iostream>
#include <sstream>

using namespace std::chrono_literals;
using pa::vec;

#define SOLVER_PANOC_LBFGS 1
#define SOLVER_PANOC_SLBFGS 2
#define SOLVER_LBFGSpp 3
#define SOLVER_LBFGSBpp 4

#if SOLVER == SOLVER_PANOC_SLBFGS
#include <panoc-alm/alm.hpp>
#include <panoc-alm/inner/panoc.hpp>
#include <panoc-alm/inner/directions/specialized-lbfgs.hpp>
using Solver = pa::ALMSolver<pa::PANOCSolver<pa::SpecializedLBFGS>>;
#elif SOLVER == SOLVER_PANOC_LBFGS
#include <panoc-alm/decl/alm.hpp>
#include <panoc-alm/inner/directions/decl/lbfgs.hpp>
#include <panoc-alm/inner/decl/panoc.hpp>
using Solver = pa::ALMSolver<>;
#elif SOLVER == SOLVER_LBFGSpp
#include <panoc-alm/alm.hpp>
#include <panoc-alm/inner/lbfgspp.hpp>
using Solver = pa::ALMSolver<pa::LBFGSSolver<>>;
#elif SOLVER == SOLVER_LBFGSBpp
#include <panoc-alm/alm.hpp>
#include <panoc-alm/inner/lbfgspp.hpp>
using Solver = pa::ALMSolver<pa::LBFGSBSolver<>>;
#endif

std::atomic<Solver *> acitve_solver{nullptr};
void signal_callback_handler(int signum) {
    if (signum == SIGINT) {
        if (auto *s = acitve_solver.load(std::memory_order_relaxed)) {
            std::atomic_signal_fence(std::memory_order_acquire);
            s->stop();
        }
    }
}

#if SOLVER == SOLVER_PANOC_SLBFGS
auto get_inner_solver() {
    pa::PANOCParams panocparams;
    panocparams.max_iter                       = 1000;
    panocparams.update_lipschitz_in_linesearch = true;
    panocparams.lbfgs_mem                      = 20;

    pa::LBFGSParams lbfgsparams;
    return Solver::InnerSolver(panocparams, lbfgsparams);
}
auto get_problem(const pa::Problem &p) { return p; }
const vec &get_y(const pa::Problem &, const vec &y) { return y; }
#elif SOLVER == SOLVER_PANOC_LBFGS
auto get_inner_solver() {
    pa::PANOCParams panocparams;
    panocparams.max_iter                       = 1000;
    panocparams.update_lipschitz_in_linesearch = true;
    panocparams.lbfgs_mem                      = 20;

    pa::LBFGSParams lbfgsparams;
    return Solver::InnerSolver(panocparams, lbfgsparams);
}
auto get_problem(const pa::Problem &p) { return p; }
const vec &get_y(const pa::Problem &, const vec &y) { return y; }
#elif SOLVER == SOLVER_LBFGSpp
auto get_inner_solver() {
    Solver::InnerSolver::Params params;
    params.max_iterations = 1000;
    params.m              = 20;
    return Solver::InnerSolver(params);
}
auto get_problem(const pa::Problem &p) { return pa::ProblemOnlyD(p); }
vec get_y(const pa::Problem &p, const vec &y) {
    vec r(p.m);
    r.topRows(p.m - p.n) = y;
    r.bottomRows(p.n).setZero();
    return r;
}
YAML::Emitter &operator<<(YAML::Emitter &out,
                          const pa::LBFGSSolver<>::Params &) {
    out << "todo";
    return out;
}
#elif SOLVER == SOLVER_LBFGSBpp
auto get_inner_solver() {
    Solver::InnerSolver::Params params;
    params.max_iterations = 1000;
    params.m              = 20;
    return Solver::InnerSolver(params);
}
auto get_problem(const pa::Problem &p) { return p; }
const vec &get_y(const pa::Problem &, const vec &y) { return y; }
YAML::Emitter &operator<<(YAML::Emitter &out,
                          const pa::LBFGSBSolver<>::Params &) {
    out << "todo";
    return out;
}
#endif

int main(int argc, char *argv[]) {
    using namespace std::string_literals;
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <problem name> <output folder>"
                  << std::endl;
        return 1;
    }

    pa::ALMParams almparams;
    almparams.max_iter        = 200;
    almparams.max_time        = 1min + 30s;
    almparams.preconditioning = false;
    // almparams.print_interval  = 1;
    // almparams.Σ₀ = 1e-2;
    // almparams.ε₀ = 1e-5;
    // almparams.Δ = 1.1;

    Solver solver{almparams, get_inner_solver()};

    if (std::strcmp(argv[1], "parameters") == 0) {
        std::ofstream f(argv[2] + "/"s + argv[1] + ".yaml");
        YAML::Emitter out(f);
        out << YAML::BeginMap;
        out << YAML::Key << "solver" << YAML::Value << solver.get_name();
        out << YAML::Key << "outer" << YAML::Value << solver.get_params();
        out << YAML::Key << "inner" << YAML::Value
            << solver.inner_solver.get_params();
        out << YAML::EndMap;
        return 0;
    }

    std::string prob_dir = "CUTEst/"s + argv[1];
    CUTEstProblem cp(prob_dir + "/libcutest-problem-" + argv[1] + ".so",
                     prob_dir + "/OUTSDIF.d");

    std::atomic_signal_fence(std::memory_order_release);
    acitve_solver.store(&solver, std::memory_order_relaxed);
    signal(SIGINT, signal_callback_handler);

    auto problem_cnt = pa::ProblemWithCounters(cp.problem);
    auto problem     = get_problem(problem_cnt);

    vec x = cp.x0;
    vec y = get_y(problem, cp.y0);

    auto status = solver(problem, y, x);
    // ??? TODO: fence
    acitve_solver.store(nullptr, std::memory_order_relaxed);
    // ??? TODO: fence
    auto report = cp.get_report();

    auto f_star = cp.problem.f(x);

    YAML::Emitter out;
    out << YAML::BeginMap;
    out << YAML::Key << "name" << YAML::Value << cp.name;
    out << YAML::Key << "n" << YAML::Value << problem.n;
    out << YAML::Key << "m" << YAML::Value << problem.m;
    out << YAML::Key << "box constraints x" << YAML::Value
        << cp.number_box_constraints;
    out << YAML::Key << "solver" << YAML::Value << solver.get_name();
    out << YAML::Key << "status" << YAML::Value << status.status;
    out << YAML::Key << "outer iterations" << YAML::Value
        << status.outer_iterations;
    out << YAML::Key << "inner iterations" << YAML::Value
        << status.inner_iterations;
    out << YAML::Key << "inner convergence failures" << YAML::Value
        << status.inner_convergence_failures;
    out << YAML::Key << "elapsed time" << YAML::Value
        << std::chrono::duration<double>(status.elapsed_time).count();
    out << YAML::Key << "ε" << YAML::Value << status.ε;
    out << YAML::Key << "δ" << YAML::Value << status.δ;
    out << YAML::Key << "f" << YAML::Value << f_star;
    out << YAML::Key << "counters" << YAML::Value << problem_cnt.evaluations;
    out << YAML::Key << "linesearch failures" << YAML::Value
        << status.inner_linesearch_failures;
    out << YAML::Key << "L-BFGS failures" << YAML::Value
        << status.inner_lbfgs_failures;
    out << YAML::Key << "L-BFGS rejected" << YAML::Value
        << status.inner_lbfgs_rejected;
    out << YAML::Key << "‖Σ‖" << YAML::Value << status.norm_penalty;
    out << YAML::Key << "‖x‖" << YAML::Value << x.norm();
    out << YAML::Key << "‖y‖" << YAML::Value << y.norm();
    // out << YAML::Key << "x" << YAML::Value << x;
    // out << YAML::Key << "y" << YAML::Value << y;
    out << YAML::EndMap;
    // out << report;

    std::ofstream(argv[2] + "/"s + argv[1] + ".yaml")
        << out.c_str() << std::endl;

    std::cout << out.c_str() << std::endl;
}
