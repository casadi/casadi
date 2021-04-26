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
#define SOLVER_PGA 5
#define SOLVER_GAAPGA 6
#define SOLVER_PANOC_ANDERSON 7
#define SOLVER_PANOC_2ND 8
#define SOLVER_PANOC_2ND_LBFGS 9

#if SOLVER == SOLVER_PANOC_SLBFGS
#include <panoc-alm/alm.hpp>
#include <panoc-alm/inner/directions/specialized-lbfgs.hpp>
#include <panoc-alm/inner/panoc.hpp>
using Solver = pa::ALMSolver<pa::PANOCSolver<pa::SpecializedLBFGS>>;
#elif SOLVER == SOLVER_PANOC_LBFGS
#include <panoc-alm/decl/alm.hpp>
#include <panoc-alm/inner/decl/panoc.hpp>
#include <panoc-alm/inner/directions/decl/lbfgs.hpp>
using Solver = pa::ALMSolver<>;
#elif SOLVER == SOLVER_PANOC_2ND
#include <panoc-alm/alm.hpp>
#include <panoc-alm/inner/second-order-panoc.hpp>
using Solver = pa::ALMSolver<pa::SecondOrderPANOCSolver>;
#elif SOLVER == SOLVER_PANOC_2ND_LBFGS
#include <panoc-alm/alm.hpp>
#include <panoc-alm/inner/second-order-panoc-lbfgs.hpp>
using Solver = pa::ALMSolver<pa::SecondOrderPANOCLBFGSSolver>;
#elif SOLVER == SOLVER_LBFGSpp
#include <panoc-alm/alm.hpp>
#include <panoc-alm/inner/lbfgspp.hpp>
using Solver = pa::ALMSolver<pa::LBFGSSolver<>>;
#elif SOLVER == SOLVER_LBFGSBpp
#include <panoc-alm/alm.hpp>
#include <panoc-alm/inner/lbfgspp.hpp>
using Solver = pa::ALMSolver<pa::LBFGSBSolver<>>;
#elif SOLVER == SOLVER_PGA
#include <panoc-alm/alm.hpp>
#include <panoc-alm/inner/pga.hpp>
using Solver = pa::ALMSolver<pa::PGA>;
#elif SOLVER == SOLVER_GAAPGA
#include <panoc-alm/alm.hpp>
#include <panoc-alm/inner/guarded-aa-pga.hpp>
using Solver = pa::ALMSolver<pa::GuardedAAPGA>;
#elif SOLVER == SOLVER_PANOC_ANDERSON
#include <panoc-alm/alm.hpp>
#include <panoc-alm/inner/directions/anderson-acceleration.hpp>
#include <panoc-alm/inner/panoc.hpp>
using Solver = pa::ALMSolver<pa::PANOCSolver<pa::AndersonAccel>>;
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
    panocparams.lbfgs_stepsize                 = panocparams.BasedOnCurvature;

    pa::LBFGSParams lbfgsparams;
    return Solver::InnerSolver(panocparams, lbfgsparams);
}
auto get_problem(const pa::Problem &p) { return p; }
const vec &get_y(const pa::Problem &, const vec &y) { return y; }
#elif SOLVER == SOLVER_PANOC_2ND
auto get_inner_solver() {
    pa::SecondOrderPANOCParams panocparams;
    panocparams.max_iter                       = 1000;
    panocparams.update_lipschitz_in_linesearch = true;
    panocparams.max_time                       = 30s;

    return Solver::InnerSolver(panocparams);
}
auto get_problem(const pa::Problem &p) { return p; }
const vec &get_y(const pa::Problem &, const vec &y) { return y; }
inline YAML::Emitter &operator<<(YAML::Emitter &out,
                                 const pa::SecondOrderPANOCParams &p) {
    out << YAML::BeginMap;
    out << YAML::Key << "Lipschitz" << YAML::Value << YAML::BeginMap;
    out << YAML::Key << "ε" << YAML::Value << p.Lipschitz.ε;
    out << YAML::Key << "δ" << YAML::Value << p.Lipschitz.δ;
    out << YAML::Key << "Lγ_factor" << YAML::Value << p.Lipschitz.Lγ_factor;
    out << YAML::EndMap;
    out << YAML::Key << "max_iter" << YAML::Value << p.max_iter;
    out << YAML::Key << "max_time" << YAML::Value << p.max_time.count();
    out << YAML::Key << "τ_min" << YAML::Value << p.τ_min;
    out << YAML::Key << "update_lipschitz_in_linesearch" << YAML::Value
        << p.update_lipschitz_in_linesearch;
    out << YAML::Key << "alternative_linesearch_cond" << YAML::Value
        << p.alternative_linesearch_cond;
    out << YAML::EndMap;
    return out;
}
#elif SOLVER == SOLVER_PANOC_2ND_LBFGS
auto get_inner_solver() {
    pa::SecondOrderPANOCLBFGSParams panocparams;
    panocparams.max_iter                       = 1000;
    panocparams.update_lipschitz_in_linesearch = true;
    panocparams.lbfgs_mem                      = 20;
    panocparams.lbfgs_stepsize = pa::LBFGSStepSize::BasedOnCurvature;
    pa::LBFGSParams lbfgsparams;

    return Solver::InnerSolver(panocparams, lbfgsparams);
}
auto get_problem(const pa::Problem &p) { return p; }
const vec &get_y(const pa::Problem &, const vec &y) { return y; }
inline YAML::Emitter &operator<<(YAML::Emitter &out,
                                 const pa::SecondOrderPANOCLBFGSParams &p) {
    out << YAML::BeginMap;
    out << YAML::Key << "Lipschitz" << YAML::Value << YAML::BeginMap;
    out << YAML::Key << "ε" << YAML::Value << p.Lipschitz.ε;
    out << YAML::Key << "δ" << YAML::Value << p.Lipschitz.δ;
    out << YAML::Key << "Lγ_factor" << YAML::Value << p.Lipschitz.Lγ_factor;
    out << YAML::EndMap;
    out << YAML::Key << "max_iter" << YAML::Value << p.max_iter;
    out << YAML::Key << "max_time" << YAML::Value << p.max_time.count();
    out << YAML::Key << "τ_min" << YAML::Value << p.τ_min;
    out << YAML::Key << "update_lipschitz_in_linesearch" << YAML::Value
        << p.update_lipschitz_in_linesearch;
    out << YAML::Key << "alternative_linesearch_cond" << YAML::Value
        << p.alternative_linesearch_cond;
    out << YAML::Key << "hessian_vec_finited_differences" << YAML::Value
        << p.hessian_vec_finited_differences;
    out << YAML::Key << "full_augmented_hessian" << YAML::Value
        << p.full_augmented_hessian;
    out << YAML::EndMap;
    return out;
}
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
#elif SOLVER == SOLVER_PGA
auto get_inner_solver() {
    pa::PGAParams params;
    params.max_iter = 1000;

    return Solver::InnerSolver(params);
}
auto get_problem(const pa::Problem &p) { return p; }
const vec &get_y(const pa::Problem &, const vec &y) { return y; }
inline YAML::Emitter &operator<<(YAML::Emitter &out, const pa::PGAParams &p) {
    out << YAML::BeginMap;
    out << YAML::Key << "Lipschitz" << YAML::Value << YAML::BeginMap;
    out << YAML::Key << "ε" << YAML::Value << p.Lipschitz.ε;
    out << YAML::Key << "δ" << YAML::Value << p.Lipschitz.δ;
    out << YAML::Key << "Lγ_factor" << YAML::Value << p.Lipschitz.Lγ_factor;
    out << YAML::EndMap;
    out << YAML::Key << "max_iter" << YAML::Value << p.max_iter;
    out << YAML::Key << "max_time" << YAML::Value << p.max_time.count();
    out << YAML::EndMap;
    return out;
}
#elif SOLVER == SOLVER_GAAPGA
auto get_inner_solver() {
    pa::GuardedAAPGAParams params;
    params.max_iter               = 1000;
    params.limitedqr_mem          = 20;
    params.full_flush_on_γ_change = false;

    return Solver::InnerSolver(params);
}
auto get_problem(const pa::Problem &p) { return p; }
const vec &get_y(const pa::Problem &, const vec &y) { return y; }
inline YAML::Emitter &operator<<(YAML::Emitter &out,
                                 const pa::GuardedAAPGAParams &p) {
    out << YAML::BeginMap;
    out << YAML::Key << "Lipschitz" << YAML::Value << YAML::BeginMap;
    out << YAML::Key << "ε" << YAML::Value << p.Lipschitz.ε;
    out << YAML::Key << "δ" << YAML::Value << p.Lipschitz.δ;
    out << YAML::Key << "Lγ_factor" << YAML::Value << p.Lipschitz.Lγ_factor;
    out << YAML::EndMap;
    out << YAML::Key << "limitedqr_mem" << YAML::Value << p.limitedqr_mem;
    out << YAML::Key << "max_iter" << YAML::Value << p.max_iter;
    out << YAML::Key << "max_time" << YAML::Value << p.max_time.count();
    out << YAML::EndMap;
    return out;
}
#elif SOLVER == SOLVER_PANOC_ANDERSON
auto get_inner_solver() {
    pa::PANOCParams params;
    params.max_iter                       = 1000;
    params.update_lipschitz_in_linesearch = true;
    params.lbfgs_mem                      = 20;

    return Solver::InnerSolver(params, {});
}
auto get_problem(const pa::Problem &p) { return p; }
const vec &get_y(const pa::Problem &, const vec &y) { return y; }
#endif

int main(int argc, char *argv[]) {
    using namespace std::string_literals;
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <problem name> <output folder>"
                  << std::endl;
        return 1;
    }

    pa::ALMParams almparams;
    almparams.max_iter = 240;
    almparams.max_time = 30s;
    // almparams.max_time        = 1min + 30s;
    almparams.preconditioning = false;
    // almparams.print_interval  = 1;
    almparams.Σ₀                      = 1;
    almparams.Δ                       = 10;
    almparams.max_num_initial_retries = 20;
    almparams.max_num_retries         = 20;
    almparams.max_total_num_retries   = 40;
    almparams.Δ_lower                 = 0.8;
    almparams.ρ_increase              = 2;
    // almparams.σ₀ = 2e-1;
    // almparams.ε₀ = 1e-5;

    Solver solver{almparams, get_inner_solver()};

    if (std::strcmp(argv[1], "parameters") == 0) {
        std::ofstream f(argv[2] + "/"s + argv[1] + ".yaml");
        YAML::Emitter out(f);
        out << YAML::BeginMap;
        out << YAML::Key << "solver" << YAML::Value << solver.get_name();
        out << YAML::Key << "outer" << YAML::Value << solver.get_params();
        out << YAML::Key << "inner" << YAML::Value
            << solver.inner_solver.get_params();
        out << YAML::Key << "directions" << YAML::Value
            << solver.inner_solver.lbfgs.get_params();
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
    out << YAML::Key << "inner convergence failures" << YAML::Value
        << status.inner_convergence_failures;
    out << YAML::Key << "initial penalty reduced" << YAML::Value
        << status.initial_penalty_reduced;
    out << YAML::Key << "penalty reduced" << YAML::Value
        << status.penalty_reduced;
    out << YAML::Key << "elapsed time" << YAML::Value
        << std::chrono::duration<double>(status.elapsed_time).count();
    out << YAML::Key << "ε" << YAML::Value << status.ε;
    out << YAML::Key << "δ" << YAML::Value << status.δ;
    out << YAML::Key << "inner" << YAML::Value << status.inner;
    out << YAML::Key << "‖Σ‖" << YAML::Value << status.norm_penalty;
    out << YAML::Key << "‖x‖" << YAML::Value << x.norm();
    out << YAML::Key << "‖y‖" << YAML::Value << y.norm();
    out << YAML::Key << "f" << YAML::Value << f_star;
    out << YAML::Key << "counters" << YAML::Value << problem_cnt.evaluations;
    // out << YAML::Key << "x" << YAML::Value << x;
    // out << YAML::Key << "y" << YAML::Value << y;
    out << YAML::EndMap;
    // out << report;

    std::ofstream(argv[2] + "/"s + argv[1] + ".yaml")
        << out.c_str() << std::endl;

    std::cout << out.c_str() << std::endl;
}
