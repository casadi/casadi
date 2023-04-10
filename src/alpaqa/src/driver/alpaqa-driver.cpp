#include <alpaqa/config/config.hpp>
#include <alpaqa/params/params.hpp>
#include <alpaqa/problem/type-erased-problem.hpp>
#include <alpaqa/util/demangled-typename.hpp>
#include <alpaqa/util/print.hpp>

#include <alpaqa/implementation/outer/alm.tpp>
#include <alpaqa/newton-tr-pantr-alm.hpp>
#include <alpaqa/panoc-alm.hpp>
#include <alpaqa/panoc-anderson-alm.hpp>
#include <alpaqa/structured-panoc-alm.hpp>
#include <alpaqa/structured-zerofpr-alm.hpp>
#include <alpaqa/zerofpr-alm.hpp>
#include <alpaqa/zerofpr-anderson-alm.hpp>

#include "output.hpp"
#include "results.hpp"
#include <atomic>
#include <type_traits>

#if WITH_IPOPT
#include <alpaqa/ipopt/ipopt-adapter.hpp>
#include <alpaqa/ipopt/ipopt-enums.hpp>
#include <IpIpoptApplication.hpp>
#endif

#if WITH_LBFGSPP
#include <alpaqa/lbfgsb-adapter.hpp>
#endif

#if WITH_LBFGSB
#include <alpaqa/lbfgsb/lbfgsb-adapter.hpp>
#endif

#include <algorithm>
#include <chrono>
#include <csignal>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iostream>
#include <iterator>
#include <mutex>
#include <random>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>
namespace fs = std::filesystem;

USING_ALPAQA_CONFIG(alpaqa::DefaultConfig);

template <class InnerSolver>
auto make_inner_solver(std::span<const std::string_view> extra_opts) {
    // Settings for the solver
    typename InnerSolver::Params solver_param;
    solver_param.max_iter       = 50'000;
    solver_param.print_interval = 0;
    solver_param.stop_crit      = alpaqa::PANOCStopCrit::ProjGradUnitNorm;
    alpaqa::params::set_params(solver_param, "solver", extra_opts);

    if constexpr (requires { typename InnerSolver::Direction; }) {
        // Settings for the direction provider
        using Direction = typename InnerSolver::Direction;
        typename Direction::DirectionParams dir_param;
        typename Direction::AcceleratorParams accel_param;
        alpaqa::params::set_params(dir_param, "dir", extra_opts);
        alpaqa::params::set_params(accel_param, "accel", extra_opts);
        return InnerSolver{
            solver_param,
            Direction{{
                .accelerator = accel_param,
                .direction   = dir_param,
            }},
        };
    } else {
        return InnerSolver{solver_param};
    }
}

#if WITH_LBFGSPP
using InnerLBFGSppSolver = alpaqa::lbfgspp::LBFGSBSolver<alpaqa::DefaultConfig>;
template <>
auto make_inner_solver<InnerLBFGSppSolver>(
    std::span<const std::string_view> extra_opts) {
    // Settings for the solver
    InnerLBFGSppSolver::Params solver_param;
    solver_param.max_iterations = 50'000;
    alpaqa::params::set_params(solver_param, "solver", extra_opts);
    return InnerLBFGSppSolver{solver_param};
}
#endif

#if WITH_LBFGSB
using InnerLBFGSBSolver = alpaqa::lbfgsb::LBFGSBSolver;
template <>
auto make_inner_solver<InnerLBFGSBSolver>(
    std::span<const std::string_view> extra_opts) {
    // Settings for the solver
    InnerLBFGSBSolver::Params solver_param;
    solver_param.max_iter       = 50'000;
    solver_param.print_interval = 0;
    solver_param.stop_crit      = alpaqa::PANOCStopCrit::ProjGradUnitNorm;
    alpaqa::params::set_params(solver_param, "solver", extra_opts);
    return InnerLBFGSBSolver{solver_param};
}
#endif

template <class InnerSolver>
auto make_solver(const auto &extra_opts) {
    // Settings for the ALM solver
    using ALMSolver = alpaqa::ALMSolver<InnerSolver>;
    typename ALMSolver::Params alm_param;
    alm_param.max_iter        = 200;
    alm_param.tolerance       = 1e-8;
    alm_param.dual_tolerance  = 1e-8;
    alm_param.print_interval  = 1;
    alm_param.print_precision = 1;
    alpaqa::params::set_params(alm_param, "alm", extra_opts);
    return ALMSolver{alm_param, make_inner_solver<InnerSolver>(extra_opts)};
}

template <class Solver>
BenchmarkResults
do_experiment_impl(LoadedProblem &problem, Solver &solver, std::ostream &os,
                   std::span<const std::string_view> extra_opts) {
    auto evals = problem.evaluations;

    // Initial guess
    vec x = vec::Zero(problem.problem.get_n()),
        y = vec::Zero(problem.problem.get_m());

    // Solve the problem
    auto stats = solver(problem.problem, x, y);
    evals      = std::make_shared<alpaqa::EvalCounter>(*evals);

    // Solve the problems again to average runtimes
    auto avg_duration = stats.elapsed_time;
    unsigned N_exp    = 0;
    alpaqa::params::set_params(N_exp, "num_exp", extra_opts);
    os.setstate(std::ios_base::badbit);
    for (unsigned i = 0; i < N_exp; ++i) {
        x.setZero();
        y.setZero();
        avg_duration += solver(problem.problem, x, y).elapsed_time;
    }
    os.clear();
    avg_duration /= (N_exp + 1);

    // Results
    auto now_ms    = timestamp_ms<std::chrono::system_clock>();
    auto f         = problem.problem.eval_f(x);
    auto kkt_err   = compute_kkt_error(problem.problem, x, y);
    real_t final_γ = 0;
    if constexpr (requires { stats.inner.final_γ; })
        final_γ = stats.inner.final_γ;
    decltype(BenchmarkResults::extra) extra{};
    if constexpr (requires { stats.inner.linesearch_failures; })
        extra.emplace_back(
            "linesearch_failures",
            static_cast<index_t>(stats.inner.linesearch_failures));
    if constexpr (requires { stats.inner.linesearch_backtracks; })
        extra.emplace_back(
            "linesearch_backtracks",
            static_cast<index_t>(stats.inner.linesearch_backtracks));
    if constexpr (requires { stats.inner.stepsize_backtracks; })
        extra.emplace_back(
            "stepsize_backtracks",
            static_cast<index_t>(stats.inner.stepsize_backtracks));
    if constexpr (requires { stats.inner.lbfgs_failures; })
        extra.emplace_back("lbfgs_failures",
                           static_cast<index_t>(stats.inner.lbfgs_failures));
    if constexpr (requires { stats.inner.lbfgs_rejected; })
        extra.emplace_back("lbfgs_rejected",
                           static_cast<index_t>(stats.inner.lbfgs_rejected));
    return BenchmarkResults{
        .problem          = problem,
        .evals            = *evals,
        .duration         = avg_duration,
        .status           = enum_name(stats.status),
        .success          = stats.status == alpaqa::SolverStatus::Converged,
        .solver           = solver.get_name(),
        .f                = f,
        .δ                = stats.δ,
        .ε                = stats.ε,
        .γ                = final_γ,
        .Σ                = stats.norm_penalty,
        .stationarity     = kkt_err.stationarity,
        .constr_violation = kkt_err.constr_violation,
        .complementarity  = kkt_err.complementarity,
        .solution         = x,
        .multipliers      = y,
        .outer_iter       = static_cast<index_t>(stats.outer_iterations),
        .inner_iter       = static_cast<index_t>(stats.inner.iterations),
        .extra            = std::move(extra),
        .options          = extra_opts,
        .timestamp        = now_ms.count(),
    };
}

#if WITH_IPOPT

auto make_ipopt_solver(std::span<const std::string_view> extra_opts) {
    using namespace Ipopt;

    // We are using the factory, since this allows us to compile this
    // example with an Ipopt Windows DLL
    SmartPtr<IpoptApplication> app = IpoptApplicationFactory();
    app->RethrowNonIpoptException(true);

    app->Options()->SetNumericValue("tol", 1e-8);
    app->Options()->SetNumericValue("constr_viol_tol", 1e-8);
    app->Options()->SetStringValue("linear_solver", "mumps");
    // app->Options()->SetStringValue("print_timing_statistics", "yes");
    // app->Options()->SetStringValue("timing_statistics", "yes");
    app->Options()->SetStringValue("hessian_approximation", "exact");

    alpaqa::params::set_params(*app, "solver", extra_opts);

    // Initialize the IpoptApplication and process the options
    ApplicationReturnStatus status = app->Initialize();
    if (status != Solve_Succeeded)
        throw std::runtime_error("Error during Ipopt initialization: " +
                                 std::string(enum_name(status)));

    return app;
}

BenchmarkResults do_experiment_impl(
    LoadedProblem &problem, Ipopt::SmartPtr<Ipopt::IpoptApplication> &solver,
    std::ostream &os, std::span<const std::string_view> extra_opts) {
    // Ipopt problem adapter
    using Problem                    = alpaqa::IpoptAdapter;
    Ipopt::SmartPtr<Ipopt::TNLP> nlp = new Problem(problem.problem);
    auto *my_nlp                     = dynamic_cast<Problem *>(GetRawPtr(nlp));

    // Initial guess
    my_nlp->initial_guess = vec::Zero(problem.problem.get_n());

    // Solve the problem
    auto t0     = std::chrono::steady_clock::now();
    auto status = solver->OptimizeTNLP(nlp);
    auto t1     = std::chrono::steady_clock::now();
    auto evals  = *problem.evaluations;

    // Solve the problems again to average runtimes
    using ns          = std::chrono::nanoseconds;
    auto avg_duration = duration_cast<ns>(t1 - t0);
    unsigned N_exp    = 0;
    alpaqa::params::set_params(N_exp, "num_exp", extra_opts);
    os.setstate(std::ios_base::badbit);
    for (unsigned i = 0; i < N_exp; ++i) {
        auto t0 = std::chrono::steady_clock::now();
        solver->OptimizeTNLP(nlp);
        auto t1 = std::chrono::steady_clock::now();
        avg_duration += duration_cast<ns>(t1 - t0);
    }
    os.clear();
    avg_duration /= (N_exp + 1);

    // Results
    auto &nlp_res = my_nlp->results;
    auto kkt_err  = compute_kkt_error(problem.problem, nlp_res.solution_x,
                                      nlp_res.solution_y);
    auto now_ms   = timestamp_ms<std::chrono::system_clock>();
    BenchmarkResults results{
        .problem          = problem,
        .evals            = evals,
        .duration         = avg_duration,
        .status           = enum_name(status),
        .success          = status == Ipopt::Solve_Succeeded,
        .solver           = "Ipopt",
        .f                = nlp_res.solution_f,
        .δ                = nlp_res.infeasibility,
        .ε                = nlp_res.nlp_error,
        .γ                = 0,
        .Σ                = 0,
        .stationarity     = kkt_err.stationarity,
        .constr_violation = kkt_err.constr_violation,
        .complementarity  = kkt_err.complementarity,
        .solution         = nlp_res.solution_x,
        .multipliers      = nlp_res.solution_y,
        .outer_iter       = nlp_res.iter_count,
        .inner_iter       = nlp_res.iter_count,
        .extra            = {},
        .options          = extra_opts,
        .timestamp        = now_ms.count(),
    };
    return results;
}

#endif

template <auto &solver_to_stop>
auto attach_cancellation(auto &solver) {
    if constexpr (requires { solver.stop(); }) {
        auto *old = solver_to_stop.exchange(&solver, std::memory_order_release);
        if (old)
            throw std::runtime_error(
                "alpaqa-driver:do_experiment is not reentrant");
        struct sigaction action;
        action.sa_handler = [](int) {
            if (auto *s = solver_to_stop.load(std::memory_order::acquire))
                s->stop();
        };
        sigemptyset(&action.sa_mask);
        action.sa_flags = 0;
        sigaction(SIGINT, &action, nullptr);
        sigaction(SIGTERM, &action, nullptr);
    }
    using solver_to_stop_t = std::remove_reference_t<decltype(solver_to_stop)>;
    auto detach_solver     = [](solver_to_stop_t *p) {
        p->store(nullptr, std::memory_order_relaxed);
    };
    return std::unique_ptr<solver_to_stop_t, decltype(detach_solver)>{
        &solver_to_stop};
}

void do_experiment(LoadedProblem &problem, auto &&solver, std::ostream &os,
                   std::span<const std::string_view> extra_opts) {
    // Try stopping gracefully when SIGINT (Ctrl+C) or SIGTERM is received.
    using stop_solver_t = std::atomic<decltype(&solver)>;
    static stop_solver_t solver_to_stop{nullptr};
    auto cancellation = attach_cancellation<solver_to_stop>(solver);
    if constexpr (requires { solver.os; })
        solver.os = &os;

    // Run experiment
    BenchmarkResults results =
        do_experiment_impl(problem, solver, os, extra_opts);

    // Print results to output
    print_results(os, results);

    // Write results to file
    auto timestamp_str = std::to_string(results.timestamp);
    auto rnd_str       = random_hex_string(std::random_device());
    auto suffix        = timestamp_str + '_' + rnd_str;
    auto results_name  = "results_" + suffix;
    os << "results: " << suffix << std::endl;
    std::ofstream res_file{results_name + ".py"};
    write_results(res_file, results);
}

std::string format_string_list(
    const auto &container,
    const auto &proj = [](const auto &x) -> decltype(auto) { return x; }) {
    if (container.empty())
        return std::string{};
    auto penult       = std::prev(container.end());
    auto quote_concat = [&](std::string &&a, const auto &b) {
        return a + "'" + std::string(proj(b)) + "', ";
    };
    return std::accumulate(container.begin(), penult, std::string{},
                           quote_concat) +
           "'" + std::string(proj(*penult)) + "'";
}

template <template <class> class Solver>
void do_experiment_with_direction(
    std::string_view direction, LoadedProblem &problem, std::ostream &os,
    std::span<const std::string_view> extra_opts) {
    if (direction.empty())
        direction = "lbfgs";
    // Available solvers
    std::map<std::string_view, std::function<void()>> directions{
        {"lbfgs",
         [&] {
             using InnerSolver = Solver<alpaqa::LBFGSDirection<config_t>>;
             do_experiment(problem, make_solver<InnerSolver>(extra_opts), os,
                           extra_opts);
         }},
        {"struclbfgs",
         [&] {
             using InnerSolver =
                 Solver<alpaqa::StructuredLBFGSDirection<config_t>>;
             do_experiment(problem, make_solver<InnerSolver>(extra_opts), os,
                           extra_opts);
         }},
        {"anderson",
         [&] {
             using InnerSolver = Solver<alpaqa::AndersonDirection<config_t>>;
             do_experiment(problem, make_solver<InnerSolver>(extra_opts), os,
                           extra_opts);
         }},
    };
    // Run experiment
    auto dir_it = directions.find(direction);
    if (dir_it != directions.end())
        dir_it->second();
    else
        throw std::invalid_argument(
            "Unknown direction '" + std::string(direction) + "'\n" +
            "  Available directions: " +
            format_string_list(directions,
                               [](const auto &x) { return x.first; }));
}

template <>
void do_experiment_with_direction<alpaqa::PANTRSolver>(
    std::string_view direction, LoadedProblem &problem, std::ostream &os,
    std::span<const std::string_view> extra_opts) {
    if (direction.empty())
        direction = "newtontr";
    // Available solvers
    std::map<std::string_view, std::function<void()>> directions{
        {"newtontr",
         [&] {
             using InnerSolver =
                 alpaqa::PANTRSolver<alpaqa::NewtonTRDirection<config_t>>;
             do_experiment(problem, make_solver<InnerSolver>(extra_opts), os,
                           extra_opts);
         }},
    };
    // Run experiment
    auto dir_it = directions.find(direction);
    if (dir_it != directions.end())
        dir_it->second();
    else
        throw std::invalid_argument(
            "Unknown direction '" + std::string(direction) + "'\n" +
            "  Available directions: " +
            format_string_list(directions,
                               [](const auto &x) { return x.first; }));
}

void print_usage(const char *a0) {
    std::cout << "Usage: " << a0
              << " [<path>/]<name> [method=<solver>] "
                 "[<key>=<value>...]\n";
}

/// Split the string @p s on the first occurrence of @p tok.
/// Returns ("", s) if tok was not found.
auto split_once(std::string_view s, char tok = '.') {
    auto tok_pos = s.find(tok);
    if (tok_pos == s.npos)
        return std::make_tuple(std::string_view{}, s);
    std::string_view key{s.begin(), s.begin() + tok_pos};
    std::string_view rem{s.begin() + tok_pos + 1, s.end()};
    return std::make_tuple(key, rem);
}

int main(int argc, char *argv[]) try {
    // Find the problem to load
    if (argc < 1)
        return -1;
    if (argc == 1)
        return print_usage(argv[0]), 0;
    if (argc < 2)
        return print_usage(argv[0]), -1;

    // Problem path
    bool rel_to_exe              = argv[1][0] == '^';
    std::string_view prob_path_s = argv[1] + static_cast<ptrdiff_t>(rel_to_exe);
    std::string_view prob_type;
    std::tie(prob_type, prob_path_s) = split_once(prob_path_s, ':');
    fs::path prob_path{prob_path_s};
    if (rel_to_exe)
        prob_path = fs::canonical(fs::path(argv[0])).parent_path() / prob_path;

    std::vector<std::string_view> extra_opts;
    std::copy(argv + 2, argv + argc, std::back_inserter(extra_opts));

    // Check which solver to use
    std::string_view solver_name = "panoc", direction_name;
    alpaqa::params::set_params(solver_name, "method", extra_opts);
    std::tie(solver_name, direction_name) =
        alpaqa::params::split_key(solver_name, '.');

    std::string out_path = "-";
    alpaqa::params::set_params(out_path, "out", extra_opts);
    std::ofstream out_fstream;
    if (out_path != "-")
        if (out_fstream.open(out_path); !out_fstream)
            throw std::runtime_error("Unable to open '" + out_path + "'");
    std::ostream &os = out_fstream.is_open() ? out_fstream : std::cout;

    // Load the problem
    auto problem = load_problem(prob_type, prob_path.parent_path(),
                                prob_path.filename(), extra_opts);
    os << "Loaded problem " << problem.path.stem().c_str() << " from "
       << problem.path << "\nProvided functions:\n";
    alpaqa::print_provided_functions(os, problem.problem);
    os << std::endl;

    // Available solvers
    std::map<std::string_view, std::function<void(std::string_view)>> solvers {
        {"pantr",
         [&](std::string_view direction) {
             do_experiment_with_direction<alpaqa::PANTRSolver>(
                 direction, problem, os, extra_opts);
         }},
            {"panoc",
             [&](std::string_view direction) {
                 do_experiment_with_direction<alpaqa::PANOCSolver>(
                     direction, problem, os, extra_opts);
             }},
            {"zerofpr",
             [&](std::string_view direction) {
                 do_experiment_with_direction<alpaqa::ZeroFPRSolver>(
                     direction, problem, os, extra_opts);
             }},
#if WITH_IPOPT
            {"ipopt",
             [&](std::string_view direction) {
                 if (!direction.empty())
                     throw std::runtime_error(
                         "Ipopt does not support directions.");
                 auto solver = make_ipopt_solver(extra_opts);
                 do_experiment(problem, solver, os, extra_opts);
             }},
#endif
#if WITH_LBFGSPP
            {"lbfgsbpp",
             [&](std::string_view direction) {
                 if (!direction.empty())
                     throw std::runtime_error(
                         "LBFGSB++ does not support directions.");
                 using InnerSolver = InnerLBFGSppSolver;
                 auto solver       = make_solver<InnerSolver>(extra_opts);
                 do_experiment(problem, solver, os, extra_opts);
             }},
#endif
#if WITH_LBFGSB
            {"lbfgsb", [&](std::string_view direction) {
                 if (!direction.empty())
                     throw std::runtime_error(
                         "Ipopt does not support directions.");
                 using InnerSolver = InnerLBFGSBSolver;
                 auto solver       = make_solver<InnerSolver>(extra_opts);
                 do_experiment(problem, solver, os, extra_opts);
             }},
#endif
    };

    // Run experiment
    auto solver_it = solvers.find(solver_name);
    if (solver_it != solvers.end())
        solver_it->second(direction_name);
    else
        throw std::invalid_argument(
            "Unknown solver '" + std::string(solver_name) + "'\n" +
            "  Available solvers: " +
            format_string_list(solvers, [](const auto &x) { return x.first; }));

} catch (std::exception &e) {
    std::cerr << "Error: " << demangled_typename(typeid(e)) << ":\n  "
              << e.what() << std::endl;
    return -1;
}
