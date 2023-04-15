#include <alpaqa/config/config.hpp>
#include <alpaqa/problem/kkt-error.hpp>
#include <alpaqa/util/demangled-typename.hpp>
#include <alpaqa/util/print.hpp>

#include "ipopt-driver.hpp"
#include "lbfgsb-driver.hpp"
#include "options.hpp"
#include "panoc-driver.hpp"
#include "pantr-driver.hpp"
#include "results.hpp"
#include "solver-driver.hpp"
#include "util.hpp"

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <random>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <tuple>
#include <type_traits>
namespace fs = std::filesystem;

USING_ALPAQA_CONFIG(alpaqa::DefaultConfig);

void print_usage(const char *a0) {
    std::cout << "Usage: " << a0
              << " [<problem-type>:][<path>/]<name> [method=<solver>] "
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

std::ostream &get_output_stream(Options &opts, std::ofstream &out_fstream) {
    std::string out_path = "-";
    set_params(out_path, "out", opts);
    if (out_path != "-")
        if (out_fstream.open(out_path); !out_fstream)
            throw std::runtime_error("Unable to open '" + out_path + "'");
    return out_fstream.is_open() ? out_fstream : std::cout;
}

std::string_view get_output_paths(Options &opts) {
    std::string_view sol_path;
    set_params(sol_path, "sol", opts);
    return sol_path;
}

auto get_problem_path(const char *const *argv) {
    bool rel_to_exe              = argv[1][0] == '^';
    std::string_view prob_path_s = argv[1] + static_cast<ptrdiff_t>(rel_to_exe);
    std::string_view prob_type;
    std::tie(prob_type, prob_path_s) = split_once(prob_path_s, ':');
    fs::path prob_path{prob_path_s};
    if (rel_to_exe)
        prob_path = fs::canonical(fs::path(argv[0])).parent_path() / prob_path;
    return std::make_tuple(std::move(prob_path), prob_type);
}

auto get_solver_builder(Options &opts) {
    std::string_view method = "panoc", direction;
    set_params(method, "method", opts);
    std::tie(method, direction) = alpaqa::params::split_key(method, '.');
    // Dictionary of available solver builders
    std::map<std::string_view, solver_builder_func_t> solvers{
        {"panoc", make_panoc_driver},   {"zerofpr", make_zerofpr_driver},
        {"pantr", make_pantr_driver},
#ifdef WITH_LBFGSB
        {"lbfgsb", make_lbfgsb_driver},
#endif
#ifdef WITH_IPOPT
        {"ipopt", make_ipopt_driver},
#endif
    };
    // Find the selected solver builder
    auto solver_it = solvers.find(method);
    if (solver_it == solvers.end())
        throw std::invalid_argument(
            "Unknown solver '" + std::string(method) + "'\n" +
            "  Available solvers: " +
            format_string_list(solvers, [](const auto &x) { return x.first; }));
    return std::make_tuple(std::move(solver_it->second), direction);
}

void store_solution(const fs::path &sol_output_dir, std::ostream &os,
                    BenchmarkResults &results, std::span<const char *> argv) {
    const auto &sol_res = results.solver_results;
    auto timestamp_str  = std::to_string(results.timestamp);
    auto rnd_str        = random_hex_string(std::random_device());
    auto name           = results.problem.path.stem().string();
    auto suffix         = '_' + name + '_' + timestamp_str + '_' + rnd_str;
    fs::create_directories(sol_output_dir);
    std::array solutions{
        std::tuple{"solution", "sol_x", &sol_res.solution},
        std::tuple{"multipliers for g", "mul_g", &sol_res.multipliers},
        std::tuple{"multipliers for x", "mul_x", &sol_res.multipliers_bounds},
    };
    for (auto [name, fname, value] : solutions) {
        if (value->size() == 0)
            continue;
        auto pth = sol_output_dir / (std::string(fname) + suffix + ".csv");
        os << "Writing " << name << " to " << pth << std::endl;
        std::ofstream output_file(pth);
        alpaqa::print_csv(output_file, *value);
    }
    auto pth = sol_output_dir / ("cmdline" + suffix + ".txt");
    os << "Writing arguments to " << pth << std::endl;
    std::ofstream output_file(pth);
    for (const char *arg : argv)
        output_file << std::quoted(arg, '\'') << ' ';
    output_file << '\n';
}

int main(int argc, const char *argv[]) try {
    // Check command line options
    if (argc < 1)
        return -1;
    if (argc == 1)
        return print_usage(argv[0]), 0;
    if (argc < 2)
        return print_usage(argv[0]), -1;
    std::span args{argv, static_cast<size_t>(argc)};
    Options opts{argc - 2, argv + 2};

    // Check where to write the output to
    std::ofstream out_fstream;
    std::ostream &os = get_output_stream(opts, out_fstream);

    // Check which problem to load
    auto [prob_path, prob_type] = get_problem_path(argv);

    // Check which solver to use
    auto [solver_builder, direction] = get_solver_builder(opts);

    // Check output paths
    fs::path sol_output_dir = get_output_paths(opts);

    // Build solver
    auto solver = solver_builder(direction, opts);

    // Load problem
    os << "Loading problem " << prob_path << std::endl;
    auto problem = load_problem(prob_type, prob_path.parent_path(),
                                prob_path.filename(), opts);
    os << "Loaded problem " << problem.path.stem().c_str() << " from "
       << problem.path << "\nProvided functions:\n";
    alpaqa::print_provided_functions(os, problem.problem);
    os << std::endl;

    // Check options
    auto used       = opts.used();
    auto unused_opt = std::find(used.begin(), used.end(), false);
    auto unused_idx = static_cast<size_t>(unused_opt - used.begin());
    if (unused_opt != used.end())
        throw std::invalid_argument("Unused option: " +
                                    std::string(opts.options()[unused_idx]));

    // Solve
    auto solver_results = solver(problem, os);

    // Compute more statistics
    real_t f     = problem.problem.eval_f(solver_results.solution);
    auto kkt_err = alpaqa::compute_kkt_error(
        problem.problem, solver_results.solution, solver_results.multipliers);
    BenchmarkResults results{
        .problem          = problem,
        .solver_results   = solver_results,
        .objective        = f + solver_results.h,
        .smooth_objective = f,
        .error            = kkt_err,
        .options          = opts.options(),
        .timestamp        = timestamp_ms<std::chrono::system_clock>().count(),
    };

    // Print results
    print_results(os, results);

    // Store solution
    if (!sol_output_dir.empty())
        store_solution(sol_output_dir, os, results, args);

} catch (std::exception &e) {
    std::cerr << "Error: " << demangled_typename(typeid(e)) << ":\n  "
              << e.what() << std::endl;
    return -1;
}
