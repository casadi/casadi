#include <alpaqa/dl/dl-problem.hpp>
#include <alpaqa/params/params.hpp>
#include <alpaqa/problem/problem-with-counters.hpp>
#include <alpaqa/problem/type-erased-problem.hpp>
#include <alpaqa/util/io/csv.hpp>
#if ALPAQA_HAVE_CASADI
#include <alpaqa/casadi/CasADiProblem.hpp>
#endif

#include <filesystem>
#include <fstream>
#include <mutex>
#include <optional>
#include <span>
#include <stdexcept>
#include <string>

#include "options.hpp"
#include "problem.hpp"

namespace fs = std::filesystem;

namespace {

USING_ALPAQA_CONFIG(alpaqa::DefaultConfig);

std::string get_prefix_option(std::span<const std::string_view> prob_opts) {
    std::string prefix          = "alpaqa_problem";
    std::string_view prefix_key = "prefix=";
    auto prefix_it              = std::find_if(
        prob_opts.rbegin(), prob_opts.rend(),
        [&](std::string_view opt) { return opt.starts_with(prefix_key); });
    if (prefix_it != prob_opts.rend())
        prefix = prefix_it->substr(prefix_key.size());
    return prefix;
}

} // namespace
namespace {
void load_initial_guess(Options &opts, LoadedProblem &problem) {
    const auto n = problem.problem.get_n(), m = problem.problem.get_m();
    alpaqa::params::vec_from_file<config_t> x0{n}, y0{m}, w0{2 * n};
    set_params(x0, "x0", opts);
    if (x0.value)
        problem.initial_guess_x = std::move(*x0.value);
    set_params(y0, "mul_g0", opts);
    if (y0.value)
        problem.initial_guess_y = std::move(*y0.value);
    set_params(w0, "mul_x0", opts);
    if (w0.value)
        problem.initial_guess_w = std::move(*w0.value);
}
} // namespace

LoadedProblem load_problem(std::string_view type, const fs::path &dir,
                           const fs::path &file, Options &opts) {
    USING_ALPAQA_CONFIG(alpaqa::DefaultConfig);
    // Isolate problem-specific options
    std::vector<std::string_view> prob_opts;
    std::string_view prob_prefix = "problem.";
    auto options                 = opts.options();
    auto used                    = opts.used();
    for (auto opt = options.begin(); opt != options.end(); ++opt) {
        if (opt->starts_with(prob_prefix)) {
            prob_opts.push_back(opt->substr(prob_prefix.size()));
            used.begin()[opt - options.begin()] = true;
        }
    }
    // Load problem
    auto full_path = dir / file;
    if (type == "dl" || type.empty()) {
        using TEProblem  = alpaqa::TypeErasedProblem<config_t>;
        using DLProblem  = alpaqa::dl::DLProblem;
        using CntProblem = alpaqa::ProblemWithCounters<DLProblem>;
        auto prefix      = get_prefix_option(prob_opts);
        std::any dl_opt  = std::span{prob_opts};
        LoadedProblem problem{
            .problem = TEProblem::make<CntProblem>(
                std::in_place, full_path.c_str(), prefix, &dl_opt),
            .abs_path = fs::absolute(full_path),
            .path     = full_path,
        };
        auto &cnt_problem   = problem.problem.as<CntProblem>();
        problem.evaluations = cnt_problem.evaluations;
        load_initial_guess(opts, problem);
        return problem;
    } else if (type == "cs") {
#if ALPAQA_HAVE_CASADI
        static std::mutex mtx;
        std::unique_lock lck{mtx};
        using TEProblem  = alpaqa::TypeErasedProblem<config_t>;
        using CsProblem  = alpaqa::CasADiProblem<config_t>;
        using CntProblem = alpaqa::ProblemWithCounters<CsProblem>;
        LoadedProblem problem{
            .problem =
                TEProblem::make<CntProblem>(std::in_place, full_path.c_str()),
            .abs_path = fs::absolute(full_path),
            .path     = full_path,
        };
        lck.unlock();
        auto &cnt_problem   = problem.problem.as<CntProblem>();
        auto &cs_problem    = cnt_problem.problem;
        problem.evaluations = cnt_problem.evaluations;
        auto param_size     = cs_problem.param.size();
        alpaqa::params::set_params(cs_problem.param, "param", prob_opts);
        if (cs_problem.param.size() != param_size)
            throw alpaqa::params::invalid_param(
                "Incorrect problem parameter size: got " +
                std::to_string(cs_problem.param.size()) + ", should be " +
                std::to_string(param_size));
        load_initial_guess(opts, problem);
        return problem;
#else
        throw std::logic_error(
            "This version of alpaqa was compiled without CasADi support");
#endif
    }
    throw std::invalid_argument("Unknown problem type '" + std::string(type) +
                                "'");
}
