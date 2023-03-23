#pragma once

#include <alpaqa/casadi/CasADiProblem.hpp>
#include <alpaqa/config/config.hpp>
#include <alpaqa/dl/dl-problem.hpp>
#include <alpaqa/problem/type-erased-problem.hpp>
#include "alpaqa/params/params.hpp"

#include <filesystem>
#include <iostream>
#include <iterator>
#include <mutex>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>
namespace fs = std::filesystem;

struct LoadedProblem {
    alpaqa::TypeErasedProblem<alpaqa::DefaultConfig> problem;
    fs::path abs_path;
    fs::path path;
    std::shared_ptr<alpaqa::EvalCounter> evaluations = nullptr;
};

inline LoadedProblem
load_problem(std::string_view type, const fs::path &dir, const fs::path &file,
             std::span<const std::string_view> extra_opts) {
    USING_ALPAQA_CONFIG(alpaqa::DefaultConfig);
    // Isolate problem-specific options
    std::vector<std::string_view> prob_opts;
    std::string_view prob_prefix = "problem.";
    std::copy_if(
        extra_opts.begin(), extra_opts.end(), std::back_inserter(prob_opts),
        [=](std::string_view opt) { return opt.starts_with(prob_prefix); });
    std::for_each(
        prob_opts.begin(), prob_opts.end(),
        [=](std::string_view &opt) { opt = opt.substr(prob_prefix.size()); });
    // Load problem
    auto full_path = dir / file;
    if (type == "cs") {
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
        return problem;
    } else if (type == "dl") {
        using TEProblem  = alpaqa::TypeErasedProblem<config_t>;
        using DLProblem  = alpaqa::dl::DLProblem;
        using CntProblem = alpaqa::ProblemWithCounters<DLProblem>;
        std::any dl_opt  = std::span{prob_opts};
        LoadedProblem problem{
            .problem = TEProblem::make<CntProblem>(
                std::in_place, full_path.c_str(), "benchmark_problem", &dl_opt),
            .abs_path = fs::absolute(full_path),
            .path     = full_path,
        };
        auto &cnt_problem   = problem.problem.as<CntProblem>();
        problem.evaluations = cnt_problem.evaluations;
        return problem;
    }
    throw std::invalid_argument("Unknown problem type '" + std::string(type) +
                                "'");
}
