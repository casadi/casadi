#pragma once

#include <alpaqa/config/config.hpp>
#include <alpaqa/problem/problem-counters.hpp>
#include <alpaqa/problem/type-erased-problem.hpp>

#include "options.hpp"

#include <filesystem>
#include <memory>
namespace fs = std::filesystem;

struct LoadedProblem {
    USING_ALPAQA_CONFIG(alpaqa::DefaultConfig);
    alpaqa::TypeErasedProblem<config_t> problem;
    fs::path abs_path;
    fs::path path;
    std::shared_ptr<alpaqa::EvalCounter> evaluations = nullptr;
    vec initial_guess_x = vec::Zero(problem.get_n()); /// Unknowns
    vec initial_guess_y = vec::Zero(problem.get_m()); /// Multipliers g
    vec initial_guess_w = alpaqa::null_vec<config_t>; /// Multipliers bounds
};

LoadedProblem load_problem(std::string_view type, const fs::path &dir,
                           const fs::path &file, Options &opts);
