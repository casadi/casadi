#pragma once

#include <alpaqa/config/config.hpp>

#include "options.hpp"
#include "problem.hpp"
#include "results.hpp"
#include "solver-driver.hpp"

solver_func_t make_panoc_driver(std::string_view direction, Options &opts);
solver_func_t make_zerofpr_driver(std::string_view direction, Options &opts);
