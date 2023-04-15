#pragma once

#include "options.hpp"
#include "results.hpp"
#include <functional>

using solver_free_func_t = SolverResults(LoadedProblem &, std::ostream &);
using solver_func_t      = std::function<solver_free_func_t>;
using solver_builder_func_t =
    std::function<solver_func_t(std::string_view, Options &)>;
