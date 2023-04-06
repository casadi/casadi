#pragma once

#include <alpaqa/config/config.hpp>

#include <chrono>
#include <optional>

namespace alpaqa {

template <Config Conf>
struct InnerSolveOptions {
    USING_ALPAQA_CONFIG(Conf);
    /// Return the final iterate and multipliers, even if the solver did not
    /// converge.
    bool always_overwrite_results = true;
    /// Maximum run time (in addition to the inner solver's own timeout).
    /// Zero means no timeout.
    std::optional<std::chrono::nanoseconds> max_time = std::nullopt;
    /// Desired tolerance (overrides the solver's own tolerance).
    /// Zero means no tolerance (use solver's own tolerance).
    real_t tolerance = 0;
    /// Output stream to print to.
    std::ostream *os = nullptr;
    /// The current iteration of the outer solver.
    unsigned outer_iter = 0;
    /// Call @ref TypeErasedProblem::check() before starting to solve.
    bool check = true;
};

} // namespace alpaqa
