#pragma once
#include <iosfwd>

namespace alpaqa {

/// Exit status of a numerical solver such as ALM or PANOC.
enum class SolverStatus {
    Unknown = 0, ///< Initial value.
    Converged,   ///< Converged and reached given tolerance.
    MaxTime,     ///< Maximum allowed execution time exceeded.
    MaxIter,     ///< Maximum number of iterations exceeded.
    NotFinite,   ///< Intermediate results were infinite or not-a-number.
    NoProgress,  ///< No progress was made in the last iteration.
    Interrupted, ///< Solver was interrupted by the user.
};

/// @related    SolverStatus
const char *enum_name(SolverStatus);
/// @related    SolverStatus
std::ostream &operator<<(std::ostream &, SolverStatus);

} // namespace alpaqa