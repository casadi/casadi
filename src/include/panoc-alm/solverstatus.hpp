#pragma once
#include <iosfwd>

namespace pa {

enum class SolverStatus {
    Unknown = 0,
    Converged,
    MaxTime,
    MaxIter,
    NotFinite,
    Interrupted,
};

/// @related    SolverStatus
const char *enum_name(SolverStatus);
/// @related    SolverStatus
std::ostream &operator<<(std::ostream &, SolverStatus);

} // namespace pa