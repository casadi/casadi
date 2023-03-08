#pragma once

#include <alpaqa/export.h>

#include <iosfwd>
#include <stdexcept>

namespace alpaqa {

/// Exit status of a numerical solver such as ALM or PANOC.
enum class SolverStatus {
    Busy = 0,    ///< In progress.
    Converged,   ///< Converged and reached given tolerance.
    MaxTime,     ///< Maximum allowed execution time exceeded.
    MaxIter,     ///< Maximum number of iterations exceeded.
    NotFinite,   ///< Intermediate results were infinite or not-a-number.
    NoProgress,  ///< No progress was made in the last iteration.
    Interrupted, ///< Solver was interrupted by the user.
    Exception,   ///< An unexpected exception was thrown.
};

/// @related    SolverStatus
inline constexpr const char *enum_name(SolverStatus s) {
    using Status = SolverStatus;
    switch (s) {
        case Status::Busy: return "Busy";
        case Status::Converged: return "Converged";
        case Status::MaxTime: return "MaxTime";
        case Status::MaxIter: return "MaxIter";
        case Status::NotFinite: return "NotFinite";
        case Status::NoProgress: return "NoProgress";
        case Status::Interrupted: return "Interrupted";
        case Status::Exception: return "Exception";
        default:;
    }
    throw std::out_of_range("invalid value for alpaqa::SolverStatus");
}

/// @related    SolverStatus
ALPAQA_EXPORT std::ostream &operator<<(std::ostream &, SolverStatus);

} // namespace alpaqa