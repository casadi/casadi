#include <alpaqa/util/solverstatus.hpp>

#include <ostream>
#include <stdexcept>

namespace alpaqa {

const char *enum_name(SolverStatus s) {
    using Status = SolverStatus;
    switch (s) {
        case Status::Unknown: return "Unknown";
        case Status::Converged: return "Converged";
        case Status::MaxTime: return "MaxTime";
        case Status::MaxIter: return "MaxIter";
        case Status::NotFinite: return "NotFinite";
        case Status::NoProgress: return "NoProgress";
        case Status::Interrupted: return "Interrupted";
    }
    throw std::out_of_range("invalid value for alpaqa::SolverStatus");
}

std::ostream &operator<<(std::ostream &os, SolverStatus s) {
    return os << enum_name(s);
}

} // namespace alpaqa