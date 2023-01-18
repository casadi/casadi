#include <alpaqa/inner/internal/solverstatus.hpp>

#include <ostream>

namespace alpaqa {

std::ostream &operator<<(std::ostream &os, SolverStatus s) {
    return os << enum_name(s);
}

} // namespace alpaqa