#include <alpaqa/inner/internal/panoc-stop-crit.hpp>

#include <ostream>

namespace alpaqa {

std::ostream &operator<<(std::ostream &os, PANOCStopCrit s) {
    return os << enum_name(s);
}

} // namespace alpaqa