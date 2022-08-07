#pragma once

#include <stdexcept>

namespace alpaqa {

struct not_implemented_error : std::logic_error {
    using std::logic_error::logic_error;
};

} // namespace alpaqa
