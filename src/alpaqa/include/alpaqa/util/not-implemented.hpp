#pragma once

#include <alpaqa/export.h>
#include <stdexcept>

namespace alpaqa {

struct ALPAQA_EXPORT not_implemented_error : std::logic_error {
    using std::logic_error::logic_error;
};

} // namespace alpaqa
