#pragma once

#ifdef ALPAQA_WITH_QUAD_PRECISION

#include <alpaqa/util/quadmath/quadmath.hpp>

#include <cassert>
#include <ostream>

namespace std {
inline ostream &operator<<(ostream &os, __float128 f) {
    char buf[128];
    auto precision         = static_cast<int>(os.precision());
    [[maybe_unused]] int n = quadmath_snprintf(buf, sizeof(buf), "%#.*Qg", precision, f);
    assert((size_t)n < sizeof buf);
    return os << buf;
}
} // namespace std

#endif