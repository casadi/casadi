#pragma once

#include <alpaqa/util/vec.hpp>

#include <functional>

namespace pa_ref {

using alpaqa::crvec;
using alpaqa::real_t;
using alpaqa::vec;

vec finite_diff(std::function<real_t(crvec)> f, crvec x);

} // namespace pa_ref