#pragma once

#include <panoc-alm/util/vec.hpp>

#include <functional>

namespace pa_ref {

using pa::crvec;
using pa::real_t;
using pa::vec;

vec finite_diff(std::function<real_t(crvec)> f, crvec x);

} // namespace pa_ref