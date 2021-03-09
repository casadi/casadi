#pragma once

#include <panoc-alm/util/vec.hpp>

#include <functional>

namespace pa_ref {

using pa::real_t;
using pa::vec;

vec finite_diff(std::function<real_t(const vec &)> f, const vec &x);

} // namespace pa_ref