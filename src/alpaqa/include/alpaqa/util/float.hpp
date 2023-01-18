#pragma once

#include <alpaqa/util/quadmath/quadmath.hpp>

#include <cmath>
#include <complex>
#include <concepts>
#include <type_traits>

namespace alpaqa {

template <typename T>
struct is_complex_float : std::false_type {};

template <std::floating_point T>
struct is_complex_float<std::complex<T>> : std::true_type {};

template <class T>
inline constexpr bool is_complex_float_v = is_complex_float<T>::value;

template <typename T>
concept float_or_complex_float =
    std::floating_point<T> || is_complex_float_v<T>;

} // namespace alpaqa