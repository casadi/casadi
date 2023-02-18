#pragma once

#include <alpaqa/export.hpp>
#include <alpaqa/util/float.hpp>

#include <iosfwd>
#include <limits>
#include <string>

#include <Eigen/Core>

namespace alpaqa {

template <std::floating_point F>
ALPAQA_EXPORT std::string
float_to_str(F value, int precision = std::numeric_limits<F>::max_digits10);

template <float_or_complex_float F>
ALPAQA_EXPORT std::ostream &print_matlab(std::ostream &os,
                                         const Eigen::MatrixX<F> &M);
template <float_or_complex_float F>
ALPAQA_EXPORT std::ostream &print_matlab(std::ostream &os,
                                         const Eigen::VectorX<F> &M);

template <float_or_complex_float F>
ALPAQA_EXPORT std::ostream &print_python(std::ostream &os,
                                         const Eigen::MatrixX<F> &M,
                                         std::string_view end = "\n");
template <float_or_complex_float F>
ALPAQA_EXPORT std::ostream &print_python(std::ostream &os,
                                         const Eigen::VectorX<F> &M,
                                         std::string_view end = "\n");

} // namespace alpaqa