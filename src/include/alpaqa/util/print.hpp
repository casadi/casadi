#pragma once

#include <alpaqa/util/float.hpp>

#include <iosfwd>
#include <string>

#include <Eigen/Core>

namespace alpaqa {

std::string float_to_str(std::floating_point auto value);

template <float_or_complex_float F>
std::ostream &print_matlab(std::ostream &os, const Eigen::MatrixX<F> &M);
template <float_or_complex_float F>
std::ostream &print_matlab(std::ostream &os, const Eigen::VectorX<F> &M);

} // namespace alpaqa