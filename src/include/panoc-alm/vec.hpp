#pragma once

#include <Eigen/Core>

namespace pa {

using real_t = double; // TODO: make template?
using realvec = Eigen::Matrix<real_t, Eigen::Dynamic, 1>;
using vec = realvec;

} // namespace pa