#pragma once

#include <Eigen/Core>

namespace pa {

using real_t         = double; // TODO: make template?
using realvec        = Eigen::Matrix<real_t, Eigen::Dynamic, 1>;
using realmat        = Eigen::Matrix<real_t, Eigen::Dynamic, Eigen::Dynamic>;
using vec            = realvec;
using rvec           = Eigen::Ref<vec>;
using crvec          = Eigen::Ref<const vec>;
using mat            = realmat;
using rmat           = Eigen::Ref<mat>;
using crmat          = Eigen::Ref<const mat>;
constexpr real_t inf = std::numeric_limits<real_t>::infinity();
constexpr real_t NaN = std::numeric_limits<real_t>::quiet_NaN();

namespace vec_util {

/// Get the Σ norm squared of a given vector, with Σ a diagonal matrix.
/// @returns @f$ \langle v, \Sigma v \rangle @f$
template <class V, class M>
auto norm_squared_weighted(V &&v, M &&Σ) {
    return v.dot(Σ.asDiagonal() * v);
}

/// Get the maximum or infinity-norm of the given vector.
/// @returns @f$ \left\|v\right\|_\infty @f$
template <class Vec>
real_t norm_inf(const Vec &v) {
    return v.template lpNorm<Eigen::Infinity>();
}

} // namespace vec_util

} // namespace pa