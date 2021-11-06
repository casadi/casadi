#pragma once

#include <Eigen/Core>

namespace alpaqa {

/// Default floating point type
using real_t         = double; // TODO: make template?
/// Default type for floating point vectors.
using realvec        = Eigen::Matrix<real_t, Eigen::Dynamic, 1>;
/// Default type for floating point matrices.
using realmat        = Eigen::Matrix<real_t, Eigen::Dynamic, Eigen::Dynamic>;
/// Default type for vectors.
using vec            = realvec;
/// Default type for mutable references to vectors.
using rvec           = Eigen::Ref<vec>;
/// Default type for immutable references to vectors.
using crvec          = Eigen::Ref<const vec>;
/// Default type for matrices.
using mat            = realmat;
/// Default type for mutable references to matrices.
using rmat           = Eigen::Ref<mat>;
/// Default type for immutable references to matrices.
using crmat          = Eigen::Ref<const mat>;
/// @f$ \infty @f$
constexpr real_t inf = std::numeric_limits<real_t>::infinity();
/// Not a number.
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

/// Get the 1-norm of the given vector.
/// @returns @f$ \left\|v\right\|_1 @f$
template <class Vec>
real_t norm_1(const Vec &v) {
    return v.template lpNorm<1>();
}

} // namespace vec_util

} // namespace alpaqa