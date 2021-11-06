#pragma once

#include "vec.hpp"

namespace alpaqa {

struct Box {
    vec upperbound;
    vec lowerbound;
};

/// Project a vector onto a box.
/// @f[ \Pi_C(v) @f]
template <class Vec>
inline auto project(const Vec &v,  ///< [in] The vector to project
                    const Box &box ///< [in] The box to project onto
) {
    using binary_real_f = real_t (*)(real_t, real_t);
    return v.binaryExpr(box.lowerbound, binary_real_f(std::fmax))
        .binaryExpr(box.upperbound, binary_real_f(std::fmin));
}

/// Get the difference between the given vector and its projection.
/// @f[ v - \Pi_C(v) @f]
/// @warning    Beware catastrophic cancellation!
template <class Vec>
inline auto
projecting_difference(const Vec &v,  ///< [in] The vector to project
                      const Box &box ///< [in] The box to project onto
) {
    return v - project(v, box);
}

/// Get the distance squared between the given vector and its projection.
/// @f[ \left\| v - \Pi_C(v) \right\|_2^2 @f]
/// @warning    Beware catastrophic cancellation!
inline real_t dist_squared(crvec v,       ///< [in] The vector to project
                           const Box &box ///< [in] The box to project onto
) {
    return projecting_difference(v, box).squaredNorm();
}

/// Get the distance squared between the given vector and its projection in the
/// Σ norm.
/// @f[ \left\| v - \Pi_C(v) \right\|_\Sigma^2
/// = \left(v - \Pi_C(v)\right)^\top \Sigma \left(v - \Pi_C(v)\right) @f]
/// @warning    Beware catastrophic cancellation!
inline real_t dist_squared(crvec v,        ///< [in] The vector to project
                           const Box &box, ///< [in] The box to project onto
                           crvec Σ ///< [in] Diagonal matrix defining norm
) {
    // TODO: Does this allocate?
    //       Does it have dangling references to temporaries?
    auto d = v - project(v, box);
    return d.dot(Σ.asDiagonal() * d);
}

} // namespace alpaqa