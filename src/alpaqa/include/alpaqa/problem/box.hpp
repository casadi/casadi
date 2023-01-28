#pragma once

#include <alpaqa/config/config.hpp>
#include <alpaqa/util/float.hpp>

namespace alpaqa {

template <Config Conf = DefaultConfig>
struct Box {
    USING_ALPAQA_CONFIG(Conf);

    Box() : Box{0} {}
    Box(length_t n)
        : lowerbound{vec::Constant(n, -inf<config_t>)}, upperbound{
                                                            vec::Constant(n, +inf<config_t>)} {}

    static Box NaN(length_t n) {
        return Box{vec::Constant(n, alpaqa::NaN<config_t>),
                   vec::Constant(n, alpaqa::NaN<config_t>)};
    }
    static Box from_lower_upper(vec lower, vec upper) {
        return Box{std::move(lower), std::move(upper)};
    }

    vec lowerbound;
    vec upperbound;

  private:
    Box(vec lower, vec upper) : lowerbound{std::move(lower)}, upperbound{std::move(upper)} {}
};

/// Project a vector onto a box.
/// @f[ \Pi_C(v) @f]
template <Config Conf>
inline auto project(const auto &v,       ///< [in] The vector to project
                    const Box<Conf> &box ///< [in] The box to project onto
) {
    USING_ALPAQA_CONFIG(Conf);
    using binary_real_f = real_t (*)(real_t, real_t);
    return v.binaryExpr(box.lowerbound, binary_real_f(std::fmax))
        .binaryExpr(box.upperbound, binary_real_f(std::fmin));
}

/// Get the difference between the given vector and its projection.
/// @f[ v - \Pi_C(v) @f]
/// @warning    Beware catastrophic cancellation!
template <Config Conf>
inline auto                                //
projecting_difference(const auto &v,       ///< [in] The vector to project
                      const Box<Conf> &box ///< [in] The box to project onto
) {
    return v - project(v, box);
}

/// Get the distance squared between the given vector and its projection.
/// @f[ \left\| v - \Pi_C(v) \right\|_2^2 @f]
/// @warning    Beware catastrophic cancellation!

template <Config Conf>
inline auto dist_squared(const auto &v,       ///< [in] The vector to project
                         const Box<Conf> &box ///< [in] The box to project onto
) {
    return projecting_difference(v, box).squaredNorm();
}

/// Get the distance squared between the given vector and its projection in the
/// Σ norm.
/// @f[ \left\| v - \Pi_C(v) \right\|_\Sigma^2
/// = \left(v - \Pi_C(v)\right)^\top \Sigma \left(v - \Pi_C(v)\right) @f]
/// @warning    Beware catastrophic cancellation!
template <Config Conf>
inline auto dist_squared(const auto &v,        ///< [in] The vector to project
                         const Box<Conf> &box, ///< [in] The box to project onto
                         const auto &Σ         ///< [in] Diagonal matrix defining norm
                         ) -> real_t<Conf> {
    // TODO: Does this allocate?
    //       Does it have dangling references to temporaries?
    auto d = v - project(v, box);
    return d.dot(Σ.asDiagonal() * d);
}

} // namespace alpaqa