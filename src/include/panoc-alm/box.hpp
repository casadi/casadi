#pragma once

#include "vec.hpp"

namespace pa {

struct Box {
    vec upperbound;
    vec lowerbound;
};

template <class Vec>
inline auto project(const Vec &v,  // in
                    const Box &box // in
) {
    using binary_real_f = real_t (*)(real_t, real_t);
    return v.binaryExpr(box.lowerbound, binary_real_f(std::fmax))
        .binaryExpr(box.upperbound, binary_real_f(std::fmin));
}

inline real_t dist_squared(const vec &v, const Box &box) {
    // TODO: does this allocate?
    auto d = v - project(v, box);
    return d.squaredNorm();
}

inline real_t dist_squared(const vec &v, const Box &box, const vec &Σ) {
    // TODO: does this allocate?
    auto d  = v - project(v, box);
    auto Σd = (Σ.array() * d.array()).matrix();
    return Σd.dot(d);
}

} // namespace pa