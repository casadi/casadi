#pragma once

#include "box.hpp"

#include <functional>

namespace pa {

struct Problem {
    unsigned int n; ///< Number of decision variables, dimension of x
    unsigned int m; ///< Number of constraints, dimension of g(x) and z
    Box C;
    Box D;

    using f_sig      = real_t(const vec &); // x      [in]
    using grad_f_sig = void(const vec &,    // x      [in]
                            vec &);         // ∇f(x)  [out]
    using g_sig      = void(const vec &,    // x      [in]
                            vec &);         // g(x)   [out]
    using grad_g_sig = void(const vec &,    // x      [in]
                            const vec &,    // y      [in]
                            vec &);         // ∇g(x)y [out]

    std::function<f_sig> f;
    std::function<grad_f_sig> grad_f;
    std::function<g_sig> g;
    std::function<grad_g_sig> grad_g;

    void validate();
};

} // namespace pa
