#pragma once

#include "problem.hpp"
#include <algorithm>
#include <limits>
#include <type_traits>
#include <vector>

namespace pa {

struct PANOCParams {
    struct {
        /// Relative step size for finite difference Lipschitz estimate
        real_t ε;
        /// Minimum step size for finite difference Lipschitz estimate
        real_t δ;
    } Lipschitz;
    unsigned lbfgs_mem;
    unsigned max_iter;

    void verify();
};

class PANOCSolver {
  public:
    using Params = PANOCParams;

    PANOCSolver(Params params) : params(params) {}

    void operator()(const Problem &problem, // in
                    vec &x,                 // inout
                    vec &z,                 // out
                    vec &y,                 // inout
                    vec &err_z,             // out
                    const vec &Σ,           // in
                    real_t ε);

    void proximal_gradient_step(const vec &x_in, const vec &z_in, vec &x_out);

  private:
    Params params;
};

namespace detail {

void project_y(vec &y,          // inout
               const vec &z_lb, // in
               const vec &z_ub, // in
               real_t M         // in
);

template <class Vec>
real_t norm_inf(const Vec &v) {
    return v.template lpNorm<Eigen::Infinity>();
}

} // namespace detail

} // namespace pa
