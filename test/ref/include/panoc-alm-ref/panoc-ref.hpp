#pragma once

#include <panoc-alm/panoc.hpp>

/// Reference implementations that are more readable than the optimized 
/// implementations, used for tests as well
namespace pa_ref {

using pa::PANOCParams;
using pa::Problem;
using pa::real_t;
using pa::vec;

class PANOCSolver {
  public:
    using Params = PANOCParams;

    using Stats = pa::PANOCSolver::Stats;

    PANOCSolver(Params params) : params(params) {}

    Stats operator()(const Problem &problem, // in
                     vec &x,                 // inout
                     vec &z,                 // out
                     vec &y,                 // inout
                     vec &err_z,             // out
                     const vec &Σ,           // in
                     real_t ε);

  private:
    Params params;
};

namespace detail {

vec eval_g(const Problem &p, const vec &x);
vec eval_ẑ(const Problem &p, const vec &x, const vec &y, const vec &Σ);
vec eval_ŷ(const Problem &p, const vec &x, const vec &y, const vec &Σ);

real_t eval_ψ(const Problem &p, const vec &x, const vec &y, const vec &Σ);
vec eval_grad_ψ(const Problem &p, const vec &x, const vec &y, const vec &Σ);

vec gradient_step(const Problem &p, const vec &x, const vec &y, const vec &Σ,
                  real_t γ);
vec T_γ(const Problem &p, const vec &x, const vec &y, const vec &Σ, real_t γ);

real_t eval_φ(const Problem &p, const vec &x, const vec &y, const vec &Σ,
              real_t γ);

real_t estimate_lipschitz(const Problem &p, const vec &x, const vec &y,
                          const vec &Σ, const PANOCParams &params);

real_t calc_error_stop_crit(const Problem &p, const vec &xₖ, const vec &x̂ₖ,
                            const vec &y, const vec &Σ, real_t γ);

bool lipschitz_check(const Problem &p, const vec &xₖ, const vec &x̂ₖ,
                     const vec &y, const vec &Σ, real_t γ, real_t L);

bool linesearch_condition(const Problem &p, const vec &xₖ, const vec &xₖ₊₁,
                          const vec &rₖ, const vec &y, const vec &Σ, real_t γ,
                          real_t σ);

} // namespace detail

} // namespace pa_ref
