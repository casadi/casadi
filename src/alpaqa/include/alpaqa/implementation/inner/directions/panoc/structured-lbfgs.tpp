#pragma once

#include <alpaqa/implementation/inner/panoc-helpers.tpp>
#include <alpaqa/inner/directions/panoc/structured-lbfgs.hpp>

namespace alpaqa {

template <Config Conf>
void StructuredLBFGSDirection<Conf>::initialize(
    const Problem &problem, crvec y, crvec Σ, [[maybe_unused]] real_t γ_0,
    [[maybe_unused]] crvec x_0, [[maybe_unused]] crvec x̂_0,
    [[maybe_unused]] crvec p_0, [[maybe_unused]] crvec grad_ψx_0) {
    if (!problem.provides_eval_inactive_indices_res_lna())
        throw std::invalid_argument(
            "Structured L-BFGS requires eval_inactive_indices_res_lna()");
    if (direction_params.hessian_vec_factor != 0 &&
        !direction_params.hessian_vec_finite_differences &&
        !direction_params.full_augmented_hessian &&
        !problem.provides_eval_hess_L_prod())
        throw std::invalid_argument(
            "Structured L-BFGS requires eval_hess_L_prod(). Alternatively, set "
            "hessian_vec_factor = 0 or hessian_vec_finite_differences = true.");
    if (direction_params.hessian_vec_factor != 0 &&
        !direction_params.hessian_vec_finite_differences &&
        direction_params.full_augmented_hessian &&
        !(problem.provides_eval_hess_L_prod() ||
          problem.provides_eval_hess_ψ_prod()))
        throw std::invalid_argument(
            "Structured L-BFGS requires _eval_hess_ψ_prod() or "
            "eval_hess_L_prod(). Alternatively, set "
            "hessian_vec_factor = 0 or hessian_vec_finite_differences = true.");
    if (direction_params.hessian_vec_factor != 0 &&
        !direction_params.hessian_vec_finite_differences &&
        direction_params.full_augmented_hessian &&
        !problem.provides_eval_hess_ψ_prod() &&
        !(problem.provides_get_box_D() && problem.provides_eval_grad_gi()))
        throw std::invalid_argument(
            "Structured L-BFGS requires either eval_hess_ψ_prod() or "
            "get_box_D() and eval_grad_gi(). Alternatively, set "
            "hessian_vec_factor = 0, hessian_vec_finite_differences = true, or "
            "full_augmented_hessian = false.");
    // Store references to problem and ALM variables
    this->problem = &problem;
    this->y.emplace(y);
    this->Σ.emplace(Σ);
    // Allocate workspaces
    const auto n = problem.get_n();
    const auto m = problem.get_m();
    lbfgs.resize(n);
    J_sto.resize(n);
    HqK.resize(n);
    if (direction_params.hessian_vec_finite_differences) {
        work_n.resize(n);
        work_n2.resize(n);
        work_m.resize(m);
    } else if (direction_params.full_augmented_hessian) {
        work_n.resize(n);
        work_m.resize(m);
    }
}

template <Config Conf>
bool StructuredLBFGSDirection<Conf>::apply(real_t γₖ, crvec xₖ,
                                           [[maybe_unused]] crvec x̂ₖ, crvec pₖ,
                                           crvec grad_ψxₖ, rvec qₖ) const {
    const auto n = problem->get_n();

    // Find inactive indices J
    auto nJ = problem->eval_inactive_indices_res_lna(γₖ, xₖ, grad_ψxₖ, J_sto);
    auto J = J_sto.topRows(nJ);

    // There are no inactive indices J
    if (nJ == 0) {
        // No free variables, no Newton step possible
        return false; // Simply use the projection step
    }
    // There are inactive indices J
    if (J.size() == n) { // There are no active indices K
        // If all indices are free, we can use standard L-BFGS,
        qₖ = (real_t(1) / γₖ) * pₖ;
        return lbfgs.apply(qₖ, γₖ);
    }
    // There are active indices K
    qₖ = pₖ;
    if (direction_params.hessian_vec_factor != 0) {
        qₖ(J).setZero();
        approximate_hessian_vec_term(xₖ, grad_ψxₖ, qₖ, J);
        // Compute right-hand side of 6.1c
        qₖ(J) = (real_t(1) / γₖ) * pₖ(J) -
                direction_params.hessian_vec_factor * HqK(J);
    } else {
        qₖ(J) = (real_t(1) / γₖ) * pₖ(J);
    }
    // If there are active indices, we need the specialized version
    // that only applies L-BFGS to the inactive indices
    bool success = lbfgs.apply_masked(qₖ, γₖ, J);
    if (success)
        return true;
    // If L-BFGS application failed, qₖ(J) still contains
    // -∇ψ(x)(J) - HqK(J) or -∇ψ(x)(J), which is not a valid step.
    // A good alternative is to use H₀ = γI as an L-BFGS estimate.
    // This seems to perform better in practice than just falling back to a
    // projected gradient step.
    switch (direction_params.failure_policy) {
        case DirectionParams::FallbackToProjectedGradient: return success;
        case DirectionParams::UseScaledLBFGSInput:
            if (nJ == n)
                qₖ *= γₖ;
            else
                qₖ(J) *= γₖ;
            return true;
        default: return false;
    }
}

template <Config Conf>
void StructuredLBFGSDirection<Conf>::approximate_hessian_vec_term(
    crvec xₖ, crvec grad_ψxₖ, rvec qₖ, crindexvec J) const {
    const auto m = problem->get_m();
    // Either compute the Hessian-vector product using finite differences
    if (direction_params.hessian_vec_finite_differences) {
        Helpers::calc_augmented_lagrangian_hessian_prod_fd(
            *problem, xₖ, *y, *Σ, grad_ψxₖ, qₖ, HqK, work_n, work_n2, work_m);
    }
    // Or using an exact AD
    else {
        if (!direction_params.full_augmented_hessian) {
            // Compute the product with the Hessian of the Lagrangian
            problem->eval_hess_L_prod(xₖ, *y, 1, qₖ, HqK);
        } else {
            if (problem->provides_eval_hess_ψ_prod()) {
                // Compute the product with the Hessian of the augmented
                // Lagrangian
                problem->eval_hess_ψ_prod(xₖ, *y, *Σ, 1, qₖ, HqK);
            } else {
                // Compute the product with the Hessian of the Lagrangian
                problem->eval_hess_L_prod(xₖ, *y, 1, qₖ, HqK);
                // And then add the Hessian of the penalty terms, to get the
                // Hessian of the full augmented Lagrangian (if required)
                if (direction_params.full_augmented_hessian) {
                    assert(m == 0 || problem->provides_eval_grad_gi());
                    const auto &D = problem->get_box_D();
                    auto &g       = work_m;
                    problem->eval_g(xₖ, g);
                    for (index_t i = 0; i < m; ++i) {
                        real_t ζ = g(i) + (*y)(i) / (*Σ)(i);
                        bool inactive =
                            D.lowerbound(i) < ζ && ζ < D.upperbound(i);
                        if (not inactive) {
                            problem->eval_grad_gi(xₖ, i, work_n);
                            auto t = (*Σ)(i)*work_n.dot(qₖ);
                            // TODO: the dot product is more work than
                            //       strictly necessary (only over K)
                            for (auto j : J)
                                HqK(j) += work_n(j) * t;
                        }
                    }
                }
            }
        }
    }
}

} // namespace alpaqa
