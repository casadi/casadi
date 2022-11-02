#pragma once

#include <alpaqa/accelerators/lbfgs.hpp>
#include <alpaqa/config/config.hpp>
#include <alpaqa/inner/directions/panoc-direction-update.hpp>
#include <alpaqa/inner/src/panoc-helpers.tpp>
#include <alpaqa/problem/type-erased-problem.hpp>
#include <optional>

namespace alpaqa {

/// @ingroup grp_DirectionProviders
template <Config Conf>
struct StructuredLBFGS {
    USING_ALPAQA_CONFIG(Conf);
    using Problem = TypeErasedProblem<config_t>;
    using LBFGS   = alpaqa::LBFGS<config_t>;

    struct ExtraParams {
        bool hessian_vec                    = true;
        bool hessian_vec_finite_differences = true;
        bool full_augmented_hessian         = true;
    };
    struct Params : LBFGS::Params, ExtraParams {};

    StructuredLBFGS(const Params &params)
        : lbfgs(params), extraparams(params) {}
    StructuredLBFGS(const typename LBFGS::Params &params,
                    const ExtraParams &extraparams = {})
        : lbfgs(params), extraparams(extraparams) {}

    /// @see @ref PANOCDirection::initializer
    void initialize(const Problem &problem, crvec y, crvec Σ,
                    [[maybe_unused]] real_t γ_0, [[maybe_unused]] crvec x_0,
                    [[maybe_unused]] crvec x̂_0, [[maybe_unused]] crvec p_0,
                    [[maybe_unused]] crvec grad_ψx_0) {
        if (!(problem.provides_get_box_C() && problem.provides_get_box_D()))
            throw std::invalid_argument(
                "Structured PANOC only supports box-constrained problems");
        // Store references to problem and ALM variables
        this->problem = &problem;
        this->y.emplace(y);
        this->Σ.emplace(Σ);
        // Allocate workspaces
        const auto n = problem.get_n();
        const auto m = problem.get_m();
        lbfgs.resize(n);
        J.reserve(static_cast<size_t>(n));
        HqK.resize(n);
        if (extraparams.hessian_vec_finite_differences) {
            work_n.resize(n);
            work_n2.resize(n);
            work_m.resize(m);
        } else if (extraparams.full_augmented_hessian) {
            work_n.resize(n);
            work_m.resize(m);
        }
    }

    /// @see @ref PANOCDirection::update
    bool update([[maybe_unused]] real_t γₖ, [[maybe_unused]] real_t γₙₑₓₜ,
                crvec xₖ, crvec xₙₑₓₜ, [[maybe_unused]] crvec pₖ,
                [[maybe_unused]] crvec pₙₑₓₜ, crvec grad_ψxₖ,
                crvec grad_ψxₙₑₓₜ) {
        const bool force = true;
        return lbfgs.update(xₖ, xₙₑₓₜ, grad_ψxₖ, grad_ψxₙₑₓₜ,
                            LBFGS::Sign::Positive, force);
    }

    /// @see @ref PANOCDirection::apply
    bool apply(real_t γₖ, crvec xₖ, [[maybe_unused]] crvec x̂ₖ, crvec pₖ,
               crvec grad_ψxₖ, rvec qₖ) const {
        const auto n  = problem->get_n();
        const auto m  = problem->get_m();
        const auto un = static_cast<std::make_unsigned_t<decltype(n)>>(n);
        const auto &C = problem->get_box_C();
        const auto &D = problem->get_box_D();

        // Find inactive indices J
        J.clear();
        for (index_t i = 0; i < n; ++i) {
            real_t gd = xₖ(i) - γₖ * grad_ψxₖ(i);
            if (gd <= C.lowerbound(i)) {        // i ∊ J̲ ⊆ K
                qₖ(i) = pₖ(i);                  //
            } else if (C.upperbound(i) <= gd) { // i ∊ J̅ ⊆ K
                qₖ(i) = pₖ(i);                  //
            } else {                            // i ∊ J
                J.push_back(i);
                qₖ(i) = extraparams.hessian_vec ? 0 : -grad_ψxₖ(i);
            }
        }

        if (not J.empty()) {      // There are inactive indices J
            if (J.size() == un) { // There are no active indices K
                qₖ = -grad_ψxₖ;
            } else if (extraparams.hessian_vec) { // There are active indices K
                if (extraparams.hessian_vec_finite_differences) {
                    Helpers::calc_augmented_lagrangian_hessian_prod_fd(
                        *problem, xₖ, *y, *Σ, grad_ψxₖ, qₖ, HqK, work_n,
                        work_n2, work_m);
                } else {
                    problem->eval_hess_L_prod(xₖ, *y, qₖ, HqK);
                    if (extraparams.full_augmented_hessian) {
                        auto &g = work_m;
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

                for (auto j : J) // Compute right-hand side of 6.1c
                    qₖ(j) = -grad_ψxₖ(j) - HqK(j);
            }

            // If all indices are inactive, we can use standard L-BFGS,
            // if there are active indices, we need the specialized version
            // that only applies L-BFGS to the inactive indices
            bool success = lbfgs.apply_masked(qₖ, γₖ, J);
            // If L-BFGS application failed, qₖ(J) still contains
            // -∇ψ(x)(J) - HqK(J) or -∇ψ(x)(J), which is not a valid step.
            // A good alternative is to use H₀ = γI as an L-BFGS estimate.
            // This seems to be better than just falling back to a projected
            // gradient step.
            if (not success) {
                if (J.size() == un)
                    qₖ *= γₖ;
                else
                    for (auto j : J)
                        qₖ(j) *= γₖ;
            }
        }
        return true;
    }

    /// @see @ref PANOCDirection::changed_γ
    void changed_γ([[maybe_unused]] real_t γₖ, [[maybe_unused]] real_t old_γₖ) {
        // Nothing, Hessian approximation is independent of step size
    }

    /// @see @ref PANOCDirection::reset
    void reset() { lbfgs.reset(); }

    /// @see @ref PANOCDirection::get_name
    std::string get_name() const {
        return "StructuredLBFGS<" + std::string(config_t::get_name()) + '>';
    }

    ExtraParams get_params() const { return lbfgs.get_params(); }

  private:
    using indexstdvec = std::vector<index_t>;
    using Helpers     = detail::PANOCHelpers<config_t>;

    const Problem *problem;
    std::optional<crvec> y;
    std::optional<crvec> Σ;

    LBFGS lbfgs;
    mutable indexstdvec J;
    mutable vec HqK;
    mutable vec work_n;
    mutable vec work_n2;
    mutable vec work_m;

    ExtraParams extraparams;
};

template <Config Conf>
struct PANOCDirection<StructuredLBFGS<Conf>> : StructuredLBFGS<Conf> {
    using StructuredLBFGS<Conf>::StructuredLBFGS;
};

} // namespace alpaqa
