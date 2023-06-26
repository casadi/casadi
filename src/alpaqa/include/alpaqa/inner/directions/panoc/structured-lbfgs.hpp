#pragma once

#include <alpaqa/accelerators/lbfgs.hpp>
#include <alpaqa/config/config.hpp>
#include <alpaqa/inner/directions/panoc-direction-update.hpp>
#include <alpaqa/inner/internal/panoc-helpers.hpp>
#include <alpaqa/problem/type-erased-problem.hpp>
#include <optional>

namespace alpaqa {

/// Parameters for the @ref StructuredLBFGSDirection class.
template <Config Conf>
struct StructuredLBFGSDirectionParams {
    USING_ALPAQA_CONFIG(Conf);
    /// Set this option to a nonzero value to include the Hessian-vector product
    /// @f$ \nabla^2_{x_\mathcal{J}x_\mathcal{K}}\psi(x) q_\mathcal{K} @f$ from
    /// equation 12b in @cite pas2022alpaqa, scaled by this parameter.
    /// Set it to zero to leave out that term.
    real_t hessian_vec_factor = 0;
    /// If @ref hessian_vec_factor is nonzero, set this option to true to
    /// approximate that term using finite differences instead of using AD.
    bool hessian_vec_finite_differences = true;
    /// If @ref hessian_vec_factor is nonzero and
    /// @ref hessian_vec_finite_differences is true, set this option to true to
    /// compute the exact Hessian of the augmented Lagrangian, false to
    /// approximate it using the Hessian of the Lagrangian.
    bool full_augmented_hessian = true;
    enum FailurePolicy {
        /// If L-BFGS fails, propagate the failure and tell PANOC that no
        /// accelerated step is available, causing it to accept the projected
        /// gradient step instead.
        FallbackToProjectedGradient,
        /// If L-BFGS fails, return @f$ q_\mathcal{J} =
        /// -\gamma\nabla_{x_\mathcal{J}}\psi(x^k)
        /// -\gamma\nabla^2_{x_\mathcal{J}x_\mathcal{K}}\psi(x)
        /// q_\mathcal{K} @f$ as the accelerated step (effectively approximating
        /// @f$ \nabla_{x_\mathcal{J}x_\mathcal{J}} \approx \gamma I @f$).
        UseScaledLBFGSInput,
    }
    /// What to do when L-BFGS failed (e.g. if there were no pairs (s, y) with
    /// positive curvature).
    failure_policy = FallbackToProjectedGradient;
};

/// @ingroup grp_DirectionProviders
template <Config Conf = DefaultConfig>
struct StructuredLBFGSDirection {
    USING_ALPAQA_CONFIG(Conf);
    using Problem           = TypeErasedProblem<config_t>;
    using LBFGS             = alpaqa::LBFGS<config_t>;
    using AcceleratorParams = typename LBFGS::Params;
    using DirectionParams   = StructuredLBFGSDirectionParams<config_t>;

    struct Params {
        AcceleratorParams accelerator = {};
        DirectionParams direction     = {};
    };

    StructuredLBFGSDirection() = default;
    StructuredLBFGSDirection(const Params &params)
        : lbfgs(params.accelerator), direction_params(params.direction) {}
    StructuredLBFGSDirection(const typename LBFGS::Params &params,
                             const DirectionParams &directionparams = {})
        : lbfgs(params), direction_params(directionparams) {}
    StructuredLBFGSDirection(const LBFGS &lbfgs,
                             const DirectionParams &directionparams = {})
        : lbfgs(lbfgs), direction_params(directionparams) {}
    StructuredLBFGSDirection(LBFGS &&lbfgs,
                             const DirectionParams &directionparams = {})
        : lbfgs(std::move(lbfgs)), direction_params(directionparams) {}

    /// @see @ref PANOCDirection::initialize
    void initialize(const Problem &problem, crvec y, crvec Σ, real_t γ_0,
                    crvec x_0, crvec x̂_0, crvec p_0, crvec grad_ψx_0);

    /// @see @ref PANOCDirection::has_initial_direction
    bool has_initial_direction() const { return false; }

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
    bool apply(real_t γₖ, crvec xₖ, crvec x̂ₖ, crvec pₖ, crvec grad_ψxₖ,
               rvec qₖ) const;

    /// @see @ref PANOCDirection::changed_γ
    void changed_γ([[maybe_unused]] real_t γₖ, [[maybe_unused]] real_t old_γₖ) {
        // Nothing, Hessian approximation is independent of step size
    }

    /// @see @ref PANOCDirection::reset
    void reset() { lbfgs.reset(); }

    /// @see @ref PANOCDirection::get_name
    std::string get_name() const {
        return "StructuredLBFGSDirection<" + std::string(config_t::get_name()) +
               '>';
    }

    auto get_params() const {
        return std::tie(lbfgs.get_params(), direction_params);
    }

  private:
    using Helpers = detail::PANOCHelpers<config_t>;

    const Problem *problem = nullptr;
#ifndef _WIN32
    std::optional<crvec> y = std::nullopt;
    std::optional<crvec> Σ = std::nullopt;
#else
    std::optional<vec> y = std::nullopt;
    std::optional<vec> Σ = std::nullopt;
#endif

    LBFGS lbfgs;
    mutable indexvec J_sto;
    mutable vec HqK;
    mutable vec work_n;
    mutable vec work_n2;
    mutable vec work_m;

    void approximate_hessian_vec_term(crvec xₖ, crvec grad_ψxₖ, rvec qₖ,
                                      crindexvec J) const;

  public:
    DirectionParams direction_params;
};

ALPAQA_EXPORT_EXTERN_TEMPLATE(struct, StructuredLBFGSDirection, DefaultConfig);
ALPAQA_EXPORT_EXTERN_TEMPLATE(struct, StructuredLBFGSDirection, EigenConfigf);
ALPAQA_EXPORT_EXTERN_TEMPLATE(struct, StructuredLBFGSDirection, EigenConfigd);
ALPAQA_EXPORT_EXTERN_TEMPLATE(struct, StructuredLBFGSDirection, EigenConfigl);
#ifdef ALPAQA_WITH_QUAD_PRECISION
ALPAQA_EXPORT_EXTERN_TEMPLATE(struct, StructuredLBFGSDirection, EigenConfigq);
#endif

} // namespace alpaqa
