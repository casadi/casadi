#pragma once

#include <alpaqa/accelerators/lbfgs.hpp>
#include <alpaqa/config/config.hpp>
#include <alpaqa/inner/directions/panoc-direction-update.hpp>
#include <alpaqa/inner/internal/panoc-helpers.hpp>
#include <alpaqa/problem/type-erased-problem.hpp>
#include <optional>

namespace alpaqa {

/// Parameters for the @ref StructuredLBFGS class.
template <Config Conf>
struct StructuredLBFGSDirectionParams {
    bool hessian_vec                    = true;
    bool hessian_vec_finite_differences = true;
    bool full_augmented_hessian         = true;
};

/// @ingroup grp_DirectionProviders
template <Config Conf = DefaultConfig>
struct StructuredLBFGS {
    USING_ALPAQA_CONFIG(Conf);
    using Problem = TypeErasedProblem<config_t>;
    using LBFGS   = alpaqa::LBFGS<config_t>;

    using DirectionParams = StructuredLBFGSDirectionParams<config_t>;

    StructuredLBFGS(const typename LBFGS::Params &params    = {},
                    const DirectionParams &direction_params = {})
        : lbfgs(params), direction_params(direction_params) {}

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
        return "StructuredLBFGS<" + std::string(config_t::get_name()) + '>';
    }

    auto get_params() const {
        return std::tie(lbfgs.get_params(), direction_params);
    }

  private:
    using indexstdvec = std::vector<index_t>;
    using Helpers     = detail::PANOCHelpers<config_t>;

    const Problem *problem = nullptr;
#ifndef _WIN32
    std::optional<crvec> y = std::nullopt;
    std::optional<crvec> Σ = std::nullopt;
#else
    std::optional<vec> y = std::nullopt;
    std::optional<vec> Σ = std::nullopt;
#endif

    LBFGS lbfgs;
    mutable indexstdvec J;
    mutable vec HqK;
    mutable vec work_n;
    mutable vec work_n2;
    mutable vec work_m;

    DirectionParams direction_params;
};

template <Config Conf>
struct PANOCDirection<StructuredLBFGS<Conf>> : StructuredLBFGS<Conf> {
    using StructuredLBFGS<Conf>::StructuredLBFGS;
    PANOCDirection(const StructuredLBFGS<Conf> &s) : StructuredLBFGS<Conf>{s} {}
    PANOCDirection(StructuredLBFGS<Conf> &&s)
        : StructuredLBFGS<Conf>{std::move(s)} {}
};

ALPAQA_EXPORT_EXTERN_TEMPLATE(struct, StructuredLBFGS, DefaultConfig);
ALPAQA_EXPORT_EXTERN_TEMPLATE(struct, StructuredLBFGS, EigenConfigf);
ALPAQA_EXPORT_EXTERN_TEMPLATE(struct, StructuredLBFGS, EigenConfigd);
ALPAQA_EXPORT_EXTERN_TEMPLATE(struct, StructuredLBFGS, EigenConfigl);
#ifdef ALPAQA_WITH_QUAD_PRECISION
ALPAQA_EXPORT_EXTERN_TEMPLATE(struct, StructuredLBFGS, EigenConfigq);
#endif

} // namespace alpaqa
