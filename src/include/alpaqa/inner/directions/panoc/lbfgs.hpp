#pragma once

#include <alpaqa/accelerators/lbfgs.hpp>
#include <alpaqa/inner/directions/panoc-direction-update.hpp>
#include <alpaqa/problem/type-erased-problem.hpp>

namespace alpaqa {

/// Parameters for the @ref PANOCDirection<LBFGS<Conf>> class.
template <Config Conf>
struct LBFGSDirectionParams {
    bool rescale_when_γ_changes = false;
};

/// @ingroup grp_DirectionProviders
template <Config Conf>
struct PANOCDirection<LBFGS<Conf>> {
    USING_ALPAQA_CONFIG(Conf);

    using Problem = TypeErasedProblem<config_t>;
    using LBFGS   = alpaqa::LBFGS<config_t>;

    LBFGS lbfgs;

    using DirectionParams = LBFGSDirectionParams<config_t>;
    struct Params : LBFGS::Params, DirectionParams {};

    PANOCDirection(const Params &params)
        : lbfgs(params), direction_params(params) {}
    PANOCDirection(const typename LBFGS::Params &params,
                   const DirectionParams &directionparams = {})
        : lbfgs(params), direction_params(directionparams) {}
    PANOCDirection(const LBFGS &lbfgs,
                   const DirectionParams &directionparams = {})
        : lbfgs(lbfgs), direction_params(directionparams) {}
    PANOCDirection(LBFGS &&lbfgs, const DirectionParams &directionparams = {})
        : lbfgs(std::move(lbfgs)), direction_params(directionparams) {}

    /// @see @ref PANOCDirection::initialize
    void initialize(const Problem &problem, [[maybe_unused]] crvec y,
                    [[maybe_unused]] crvec Σ, [[maybe_unused]] real_t γ_0,
                    [[maybe_unused]] crvec x_0, [[maybe_unused]] crvec x̂_0,
                    [[maybe_unused]] crvec p_0,
                    [[maybe_unused]] crvec grad_ψx_0) {
        lbfgs.resize(problem.get_n());
    }

    /// @see @ref PANOCDirection::has_initial_direction
    bool has_initial_direction() const { return false; }

    /// @see @ref PANOCDirection::update
    bool update([[maybe_unused]] real_t γₖ, [[maybe_unused]] real_t γₙₑₓₜ,
                crvec xₖ, crvec xₙₑₓₜ, crvec pₖ, crvec pₙₑₓₜ,
                [[maybe_unused]] crvec grad_ψxₖ,
                [[maybe_unused]] crvec grad_ψxₙₑₓₜ) {
        return lbfgs.update(xₖ, xₙₑₓₜ, pₖ, pₙₑₓₜ, LBFGS::Sign::Negative);
    }

    /// @see @ref PANOCDirection::apply
    bool apply([[maybe_unused]] real_t γₖ, [[maybe_unused]] crvec xₖ,
               [[maybe_unused]] crvec x̂ₖ, crvec pₖ,
               [[maybe_unused]] crvec grad_ψxₖ, rvec qₖ) const {
        qₖ = pₖ;
        return lbfgs.apply(qₖ, γₖ);
    }

    /// @see @ref PANOCDirection::changed_γ
    void changed_γ(real_t γₖ, real_t old_γₖ) {
        if (direction_params.rescale_when_γ_changes)
            lbfgs.scale_y(γₖ / old_γₖ);
        else
            lbfgs.reset();
    }

    /// @see @ref PANOCDirection::reset
    void reset() { lbfgs.reset(); }

    /// @see @ref PANOCDirection::get_name
    std::string get_name() const { return lbfgs.get_name(); }

    auto get_params() const {
        return std::tie(lbfgs.get_params(), direction_params);
    }

    DirectionParams direction_params;
};

} // namespace alpaqa