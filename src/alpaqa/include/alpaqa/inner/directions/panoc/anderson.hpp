#pragma once

#include <alpaqa/accelerators/anderson.hpp>
#include <alpaqa/inner/directions/panoc-direction-update.hpp>
#include <alpaqa/problem/type-erased-problem.hpp>

namespace alpaqa {

/// Parameters for the @ref AndersonDirection class.
template <Config Conf>
struct AndersonDirectionParams {
    bool rescale_on_step_size_changes = false;
};

/// @ingroup grp_DirectionProviders
template <Config Conf>
struct AndersonDirection {
    USING_ALPAQA_CONFIG(Conf);

    using Problem           = TypeErasedProblem<config_t>;
    using AndersonAccel     = alpaqa::AndersonAccel<config_t>;
    using AcceleratorParams = typename AndersonAccel::Params;
    using DirectionParams   = AndersonDirectionParams<config_t>;

    mutable AndersonAccel anderson;

    struct Params {
        AcceleratorParams accelerator = {};
        DirectionParams direction     = {};
    };

    AndersonDirection() = default;
    AndersonDirection(const Params &params)
        : anderson(params.accelerator), direction_params(params.direction) {}
    AndersonDirection(const typename AndersonAccel::Params &params,
                      const DirectionParams &directionparams = {})
        : anderson(params), direction_params(directionparams) {}
    AndersonDirection(const AndersonAccel &anderson,
                      const DirectionParams &directionparams = {})
        : anderson(anderson), direction_params(directionparams) {}
    AndersonDirection(AndersonAccel &&anderson,
                      const DirectionParams &directionparams = {})
        : anderson(std::move(anderson)), direction_params(directionparams) {}

    /// @see @ref PANOCDirection::initialize
    void initialize(const Problem &problem, [[maybe_unused]] crvec y,
                    [[maybe_unused]] crvec Σ, [[maybe_unused]] real_t γ_0,
                    [[maybe_unused]] crvec x_0, [[maybe_unused]] crvec x̂_0,
                    [[maybe_unused]] crvec p_0,
                    [[maybe_unused]] crvec grad_ψx_0) {
        anderson.resize(problem.get_n());
        anderson.initialize(x̂_0, p_0);
    }

    /// @see @ref PANOCDirection::has_initial_direction
    bool has_initial_direction() const { return false; }

    /// @see @ref PANOCDirection::update
    bool update([[maybe_unused]] real_t γₖ, [[maybe_unused]] real_t γₙₑₓₜ,
                [[maybe_unused]] crvec xₖ, [[maybe_unused]] crvec xₙₑₓₜ,
                [[maybe_unused]] crvec pₖ, [[maybe_unused]] crvec pₙₑₓₜ,
                [[maybe_unused]] crvec grad_ψxₖ,
                [[maybe_unused]] crvec grad_ψxₙₑₓₜ) {
        return true;
    }

    /// @see @ref PANOCDirection::apply
    bool apply([[maybe_unused]] real_t γₖ, [[maybe_unused]] crvec xₖ,
               [[maybe_unused]] crvec x̂ₖ, crvec pₖ,
               [[maybe_unused]] crvec grad_ψxₖ, rvec qₖ) const {
        anderson.compute(x̂ₖ, pₖ, qₖ);
        qₖ -= xₖ;
        return true;
    }

    /// @see @ref PANOCDirection::changed_γ
    void changed_γ(real_t γₖ, real_t old_γₖ) {
        if (direction_params.rescale_on_step_size_changes)
            anderson.scale_R(γₖ / old_γₖ);
        else
            anderson.reset();
    }

    /// @see @ref PANOCDirection::reset
    void reset() { anderson.reset(); }

    /// @see @ref PANOCDirection::get_name
    std::string get_name() const {
        return "AndersonDirection<" + std::string(config_t::get_name()) + '>';
    }
    auto get_params() const {
        return std::tie(anderson.get_params(), direction_params);
    }

    DirectionParams direction_params;
};

} // namespace alpaqa