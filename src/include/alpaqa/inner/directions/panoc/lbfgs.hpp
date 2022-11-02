#pragma once

#include <alpaqa/accelerators/lbfgs.hpp>
#include <alpaqa/inner/directions/panoc-direction-update.hpp>
#include <alpaqa/problem/type-erased-problem.hpp>

namespace alpaqa {

/// @ingroup grp_DirectionProviders
template <Config Conf>
struct PANOCDirection<LBFGS<Conf>> {
    USING_ALPAQA_CONFIG(Conf);

    using Problem = TypeErasedProblem<Conf>;
    using LBFGS   = alpaqa::LBFGS<Conf>;
    using Params  = typename LBFGS::Params;

    LBFGS lbfgs;

    struct ExtraParams {
        bool rescale_when_γ_changes = false;
    } extraparams;

    PANOCDirection(const Params &params, const ExtraParams &extraparams = {})
        : lbfgs(params), extraparams(extraparams) {}
    PANOCDirection(const LBFGS &lbfgs, const ExtraParams &extraparams = {})
        : lbfgs(lbfgs), extraparams(extraparams) {}
    PANOCDirection(LBFGS &&lbfgs, const ExtraParams &extraparams = {})
        : lbfgs(std::move(lbfgs)), extraparams(extraparams) {}

    void initialize(const Problem &problem, [[maybe_unused]] crvec y,
                    [[maybe_unused]] crvec Σ, [[maybe_unused]] real_t γ_0,
                    [[maybe_unused]] crvec x_0, [[maybe_unused]] crvec x̂_0,
                    [[maybe_unused]] crvec p_0,
                    [[maybe_unused]] crvec grad_ψx_0) {
        lbfgs.resize(problem.get_n());
    }

    bool update([[maybe_unused]] real_t γₖ, [[maybe_unused]] real_t γₙₑₓₜ,
                crvec xₖ, crvec xₙₑₓₜ, crvec pₖ, crvec pₙₑₓₜ,
                [[maybe_unused]] crvec grad_ψxₖ,
                [[maybe_unused]] crvec grad_ψxₙₑₓₜ) {
        return lbfgs.update(xₖ, xₙₑₓₜ, pₖ, pₙₑₓₜ, LBFGS::Sign::Negative);
    }

    bool apply([[maybe_unused]] real_t γₖ, [[maybe_unused]] crvec xₖ,
               [[maybe_unused]] crvec x̂ₖ, crvec pₖ,
               [[maybe_unused]] crvec grad_ψxₖ, rvec qₖ) const {
        qₖ = pₖ;
        return lbfgs.apply(qₖ, γₖ);
    }

    void changed_γ(real_t γₖ, real_t old_γₖ) {
        if (extraparams.rescale_when_γ_changes)
            lbfgs.scale_y(γₖ / old_γₖ);
        else
            lbfgs.reset();
    }

    void reset() { lbfgs.reset(); }

    std::string get_name() const { return lbfgs.get_name(); }
    Params get_params() const { return lbfgs.get_params(); }
};

} // namespace alpaqa