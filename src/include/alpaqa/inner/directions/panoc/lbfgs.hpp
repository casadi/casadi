#pragma once

#include <alpaqa/accelerators/lbfgs.hpp>
#include <alpaqa/inner/directions/panoc-direction-update.hpp>

namespace alpaqa {

/// @ingroup grp_DirectionProviders
template <Config Conf>
struct PANOCDirection<LBFGS<Conf>> {
    USING_ALPAQA_CONFIG(Conf);

    using LBFGS  = alpaqa::LBFGS<Conf>;
    using Params = typename LBFGS::Params;

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

    void initialize(crvec x_0, crvec x̂_0, crvec p_0, crvec grad_0) {
        lbfgs.resize(x_0.size());
        (void)x̂_0;
        (void)p_0;
        (void)grad_0;
    }

    bool update(crvec xₖ, crvec xₙₑₓₜ, crvec pₖ, crvec pₙₑₓₜ, crvec grad_new,
                const Box<config_t> &C, real_t γ_new) {
        (void)grad_new;
        (void)C;
        (void)γ_new;
        return lbfgs.update(xₖ, xₙₑₓₜ, pₖ, pₙₑₓₜ, LBFGS::Sign::Negative);
    }

    bool apply(crvec xₖ, crvec x̂ₖ, crvec pₖ, crvec grad_xₖ, crvec grad_x̂ₖ,
               real_t γ, rvec qₖ) const {
        (void)xₖ;
        (void)x̂ₖ;
        (void)grad_xₖ;
        (void)grad_x̂ₖ;
        qₖ = pₖ;
        return lbfgs.apply(qₖ, γ);
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