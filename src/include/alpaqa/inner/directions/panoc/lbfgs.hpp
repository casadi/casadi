#pragma once

#include <alpaqa/accelerators/lbfgs.hpp>
#include <alpaqa/inner/directions/panoc-direction-update.hpp>

namespace alpaqa {

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

    void initialize(crvec x₀, crvec x̂₀, crvec p₀, crvec grad₀) {
        lbfgs.resize(x₀.size());
        (void)x̂₀;
        (void)p₀;
        (void)grad₀;
    }

    bool update(crvec xₖ, crvec xₖ₊₁, crvec pₖ, crvec pₖ₊₁, crvec grad_new,
                const Box<config_t> &C, real_t γ_new) {
        (void)grad_new;
        (void)C;
        (void)γ_new;
        return lbfgs.update(xₖ, xₖ₊₁, pₖ, pₖ₊₁, LBFGS::Sign::Negative);
    }

    bool apply(crvec xₖ, crvec x̂ₖ, crvec pₖ, real_t γ, rvec qₖ) const {
        (void)xₖ;
        (void)x̂ₖ;
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