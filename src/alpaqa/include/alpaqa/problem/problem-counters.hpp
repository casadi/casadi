#pragma once

#include <alpaqa/export.h>

#include <chrono>
#include <iosfwd>

namespace alpaqa {

struct EvalCounter {
    unsigned proj_diff_g{};
    unsigned proj_multipliers{};
    unsigned prox_grad_step{};
    unsigned inactive_indices_res_lna{};
    unsigned f{};
    unsigned grad_f{};
    unsigned f_grad_f{};
    unsigned f_g{};
    unsigned grad_f_grad_g_prod{};
    unsigned g{};
    unsigned grad_g_prod{};
    unsigned grad_gi{};
    unsigned jac_g{};
    unsigned grad_L{};
    unsigned hess_L_prod{};
    unsigned hess_L{};
    unsigned hess_ψ_prod{};
    unsigned hess_ψ{};
    unsigned ψ{};
    unsigned grad_ψ{};
    unsigned ψ_grad_ψ{};

    struct EvalTimer {
        std::chrono::nanoseconds proj_diff_g{};
        std::chrono::nanoseconds proj_multipliers{};
        std::chrono::nanoseconds prox_grad_step{};
        std::chrono::nanoseconds inactive_indices_res_lna{};
        std::chrono::nanoseconds f{};
        std::chrono::nanoseconds grad_f{};
        std::chrono::nanoseconds f_grad_f{};
        std::chrono::nanoseconds f_g{};
        std::chrono::nanoseconds grad_f_grad_g_prod{};
        std::chrono::nanoseconds g{};
        std::chrono::nanoseconds grad_g_prod{};
        std::chrono::nanoseconds grad_gi{};
        std::chrono::nanoseconds jac_g{};
        std::chrono::nanoseconds grad_L{};
        std::chrono::nanoseconds hess_L_prod{};
        std::chrono::nanoseconds hess_L{};
        std::chrono::nanoseconds hess_ψ_prod{};
        std::chrono::nanoseconds hess_ψ{};
        std::chrono::nanoseconds ψ{};
        std::chrono::nanoseconds grad_ψ{};
        std::chrono::nanoseconds ψ_grad_ψ{};
    } time;

    void reset() { *this = {}; }
};

ALPAQA_EXPORT std::ostream &operator<<(std::ostream &, const EvalCounter &);

inline EvalCounter::EvalTimer &operator+=(EvalCounter::EvalTimer &a,
                                          const EvalCounter::EvalTimer &b) {
    a.proj_diff_g += b.proj_diff_g;
    a.proj_multipliers += b.proj_multipliers;
    a.prox_grad_step += b.prox_grad_step;
    a.inactive_indices_res_lna += b.inactive_indices_res_lna;
    a.f += b.f;
    a.grad_f += b.grad_f;
    a.f_grad_f += b.f_grad_f;
    a.f_g += b.f_g;
    a.grad_f_grad_g_prod += b.grad_f_grad_g_prod;
    a.g += b.g;
    a.grad_g_prod += b.grad_g_prod;
    a.grad_gi += b.grad_gi;
    a.jac_g += b.jac_g;
    a.grad_L += b.grad_L;
    a.hess_L_prod += b.hess_L_prod;
    a.hess_L += b.hess_L;
    a.hess_ψ_prod += b.hess_ψ_prod;
    a.hess_ψ += b.hess_ψ;
    a.ψ += b.ψ;
    a.grad_ψ += b.grad_ψ;
    a.ψ_grad_ψ += b.ψ_grad_ψ;
    return a;
}

inline EvalCounter &operator+=(EvalCounter &a, const EvalCounter &b) {
    a.proj_diff_g += b.proj_diff_g;
    a.proj_multipliers += b.proj_multipliers;
    a.prox_grad_step += b.prox_grad_step;
    a.inactive_indices_res_lna += b.inactive_indices_res_lna;
    a.f += b.f;
    a.grad_f += b.grad_f;
    a.f_grad_f += b.f_grad_f;
    a.f_g += b.f_g;
    a.grad_f_grad_g_prod += b.grad_f_grad_g_prod;
    a.g += b.g;
    a.grad_g_prod += b.grad_g_prod;
    a.grad_gi += b.grad_gi;
    a.jac_g += b.jac_g;
    a.grad_L += b.grad_L;
    a.hess_L_prod += b.hess_L_prod;
    a.hess_L += b.hess_L;
    a.hess_ψ_prod += b.hess_ψ_prod;
    a.hess_ψ += b.hess_ψ;
    a.ψ += b.ψ;
    a.grad_ψ += b.grad_ψ;
    a.ψ_grad_ψ += b.ψ_grad_ψ;
    a.time += b.time;
    return a;
}

inline EvalCounter operator+(EvalCounter a, const EvalCounter &b) { return a += b; }

} // namespace alpaqa
