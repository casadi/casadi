#pragma once

#include <alpaqa/problem/type-erased-problem.hpp>

namespace alpaqa {

template <Config Conf>
auto ProblemVTable<Conf>::calc_ŷ_dᵀŷ(const void *self, rvec g_ŷ, crvec y, crvec Σ,
                                     const ProblemVTable &vtable) -> real_t {
    if (Σ.size() == 1) {
        // ζ = g(x) + Σ⁻¹y
        g_ŷ += (1 / Σ(0)) * y;
        // d = ζ - Π(ζ, D)
        vtable.eval_proj_diff_g(self, g_ŷ, g_ŷ);
        // dᵀŷ, ŷ = Σ d
        real_t dᵀŷ = Σ(0) * g_ŷ.dot(g_ŷ);
        g_ŷ *= Σ(0);
        return dᵀŷ;
    } else {
        // ζ = g(x) + Σ⁻¹y
        g_ŷ += Σ.asDiagonal().inverse() * y;
        // d = ζ - Π(ζ, D)
        vtable.eval_proj_diff_g(self, g_ŷ, g_ŷ);
        // dᵀŷ, ŷ = Σ d
        real_t dᵀŷ = 0;
        for (index_t i = 0; i < y.size(); ++i) {
            dᵀŷ += g_ŷ(i) * Σ(i) * g_ŷ(i); // TODO: vectorize
            g_ŷ(i) = Σ(i) * g_ŷ(i);
        }
        return dᵀŷ;
    }
}

template <Config Conf>
auto ProblemVTable<Conf>::default_eval_inactive_indices_res_lna(const void *, real_t, crvec, crvec,
                                                                rindexvec, const ProblemVTable &)
    -> index_t {
    throw not_implemented_error("eval_inactive_indices_res_lna");
}

template <Config Conf>
void ProblemVTable<Conf>::default_eval_jac_g(const void *, crvec, rindexvec, rindexvec, rvec,
                                             const ProblemVTable &vtable) {
    if (vtable.m != 0)
        throw not_implemented_error("eval_jac_g");
}

template <Config Conf>
auto ProblemVTable<Conf>::default_get_jac_g_num_nonzeros(const void *, const ProblemVTable &)
    -> length_t {
    return -1;
}

template <Config Conf>
void ProblemVTable<Conf>::default_eval_grad_gi(const void *, crvec, index_t, rvec,
                                               const ProblemVTable &) {
    throw not_implemented_error("eval_grad_gi");
}

template <Config Conf>
void ProblemVTable<Conf>::default_eval_hess_L_prod(const void *, crvec, crvec, real_t, crvec, rvec,
                                                   const ProblemVTable &) {
    throw not_implemented_error("eval_hess_L_prod");
}

template <Config Conf>
void ProblemVTable<Conf>::default_eval_hess_L(const void *, crvec, crvec, real_t, rindexvec,
                                              rindexvec, rvec, const ProblemVTable &) {
    throw not_implemented_error("eval_hess_L");
}

template <Config Conf>
auto ProblemVTable<Conf>::default_get_hess_L_num_nonzeros(const void *, const ProblemVTable &)
    -> length_t {
    return -1;
}

template <Config Conf>
void ProblemVTable<Conf>::default_eval_hess_ψ_prod(const void *self, crvec x, crvec y, crvec,
                                                   real_t scale, crvec v, rvec Hv,
                                                   const ProblemVTable &vtable) {
    if (y.size() == 0 && vtable.eval_hess_L_prod != ProblemVTable<Conf>::default_eval_hess_L_prod)
        return vtable.eval_hess_L_prod(self, x, y, scale, v, Hv, vtable);
    throw not_implemented_error("eval_hess_ψ_prod");
}

template <Config Conf>
void ProblemVTable<Conf>::default_eval_hess_ψ(const void *self, crvec x, crvec y, crvec,
                                              real_t scale, rindexvec inner_idx,
                                              rindexvec outer_ptr, rvec H_values,
                                              const ProblemVTable &vtable) {
    if (y.size() == 0 && vtable.eval_hess_L != default_eval_hess_L)
        return vtable.eval_hess_L(self, x, y, scale, inner_idx, outer_ptr, H_values, vtable);
    throw not_implemented_error("eval_hess_ψ");
}

template <Config Conf>
auto ProblemVTable<Conf>::default_get_hess_ψ_num_nonzeros(const void *, const ProblemVTable &)
    -> length_t {
    return 0;
}

/** @implementation{ProblemVTable<Conf>::default_eval_f_grad_f} */
template <Config Conf>
/* [ProblemVTable<Conf>::default_eval_f_grad_f] */
auto ProblemVTable<Conf>::default_eval_f_grad_f(const void *self, crvec x, rvec grad_fx,
                                                const ProblemVTable &vtable) -> real_t {
    vtable.eval_grad_f(self, x, grad_fx);
    return vtable.eval_f(self, x);
}
/* [ProblemVTable<Conf>::default_eval_f_grad_f] */

/** @implementation{ProblemVTable<Conf>::default_eval_f_g} */
template <Config Conf>
/* [ProblemVTable<Conf>::default_eval_f_g] */
auto ProblemVTable<Conf>::default_eval_f_g(const void *self, crvec x, rvec g,
                                           const ProblemVTable &vtable) -> real_t {
    vtable.eval_g(self, x, g);
    return vtable.eval_f(self, x);
}
/* [ProblemVTable<Conf>::default_eval_f_g] */

/** @implementation{ProblemVTable<Conf>::default_eval_grad_f_grad_g_prod} */
template <Config Conf>
/* [ProblemVTable<Conf>::default_eval_grad_f_grad_g_prod] */
void ProblemVTable<Conf>::default_eval_grad_f_grad_g_prod(const void *self, crvec x, crvec y,
                                                          rvec grad_f, rvec grad_gxy,
                                                          const ProblemVTable &vtable) {
    vtable.eval_grad_f(self, x, grad_f);
    vtable.eval_grad_g_prod(self, x, y, grad_gxy);
}
/* [ProblemVTable<Conf>::default_eval_grad_f_grad_g_prod] */

/** @implementation{ProblemVTable<Conf>::default_eval_grad_L} */
template <Config Conf>
/* [ProblemVTable<Conf>::default_eval_grad_L] */
void ProblemVTable<Conf>::default_eval_grad_L(const void *self, crvec x, crvec y, rvec grad_L,
                                              rvec work_n, const ProblemVTable &vtable) {
    if (y.size() == 0) /* [[unlikely]] */
        return vtable.eval_grad_f(self, x, grad_L);
    vtable.eval_grad_f_grad_g_prod(self, x, y, grad_L, work_n, vtable);
    grad_L += work_n;
}
/* [ProblemVTable<Conf>::default_eval_grad_L] */

/** @implementation{ProblemVTable<Conf>::default_eval_ψ} */
template <Config Conf>
/* [ProblemVTable<Conf>::default_eval_ψ] */
auto ProblemVTable<Conf>::default_eval_ψ(const void *self, crvec x, crvec y, crvec Σ, rvec ŷ,
                                         const ProblemVTable &vtable) -> real_t {
    if (y.size() == 0) /* [[unlikely]] */
        return vtable.eval_f(self, x);

    auto f   = vtable.eval_f_g(self, x, ŷ, vtable);
    auto dᵀŷ = calc_ŷ_dᵀŷ(self, ŷ, y, Σ, vtable);
    // ψ(x) = f(x) + ½ dᵀŷ
    auto ψ = f + real_t(0.5) * dᵀŷ;
    return ψ;
}
/* [ProblemVTable<Conf>::default_eval_ψ] */

/** @implementation{ProblemVTable<Conf>::default_eval_grad_ψ} */
template <Config Conf>
/* [ProblemVTable<Conf>::default_eval_grad_ψ] */
void ProblemVTable<Conf>::default_eval_grad_ψ(const void *self, crvec x, crvec y, crvec Σ,
                                              rvec grad_ψ, rvec work_n, rvec work_m,
                                              const ProblemVTable &vtable) {
    if (y.size() == 0) /* [[unlikely]] */ {
        vtable.eval_grad_f(self, x, grad_ψ);
    } else {
        vtable.eval_g(self, x, work_m);
        (void)calc_ŷ_dᵀŷ(self, work_m, y, Σ, vtable);
        vtable.eval_grad_L(self, x, work_m, grad_ψ, work_n, vtable);
    }
}
/* [ProblemVTable<Conf>::default_eval_grad_ψ] */

/** @implementation{ProblemVTable<Conf>::default_eval_ψ_grad_ψ} */
template <Config Conf>
/* [ProblemVTable<Conf>::default_eval_ψ_grad_ψ] */
auto ProblemVTable<Conf>::default_eval_ψ_grad_ψ(const void *self, crvec x, crvec y, crvec Σ,
                                                rvec grad_ψ, rvec work_n, rvec work_m,
                                                const ProblemVTable &vtable) -> real_t {
    if (y.size() == 0) /* [[unlikely]] */
        return vtable.eval_f_grad_f(self, x, grad_ψ, vtable);

    auto &ŷ = work_m;
    // ψ(x) = f(x) + ½ dᵀŷ
    auto f   = vtable.eval_f_g(self, x, ŷ, vtable);
    auto dᵀŷ = calc_ŷ_dᵀŷ(self, ŷ, y, Σ, vtable);
    auto ψ   = f + real_t(0.5) * dᵀŷ;
    // ∇ψ(x) = ∇f(x) + ∇g(x) ŷ
    vtable.eval_grad_L(self, x, ŷ, grad_ψ, work_n, vtable);
    return ψ;
}
/* [ProblemVTable<Conf>::default_eval_ψ_grad_ψ] */

template <Config Conf>
auto ProblemVTable<Conf>::default_get_box_C(const void *, const ProblemVTable &) -> const Box & {
    throw not_implemented_error("get_box_C");
}

template <Config Conf>
auto ProblemVTable<Conf>::default_get_box_D(const void *, const ProblemVTable &) -> const Box & {
    throw not_implemented_error("get_box_D");
}

template <Config Conf>
void ProblemVTable<Conf>::default_check(const void *, const ProblemVTable &) {}

} // namespace alpaqa