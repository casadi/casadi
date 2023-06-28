#pragma once

#include <alpaqa/problem/ocproblem.hpp>

namespace alpaqa {

template <Config Conf>
void ControlProblemVTable<Conf>::default_get_D_N(const void *self, Box &D,
                                                 const ControlProblemVTable &vtable) {
    vtable.get_D(self, D, vtable);
}
template <Config Conf>
void ControlProblemVTable<Conf>::default_eval_add_Q_N(const void *self, crvec x, crvec h, rmat Q,
                                                      const ControlProblemVTable &vtable) {
    vtable.eval_add_Q(self, vtable.N, x, h, Q);
}
template <Config Conf>
void ControlProblemVTable<Conf>::default_eval_add_R_prod_masked(const void *, index_t, crvec, crvec,
                                                                crindexvec, crindexvec, crvec, rvec,
                                                                rvec,
                                                                const ControlProblemVTable &) {
    throw not_implemented_error("default_eval_add_R_prod_masked");
}
template <Config Conf>
void ControlProblemVTable<Conf>::default_eval_add_S_prod_masked(const void *, index_t, crvec, crvec,
                                                                crindexvec, crvec, rvec, rvec,
                                                                const ControlProblemVTable &) {
    throw not_implemented_error("default_eval_add_S_prod_masked");
}
template <Config Conf>
auto ControlProblemVTable<Conf>::default_get_R_work_size(const void *, const ControlProblemVTable &)
    -> length_t {
    return 0;
}
template <Config Conf>
auto ControlProblemVTable<Conf>::default_get_S_work_size(const void *, const ControlProblemVTable &)
    -> length_t {
    return 0;
}
template <Config Conf>
void ControlProblemVTable<Conf>::default_eval_constr_N(const void *self, crvec x, rvec c,
                                                       const ControlProblemVTable &vtable) {
    vtable.eval_constr(self, vtable.N, x, c, vtable);
}
template <Config Conf>
void ControlProblemVTable<Conf>::default_eval_grad_constr_prod_N(
    const void *self, crvec x, crvec p, rvec grad_cx_p, const ControlProblemVTable &vtable) {
    vtable.eval_grad_constr_prod(self, vtable.N, x, p, grad_cx_p, vtable);
}
template <Config Conf>
void ControlProblemVTable<Conf>::default_eval_add_gn_hess_constr_N(
    const void *self, crvec x, crvec M, rmat out, const ControlProblemVTable &vtable) {
    vtable.eval_add_gn_hess_constr(self, vtable.N, x, M, out, vtable);
}

} // namespace alpaqa