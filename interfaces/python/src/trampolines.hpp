#pragma once

#include <alpaqa/config/config.hpp>

template <class ProblemBase>
class ProblemTrampoline : ProblemBase {
    USING_ALPAQA_CONFIG_TEMPLATE(ProblemBase::config_t);
    using ProblemBase::ProblemBase;

    // clang-format off
    real_t eval_f(crvec x) const override { PYBIND11_OVERRIDE(real_t, ProblemBase, eval_f, x); }
    void eval_grad_f(crvec x, rvec grad_fx) const override { PYBIND11_OVERRIDE(void, ProblemBase, eval_grad_f, x, grad_fx); }
    void eval_g(crvec x, rvec gx) const override { PYBIND11_OVERRIDE(void, ProblemBase, eval_g, x, gx); }
    void eval_grad_g_prod(crvec x, crvec y, rvec grad_gxy) const override { PYBIND11_OVERRIDE(void, ProblemBase, eval_grad_g_prod, x, y, grad_gxy); }
    void eval_grad_gi(crvec x, index_t i, rvec grad_gi) const override { PYBIND11_OVERRIDE(void, ProblemBase, eval_grad_gi, x, i, grad_gi); }
    void eval_hess_L_prod(crvec x, crvec y, crvec v, rvec Hv) const override { PYBIND11_OVERRIDE(void, ProblemBase, eval_hess_L_prod, x, y, v, Hv); }
    void eval_hess_L(crvec x, crvec y, rmat H) const override { PYBIND11_OVERRIDE(void, ProblemBase, eval_hess_L, x, y, H); }
    // clang-format on
};
