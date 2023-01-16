#pragma once

#include <alpaqa/interop/casadi/CasADiControlProblem.hpp>
#include <alpaqa/interop/casadi/CasADiFunctionWrapper.hpp>
#include <alpaqa/util/not-implemented.hpp>
#include "CasADiLoader-util.hpp"

#include <casadi/core/external.hpp>

#include <memory>
#include <optional>
#include <stdexcept>
#include <type_traits>

namespace alpaqa {

namespace casadi_loader {

template <Config Conf>
struct CasADiControlFunctionsWithParam {
    static constexpr bool WithParam = true;
    CasADiFunctionEvaluator<Conf, 2 + WithParam, 1> f;
    CasADiFunctionEvaluator<Conf, 2 + WithParam, 1> jac_f;
    CasADiFunctionEvaluator<Conf, 3 + WithParam, 1> grad_f_prod;
    CasADiFunctionEvaluator<Conf, 1 + WithParam, 1> l;
    CasADiFunctionEvaluator<Conf, 1 + WithParam, 1> l_N;
    CasADiFunctionEvaluator<Conf, 1 + WithParam, 1> grad_l;
    CasADiFunctionEvaluator<Conf, 1 + WithParam, 1> grad_l_N;
    CasADiFunctionEvaluator<Conf, 1 + WithParam, 1> hess_l;
    CasADiFunctionEvaluator<Conf, 1 + WithParam, 1> hess_l_N;
};

} // namespace casadi_loader

template <Config Conf>
CasADiControlProblem<Conf>::CasADiControlProblem(const std::string &so_name,
                                                 length_t N, length_t nx,
                                                 length_t nu, length_t p)
    : N{N} {
    using namespace casadi_loader;
    auto load_f_unknown_dims = [&]() -> CasADiFunctionEvaluator<Conf, 3, 1> {
        casadi::Function ffun = casadi::external("f", so_name);
        using namespace std::literals::string_literals;
        if (ffun.n_in() != 3)
            throw std::invalid_argument(
                "Invalid number of input arguments: got "s +
                std::to_string(ffun.n_in()) + ", should be 3.");
        if (ffun.n_out() != 1)
            throw std::invalid_argument(
                "Invalid number of output arguments: got "s +
                std::to_string(ffun.n_in()) + ", should be 1.");
        if (ffun.size2_in(0) != 1)
            throw std::invalid_argument(
                "First input argument should be a column vector.");
        if (ffun.size2_in(1) != 1)
            throw std::invalid_argument(
                "Second input argument should be a column vector.");
        if (ffun.size2_in(2) != 1)
            throw std::invalid_argument(
                "Third input argument should be a column vector.");
        if (ffun.n_out() == 1 && ffun.size2_out(0) != 1)
            throw std::invalid_argument(
                "First output argument should be a column vector.");
        if (nx <= 0)
            nx = static_cast<length_t>(ffun.size1_in(0));
        if (nu <= 0)
            nu = static_cast<length_t>(ffun.size1_in(1));
        if (p <= 0)
            p = static_cast<length_t>(ffun.size1_in(2));
        CasADiFunctionEvaluator<Conf, 3, 1> f{std::move(ffun)};
        f.validate_dimensions({dim(nx, 1), dim(nu, 1), dim(p, 1)},
                              {dim(nx, 1)});
        return f;
    };

    auto load_f_known_dims = [&] {
        CasADiFunctionEvaluator<Conf, 3, 1> f{
            casadi::external("f", so_name),
            {dim(nx, 1), dim(nu, 1), dim(p, 1)},
            {dim(nx, 1)}};
        return f;
    };

    CasADiFunctionEvaluator<Conf, 3, 1> f =
        (nx <= 0 || nu <= 0 || p <= 0)
            // If not all dimensions are specified, load the function "f" to
            // determine the missing dimensions.
            ? wrap_load(so_name, "f", load_f_unknown_dims)
            // Otherwise, load the function "f" and compare its dimensions to
            // the dimensions specified by the user.
            : wrap_load(so_name, "f", load_f_known_dims);

    this->nx     = nx;
    this->nu     = nu;
    this->x_init = vec::Constant(nx, alpaqa::NaN<Conf>);
    this->param  = vec::Constant(p, alpaqa::NaN<Conf>);
    this->U      = Box{nu};

    impl = std::make_unique<CasADiControlFunctionsWithParam<Conf>>(
        CasADiControlFunctionsWithParam<Conf>{
            .f     = std::move(f),
            .jac_f = wrapped_load<CasADiFunctionEvaluator<Conf, 3, 1>>(
                so_name, "jacobian_f", dims(nx, nu, p), dims(dim(nx, nx + nu))),
            .grad_f_prod = wrapped_load<CasADiFunctionEvaluator<Conf, 4, 1>>(
                so_name, "grad_f_prod", dims(nx, nu, p, nx), dims(nx + nu)),
            .l = wrapped_load<CasADiFunctionEvaluator<Conf, 2, 1>>(
                so_name, "l", dims(nx + nu, p), dims(1)),
            .l_N = wrapped_load<CasADiFunctionEvaluator<Conf, 2, 1>>(
                so_name, "l_N", dims(nx, p), dims(1)),
            .grad_l = wrapped_load<CasADiFunctionEvaluator<Conf, 2, 1>>(
                so_name, "grad_l", dims(nx + nu, p), dims(nx + nu)),
            .grad_l_N = wrapped_load<CasADiFunctionEvaluator<Conf, 2, 1>>(
                so_name, "grad_l_N", dims(nx, p), dims(nx)),
            .hess_l = wrapped_load<CasADiFunctionEvaluator<Conf, 2, 1>>(
                so_name, "hess_l", dims(nx + nu, p),
                dims(dim{nx + nu, nx + nu})),
            .hess_l_N = wrapped_load<CasADiFunctionEvaluator<Conf, 2, 1>>(
                so_name, "hess_l_N", dims(nx, p), dims(dim{nx, nx})),
        });
}

template <Config Conf>
CasADiControlProblem<Conf>::CasADiControlProblem(
    const CasADiControlProblem &o) = default;
template <Config Conf>
CasADiControlProblem<Conf> &
CasADiControlProblem<Conf>::operator=(const CasADiControlProblem &o) = default;

template <Config Conf>
CasADiControlProblem<Conf>::CasADiControlProblem(
    CasADiControlProblem &&) noexcept = default;
template <Config Conf>
CasADiControlProblem<Conf> &CasADiControlProblem<Conf>::operator=(
    CasADiControlProblem &&) noexcept = default;

template <Config Conf>
CasADiControlProblem<Conf>::~CasADiControlProblem() = default;

template <Config Conf>
void CasADiControlProblem<Conf>::eval_f(index_t, crvec x, crvec u,
                                        rvec fxu) const {
    impl->f({x.data(), u.data(), param.data()}, {fxu.data()});
}
template <Config Conf>
void CasADiControlProblem<Conf>::eval_jac_f(index_t, crvec x, crvec u,
                                            rmat J_fxu) const {
    impl->jac_f({x.data(), u.data(), param.data()}, {J_fxu.data()});
}
template <Config Conf>
void CasADiControlProblem<Conf>::eval_grad_f_prod(index_t, crvec x, crvec u,
                                                  crvec v,
                                                  rvec grad_fxu_p) const {
    impl->grad_f_prod({x.data(), u.data(), param.data(), v.data()},
                      {grad_fxu_p.data()});
}
template <Config Conf>
auto CasADiControlProblem<Conf>::eval_l(index_t, crvec h) const -> real_t {
    real_t l;
    impl->l({h.data(), param.data()}, {&l});
    return l;
}
template <Config Conf>
auto CasADiControlProblem<Conf>::eval_l_N(crvec h) const -> real_t {
    real_t l;
    impl->l_N({h.data(), param.data()}, {&l});
    return l;
}
template <Config Conf>
void CasADiControlProblem<Conf>::eval_grad_l(index_t, crvec h,
                                             rvec grad_lh) const {
    impl->grad_l({h.data(), param.data()}, {grad_lh.data()});
}
template <Config Conf>
void CasADiControlProblem<Conf>::eval_grad_l_N(crvec h, rvec grad_lh) const {
    impl->grad_l_N({h.data(), param.data()}, {grad_lh.data()});
}
template <Config Conf>
void CasADiControlProblem<Conf>::eval_hess_l(index_t, crvec h,
                                             rmat hess_lh) const {
    impl->hess_l({h.data(), param.data()}, {hess_lh.data()});
}
template <Config Conf>
void CasADiControlProblem<Conf>::eval_hess_l_N(crvec h, rmat hess_lh) const {
    impl->hess_l_N({h.data(), param.data()}, {hess_lh.data()});
}

} // namespace alpaqa