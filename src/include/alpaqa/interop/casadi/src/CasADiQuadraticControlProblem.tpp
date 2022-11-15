#pragma once

#include <alpaqa/interop/casadi/CasADiFunctionWrapper.hpp>
#include <alpaqa/interop/casadi/CasADiQuadraticControlProblem.hpp>
#include <alpaqa/util/not-implemented.hpp>
#include "CasADiLoader-util.hpp"

#include <Eigen/src/Core/CwiseTernaryOp.h>
#include <casadi/core/external.hpp>

#include <memory>
#include <optional>
#include <stdexcept>
#include <type_traits>

namespace alpaqa {

namespace casadi_loader {

template <Config Conf>
struct CasADiQuadraticControlFunctionsWithParam {
    static constexpr bool WithParam = true;
    CasADiFunctionEvaluator<Conf, 2 + WithParam, 1> f;
    CasADiFunctionEvaluator<Conf, 2 + WithParam, 1> jac_f;
    CasADiFunctionEvaluator<Conf, 3 + WithParam, 1> grad_f_prod;
};

} // namespace casadi_loader

template <Config Conf>
CasADiQuadraticControlProblem<Conf>::CasADiQuadraticControlProblem(
    const std::string &filename, length_t N, length_t nx, length_t nu,
    length_t p)
    : N{N} {
    using namespace casadi_loader;
    auto load_f_unknown_dims = [&]() -> CasADiFunctionEvaluator<Conf, 3, 1> {
        casadi::Function ffun = casadi::external("f", filename);
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
            casadi::external("f", filename),
            {dim(nx, 1), dim(nu, 1), dim(p, 1)},
            {dim(nx, 1)}};
        return f;
    };

    CasADiFunctionEvaluator<Conf, 3, 1> f =
        (nx <= 0 || nu <= 0 || p <= 0)
            // If not all dimensions are specified, load the function "f" to
            // determine the missing dimensions.
            ? wrap_load(filename, "f", load_f_unknown_dims)
            // Otherwise, load the function "f" and compare its dimensions to
            // the dimensions specified by the user.
            : wrap_load(filename, "f", load_f_known_dims);

    this->nx     = nx;
    this->nu     = nu;
    this->x_init = vec::Constant(nx, alpaqa::NaN<Conf>);
    this->param  = vec::Constant(p, alpaqa::NaN<Conf>);
    this->U      = {vec::Constant(nu, +alpaqa::inf<Conf>),
                    vec::Constant(nu, -alpaqa::inf<Conf>)};
    this->D      = {vec::Constant(nx, +alpaqa::inf<Conf>),
                    vec::Constant(nx, -alpaqa::inf<Conf>)};
    this->D_N    = {vec::Constant(nx, +alpaqa::inf<Conf>),
                    vec::Constant(nx, -alpaqa::inf<Conf>)};
    this->Q      = vec::Constant(nx, alpaqa::NaN<Conf>);
    this->R      = vec::Constant(nu, alpaqa::NaN<Conf>);
    this->Q_N    = vec::Constant(nx, alpaqa::NaN<Conf>);
    this->x_ref  = vec::Constant(nx, alpaqa::NaN<Conf>);
    this->u_ref  = vec::Constant(nu, alpaqa::NaN<Conf>);
    this->μ      = mat::Constant(nx, N + 1, alpaqa::NaN<Conf>);

    impl = std::make_unique<CasADiQuadraticControlFunctionsWithParam<Conf>>(
        CasADiQuadraticControlFunctionsWithParam<Conf>{
            .f     = std::move(f),
            .jac_f = wrapped_load<CasADiFunctionEvaluator<Conf, 3, 1>>(
                filename, "jac_f", dims(nx, nu, p), dims(dim(nx, nx + nu))),
            .grad_f_prod = wrapped_load<CasADiFunctionEvaluator<Conf, 4, 1>>(
                filename, "grad_f_prod", dims(nx, nu, p, nx), dims(nx + nu)),
        });
}

template <Config Conf>
CasADiQuadraticControlProblem<Conf>::CasADiQuadraticControlProblem(
    const CasADiQuadraticControlProblem &o) = default;
template <Config Conf>
CasADiQuadraticControlProblem<Conf> &
CasADiQuadraticControlProblem<Conf>::operator=(
    const CasADiQuadraticControlProblem &o) = default;

template <Config Conf>
CasADiQuadraticControlProblem<Conf>::CasADiQuadraticControlProblem(
    CasADiQuadraticControlProblem &&) noexcept = default;
template <Config Conf>
CasADiQuadraticControlProblem<Conf> &
CasADiQuadraticControlProblem<Conf>::operator=(
    CasADiQuadraticControlProblem &&) noexcept = default;

template <Config Conf>
CasADiQuadraticControlProblem<Conf>::~CasADiQuadraticControlProblem() = default;

template <Config Conf>
void CasADiQuadraticControlProblem<Conf>::eval_f(index_t, crvec x, crvec u,
                                                 rvec fxu) const {
    impl->f({x.data(), u.data(), param.data()}, {fxu.data()});
}
template <Config Conf>
void CasADiQuadraticControlProblem<Conf>::eval_jac_f(index_t, crvec x, crvec u,
                                                     rmat J_fxu) const {
    impl->jac_f({x.data(), u.data(), param.data()}, {J_fxu.data()});
}
template <Config Conf>
void CasADiQuadraticControlProblem<Conf>::eval_grad_f_prod(
    index_t, crvec x, crvec u, crvec v, rvec grad_fxu_p) const {
    impl->grad_f_prod({x.data(), u.data(), param.data(), v.data()},
                      {grad_fxu_p.data()});
}
template <Config Conf>
auto CasADiQuadraticControlProblem<Conf>::eval_l(index_t t, crvec h) const
    -> real_t {
    auto x       = h.segment(0, nx);
    auto u       = h.segment(nx, nu);
    auto x_ref_t = x_ref.col(t >= x_ref.cols() ? 0 : t);
    auto u_ref_t = u_ref.col(t >= u_ref.cols() ? 0 : t);
    real_t lx    = (x - x_ref_t).cwiseAbs2().cwiseProduct(Q).sum();
    real_t lu    = (u - u_ref_t).cwiseAbs2().cwiseProduct(R).sum();
    auto x_proj  = D.lowerbound.cwiseMax(x).cwiseMin(D.upperbound);
    real_t lz    = μ(t) * (x - x_proj).squaredNorm();
    return real_t{0.5} * (lx + lu + lz);
}
template <Config Conf>
auto CasADiQuadraticControlProblem<Conf>::eval_l_N(crvec h) const -> real_t {
    auto x_ref_t = x_ref.col(N >= x_ref.cols() ? 0 : N);
    real_t lx    = (h - x_ref_t).cwiseAbs2().cwiseProduct(Q_N).sum();
    auto h_proj  = D_N.lowerbound.cwiseMax(h).cwiseMin(D_N.upperbound);
    real_t lz    = μ(N) * (h - h_proj).squaredNorm();
    return real_t{0.5} * (lx + lz);
}
template <Config Conf>
void dist_sq_hess(crvec<Conf> z, const Box<Conf> &D, crvec<Conf> μ,
                  rvec<Conf> hess) {
    for (index_t<Conf> i = 0; i < z.size(); ++i)
        hess(i) = μ(i) * (z(i) <= D.lowerbound(i) || z(i) >= D.upperbound(i));
}
template <Config Conf>
void CasADiQuadraticControlProblem<Conf>::eval_grad_l(index_t t, crvec h,
                                                      rvec grad_lh) const {
    auto x                  = h.segment(0, nx);
    auto u                  = h.segment(nx, nu);
    auto x_ref_t            = x_ref.col(t >= x_ref.cols() ? 0 : t);
    auto u_ref_t            = u_ref.col(t >= u_ref.cols() ? 0 : t);
    grad_lh.segment(0, nx)  = Q.cwiseProduct(x - x_ref_t);
    grad_lh.segment(nx, nu) = R.cwiseProduct(u - u_ref_t);
    auto x_proj             = D.lowerbound.cwiseMax(x).cwiseMin(D.upperbound);
    grad_lh.segment(0, nx) += μ.col(t).cwiseProduct(x - x_proj);
}
template <Config Conf>
void CasADiQuadraticControlProblem<Conf>::eval_grad_l_N(crvec h,
                                                        rvec grad_lh) const {
    auto x_ref_t = x_ref.col(N >= x_ref.cols() ? 0 : N);
    auto h_proj  = D_N.lowerbound.cwiseMax(h).cwiseMin(D_N.upperbound);
    grad_lh      = Q_N.cwiseProduct(h - x_ref_t);
    grad_lh += μ.col(N).cwiseProduct(h - h_proj);
}
template <Config Conf>
void CasADiQuadraticControlProblem<Conf>::eval_hess_l(index_t t, crvec h,
                                                      rmat hess_lh) const {
    auto x = h.segment(0, nx);
    dist_sq_hess(x, D, μ.col(t), hess_lh.col(0).segment(0, nx));
    hess_lh.col(0).segment(0, nx) += Q;
    hess_lh.col(0).segment(nx, nu) = R;
}
template <Config Conf>
void CasADiQuadraticControlProblem<Conf>::eval_hess_l_N(crvec h,
                                                        rmat hess_lh) const {
    dist_sq_hess(h, D_N, μ.col(N), hess_lh.col(0));
    hess_lh.col(0) += Q_N;
}

} // namespace alpaqa