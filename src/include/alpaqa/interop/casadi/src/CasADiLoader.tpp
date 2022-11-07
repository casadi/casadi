#pragma once

#include <alpaqa/interop/casadi/CasADiFunctionWrapper.hpp>
#include <alpaqa/interop/casadi/CasADiLoader.hpp>
#include <alpaqa/util/not-implemented.hpp>

#include <casadi/core/external.hpp>

#include <memory>
#include <optional>
#include <stdexcept>
#include <type_traits>

namespace alpaqa {

namespace casadi_loader {

template <class F>
auto wrap_load(const std::string &so_name, const char *name, F f) {
    try {
        return f();
    } catch (const std::invalid_argument &e) {
        throw std::invalid_argument("Unable to load function '" + so_name +
                                    ":" + name + "': " + e.what());
    }
}

template <class T, class... Args>
auto wrapped_load(const std::string &so_name, const char *name,
                  Args &&...args) {
    return wrap_load(so_name, name, [&] {
        return T(casadi::external(name, so_name), std::forward<Args>(args)...);
    });
}

using dim = std::pair<casadi_int, casadi_int>;
inline constexpr auto dims(auto... a) {
    if constexpr ((... && std::is_constructible_v<dim, decltype(a)>))
        return std::array{a...};
    else
        return std::array{dim{a, 1}...};
}

template <Config Conf>
struct CasADiFunctionsWithParam {
    static constexpr bool WithParam = true;
    CasADiFunctionEvaluator<Conf, 1 + WithParam, 1> f;
    // CasADiFunctionEvaluator<5 + WithParam, 1> grad_ψ;
    CasADiFunctionEvaluator<Conf, 5 + WithParam, 2> ψ_grad_ψ;
    struct ConstrFun {
        CasADiFunctionEvaluator<Conf, 1 + WithParam, 1> g;
        CasADiFunctionEvaluator<Conf, 2 + WithParam, 1> grad_L;
        CasADiFunctionEvaluator<Conf, 5 + WithParam, 2> ψ;
    };
    std::optional<ConstrFun> constr;
    struct HessFun {
        CasADiFunctionEvaluator<Conf, 3 + WithParam, 1> hess_L_prod;
        CasADiFunctionEvaluator<Conf, 2 + WithParam, 1> hess_L;
    };
    std::optional<HessFun> hess;
};

template <Config Conf>
struct CasADiControlFunctionsWithParam {
    static constexpr bool WithParam = true;
    CasADiFunctionEvaluator<Conf, 2 + WithParam, 1> f;
    CasADiFunctionEvaluator<Conf, 2 + WithParam, 1> jac_f;
    CasADiFunctionEvaluator<Conf, 1 + WithParam, 1> l;
    CasADiFunctionEvaluator<Conf, 1 + WithParam, 1> l_N;
    CasADiFunctionEvaluator<Conf, 1 + WithParam, 1> grad_l;
    CasADiFunctionEvaluator<Conf, 1 + WithParam, 1> grad_l_N;
    CasADiFunctionEvaluator<Conf, 1 + WithParam, 1> hess_l;
    CasADiFunctionEvaluator<Conf, 1 + WithParam, 1> hess_l_N;
};

} // namespace casadi_loader

template <Config Conf>
CasADiProblem<Conf>::CasADiProblem(const std::string &so_name, length_t n,
                                   length_t m, length_t p, bool second_order)
    : BoxConstrProblem<Conf>{n, m} {
    using namespace casadi_loader;
    auto load_g_unknown_dims =
        [&]() -> std::optional<CasADiFunctionEvaluator<Conf, 2, 1>> {
        casadi::Function gfun = casadi::external("g", so_name);
        using namespace std::literals::string_literals;
        if (gfun.n_in() != 2)
            throw std::invalid_argument(
                "Invalid number of input arguments: got "s +
                std::to_string(gfun.n_in()) + ", should be 2.");
        if (gfun.n_out() > 1)
            throw std::invalid_argument(
                "Invalid number of output arguments: got "s +
                std::to_string(gfun.n_in()) + ", should be 0 or 1.");
        if (gfun.size2_in(0) != 1)
            throw std::invalid_argument(
                "First input argument should be a column vector.");
        if (gfun.size2_in(1) != 1)
            throw std::invalid_argument(
                "Second input argument should be a column vector.");
        if (gfun.n_out() == 1 && gfun.size2_out(0) != 1)
            throw std::invalid_argument(
                "First output argument should be a column vector.");
        if (n <= 0)
            n = static_cast<length_t>(gfun.size1_in(0));
        if (m <= 0 && gfun.n_out() == 1)
            m = static_cast<length_t>(gfun.size1_out(0));
        if (p <= 0)
            p = static_cast<length_t>(gfun.size1_in(1));
        if (gfun.n_out() == 0) {
            if (m != 0)
                throw std::invalid_argument(
                    "Function g has no outputs but m != 0");
            return std::nullopt;
        }
        CasADiFunctionEvaluator<Conf, 2, 1> g{std::move(gfun)};
        g.validate_dimensions({dim(n, 1), dim(p, 1)}, {dim(m, 1)});
        return std::make_optional(std::move(g));
    };

    auto load_g_known_dims = [&] {
        CasADiFunctionEvaluator<Conf, 2, 1> g{casadi::external("g", so_name),
                                              {dim(n, 1), dim(p, 1)},
                                              {dim(m, 1)}};
        return g;
    };

    std::optional<CasADiFunctionEvaluator<Conf, 2, 1>> g =
        (n <= 0 || m <= 0 || p <= 0)
            // If not all dimensions are specified, load the function "g" to
            // determine the missing dimensions.
            ? wrap_load(so_name, "g", load_g_unknown_dims)
            // Otherwise, load the function "g" and compare its dimensions to
            // the dimensions specified by the user.
            : wrap_load(so_name, "g", load_g_known_dims);

    this->n     = n;
    this->m     = m;
    this->param = vec::Constant(p, alpaqa::NaN<Conf>);
    this->C     = {vec::Constant(n, +alpaqa::inf<Conf>),
                   vec::Constant(n, -alpaqa::inf<Conf>)};
    this->D     = {vec::Constant(m, +alpaqa::inf<Conf>),
                   vec::Constant(m, -alpaqa::inf<Conf>)};

    impl = std::make_unique<CasADiFunctionsWithParam<Conf>>(
        CasADiFunctionsWithParam<Conf>{
            wrapped_load<CasADiFunctionEvaluator<Conf, 2, 1>>( //
                so_name, "f", dims(n, p), dims(1)),
            // wrapped_load<CasADiFunctionEvaluator<6, 1>>( //
            //     so_name, "grad_psi", dims(n, p, m, m, m, m), dims(n)),
            wrapped_load<CasADiFunctionEvaluator<Conf, 6, 2>>( //
                so_name, "psi_grad_psi", dims(n, p, m, m, m, m), dims(1, n)),
            std::nullopt,
            std::nullopt,
        });

    if (g)
        impl->constr = std::make_optional(
            typename CasADiFunctionsWithParam<Conf>::ConstrFun{
                std::move(*g),
                wrapped_load<CasADiFunctionEvaluator<Conf, 3, 1>>( //
                    so_name, "grad_L", dims(n, p, m), dims(n)),
                wrapped_load<CasADiFunctionEvaluator<Conf, 6, 2>>( //
                    so_name, "psi", dims(n, p, m, m, m, m), dims(1, m)),
            });

    if (second_order)
        impl->hess =
            std::make_optional(typename CasADiFunctionsWithParam<Conf>::HessFun{
                wrapped_load<CasADiFunctionEvaluator<Conf, 4, 1>>( //
                    so_name, "hess_L_prod", dims(n, p, m, n), dims(n)),
                wrapped_load<CasADiFunctionEvaluator<Conf, 3, 1>>( //
                    so_name, "hess_L", dims(n, p, m), dims(dim(n, n))),
            });
}

template <Config Conf>
CasADiProblem<Conf>::CasADiProblem(const CasADiProblem &o)
    : BoxConstrProblem<config_t>{o}, param{o.param},
      impl(o.impl ? std::make_unique<Functions>(*o.impl) : nullptr) {}
// Note: impl can be null if the CasADiProblem was moved from
template <Config Conf>
CasADiProblem<Conf> &CasADiProblem<Conf>::operator=(const CasADiProblem &o) {
    BoxConstrProblem<config_t>::operator=(o);
    this->param = o.param;
    this->impl  = o.impl ? std::make_unique<Functions>(*o.impl) : nullptr;
    return *this;
}

template <Config Conf>
CasADiProblem<Conf>::CasADiProblem(CasADiProblem &&) = default;
template <Config Conf>
CasADiProblem<Conf> &CasADiProblem<Conf>::operator=(CasADiProblem &&) = default;

template <Config Conf>
CasADiProblem<Conf>::~CasADiProblem() = default;

template <Config Conf>
auto CasADiProblem<Conf>::eval_f(crvec x) const -> real_t {
    real_t f;
    impl->f({x.data(), param.data()}, {&f});
    return f;
}

template <Config Conf>
void CasADiProblem<Conf>::eval_grad_f(crvec, rvec) const {
    throw not_implemented_error("CasADiProblem::eval_grad_f"); // TODO
}

template <Config Conf>
void CasADiProblem<Conf>::eval_g(crvec x, rvec g) const {
    if (impl->constr)
        impl->constr->g({x.data(), param.data()}, {g.data()});
    else
        throw std::logic_error("No constraints function g"); // TODO
}

template <Config Conf>
void CasADiProblem<Conf>::eval_grad_g_prod(crvec, crvec, rvec) const {
    throw not_implemented_error("CasADiProblem::eval_grad_g_prod"); // TODO
}

template <Config Conf>
void CasADiProblem<Conf>::eval_grad_ψ(crvec x, crvec y, crvec Σ, rvec grad_ψ,
                                      rvec, rvec) const {
#if 0
    impl->grad_ψ({x.data(), param.data(), y.data(), Σ.data(),
                  this->D.lowerbound.data(), this->D.upperbound.data()},
                 {grad_ψ.data()});
#else
    // This seems to be faster than having a specialized function. Possibly
    // cache-related?
    real_t ψ;
    impl->ψ_grad_ψ({x.data(), param.data(), y.data(), Σ.data(),
                    this->D.lowerbound.data(), this->D.upperbound.data()},
                   {&ψ, grad_ψ.data()});
#endif
}

template <Config Conf>
typename CasADiProblem<Conf>::real_t
CasADiProblem<Conf>::eval_ψ_grad_ψ(crvec x, crvec y, crvec Σ, rvec grad_ψ, rvec,
                                   rvec) const {
    real_t ψ;
    impl->ψ_grad_ψ({x.data(), param.data(), y.data(), Σ.data(),
                    this->D.lowerbound.data(), this->D.upperbound.data()},
                   {&ψ, grad_ψ.data()});
    return ψ;
}

template <Config Conf>
void CasADiProblem<Conf>::eval_grad_L(crvec x, crvec y, rvec grad_L,
                                      rvec) const {
    if (impl->constr)
        impl->constr->grad_L({x.data(), param.data(), y.data()},
                             {grad_L.data()});
    else
        throw std::logic_error("No function grad_L");
}

template <Config Conf>
typename CasADiProblem<Conf>::real_t
CasADiProblem<Conf>::eval_ψ(crvec x, crvec y, crvec Σ, rvec ŷ) const {
    real_t ψ;
    if (impl->constr)
        impl->constr->ψ({x.data(), param.data(), y.data(), Σ.data(),
                         this->D.lowerbound.data(), this->D.upperbound.data()},
                        {&ψ, ŷ.data()});
    else
        impl->f({x.data(), param.data()}, {&ψ});
    return ψ;
}

template <Config Conf>
void CasADiProblem<Conf>::eval_grad_ψ_from_ŷ(crvec x, crvec ŷ, rvec grad_ψ,
                                             rvec) const {
    if (this->m == 0) {
        real_t ψ;
        impl->ψ_grad_ψ(
            {x.data(), param.data(), nullptr, nullptr, nullptr, nullptr},
            {&ψ, grad_ψ.data()});
    } else {
        impl->constr->grad_L({x.data(), param.data(), ŷ.data()},
                             {grad_ψ.data()});
    }
}

template <Config Conf>
void CasADiProblem<Conf>::eval_grad_gi(crvec, index_t, rvec) const {
    throw not_implemented_error("CasADiProblem::eval_grad_gi"); // TODO
}

template <Config Conf>
void CasADiProblem<Conf>::eval_hess_L_prod(crvec x, crvec y, crvec v,
                                           rvec Hv) const {
    impl->hess->hess_L_prod({x.data(), param.data(), y.data(), v.data()},
                            {Hv.data()});
}
template <Config Conf>
void CasADiProblem<Conf>::eval_hess_L(crvec x, crvec y, rmat H) const {
    impl->hess->hess_L({x.data(), param.data(), y.data()}, {H.data()});
}

template <Config Conf>
bool CasADiProblem<Conf>::provides_eval_grad_gi() const {
    return false; // TODO
}
template <Config Conf>
bool CasADiProblem<Conf>::provides_eval_hess_L_prod() const {
    return impl->hess.has_value();
}
template <Config Conf>
bool CasADiProblem<Conf>::provides_eval_hess_L() const {
    return impl->hess.has_value();
}

template <Config Conf>
auto load_CasADi_problem_with_param(const std::string &filename,
                                    typename Conf::length_t n,
                                    typename Conf::length_t m,
                                    typename Conf::length_t p,
                                    bool second_order) {
    return CasADiProblem<Conf>{filename, n, m, p, second_order};
}

// -------------------------------------------------------------------------- //

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
    this->U      = {vec::Constant(nu, +alpaqa::inf<Conf>),
                    vec::Constant(nu, -alpaqa::inf<Conf>)};

    impl = std::make_unique<CasADiControlFunctionsWithParam<Conf>>(
        CasADiControlFunctionsWithParam<Conf>{
            .f     = std::move(f),
            .jac_f = wrapped_load<CasADiFunctionEvaluator<Conf, 3, 1>>(
                so_name, "jac_f", dims(nx, nu, p), dims(dim(nx, nx + nu))),
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
CasADiControlProblem<Conf>::CasADiControlProblem(CasADiControlProblem &&) =
    default;
template <Config Conf>
CasADiControlProblem<Conf> &
CasADiControlProblem<Conf>::operator=(CasADiControlProblem &&) = default;

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