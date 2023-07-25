#pragma once

#include <alpaqa/casadi/CasADiFunctionWrapper.hpp>
#include <alpaqa/casadi/CasADiProblem.hpp>
#include <alpaqa/util/io/csv.hpp>
#include <alpaqa/util/not-implemented.hpp>
#include "CasADiLoader-util.hpp"

#include <casadi/core/external.hpp>

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <memory>
#include <optional>
#include <stdexcept>
#include <type_traits>

namespace alpaqa {

namespace fs = std::filesystem;

namespace casadi_loader {

template <Config Conf>
struct CasADiFunctionsWithParam {
    USING_ALPAQA_CONFIG(Conf);
    CasADiFunctionEvaluator<Conf, 2, 1> f;
    CasADiFunctionEvaluator<Conf, 2, 2> f_grad_f;
    // CasADiFunctionEvaluator<6, 1> grad_ψ;
    CasADiFunctionEvaluator<Conf, 6, 2> ψ_grad_ψ;
    struct ConstrFun {
        CasADiFunctionEvaluator<Conf, 2, 1> g;
        CasADiFunctionEvaluator<Conf, 3, 1> grad_L;
        CasADiFunctionEvaluator<Conf, 6, 2> ψ;
    };
    std::optional<ConstrFun> constr = std::nullopt;
    std::optional<CasADiFunctionEvaluator<Conf, 5, 1>> hess_L_prod =
        std::nullopt;
    std::optional<CasADiFunctionEvaluator<Conf, 4, 1>> hess_L = std::nullopt;
    std::optional<CasADiFunctionEvaluator<Conf, 8, 1>> hess_ψ_prod =
        std::nullopt;
    std::optional<CasADiFunctionEvaluator<Conf, 7, 1>> hess_ψ = std::nullopt;
    std::optional<CasADiFunctionEvaluator<Conf, 2, 1>> jac_g  = std::nullopt;
    indexvec inner_H{}, outer_H{};
    indexvec inner_J{}, outer_J{};
};

} // namespace casadi_loader

namespace detail {

template <Config Conf>
auto casadi_to_index(casadi_int i) -> index_t<Conf> {
    return static_cast<index_t<Conf>>(i);
}
} // namespace detail

template <Config Conf>
CasADiProblem<Conf>::CasADiProblem(const std::string &so_name)
    : BoxConstrProblem<Conf>{0, 0} {
    using namespace casadi_loader;
    length_t n = 0, m = 0, p = 0;
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
        n = static_cast<length_t>(gfun.size1_in(0));
        if (gfun.n_out() == 1)
            m = static_cast<length_t>(gfun.size1_out(0));
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

    auto g = wrap_load(so_name, "g", load_g_unknown_dims);

    this->n     = n;
    this->m     = m;
    this->param = vec::Constant(p, alpaqa::NaN<Conf>);
    this->C     = Box<config_t>{n};
    this->D     = Box<config_t>{m};

    impl = std::make_unique<CasADiFunctionsWithParam<Conf>>(
        CasADiFunctionsWithParam<Conf>{
            .f        = wrapped_load<CasADiFunctionEvaluator<Conf, 2, 1>>( //
                so_name, "f", dims(n, p), dims(1)),
            .f_grad_f = wrapped_load<CasADiFunctionEvaluator<Conf, 2, 2>>( //
                so_name, "f_grad_f", dims(n, p), dims(1, n)),
            // .grad_ψ = wrapped_load<CasADiFunctionEvaluator<6, 1>>( //
            //     so_name, "grad_psi", dims(n, p, m, m, m, m), dims(n)),
            .ψ_grad_ψ = wrapped_load<CasADiFunctionEvaluator<Conf, 6, 2>>( //
                so_name, "psi_grad_psi", dims(n, p, m, m, m, m), dims(1, n)),
        });

    if (g)
        impl->constr = {
            std::move(*g),
            wrapped_load<CasADiFunctionEvaluator<Conf, 3, 1>>( //
                so_name, "grad_L", dims(n, p, m), dims(n)),
            wrapped_load<CasADiFunctionEvaluator<Conf, 6, 2>>( //
                so_name, "psi", dims(n, p, m, m, m, m), dims(1, m)),
        };

    impl->hess_L_prod = try_load<CasADiFunctionEvaluator<Conf, 5, 1>>( //
        so_name, "hess_L_prod", dims(n, p, m, 1, n), dims(n));
    impl->hess_L      = try_load<CasADiFunctionEvaluator<Conf, 4, 1>>( //
        so_name, "hess_L", dims(n, p, m, 1), dims(dim(n, n)));
    impl->hess_ψ_prod = try_load<CasADiFunctionEvaluator<Conf, 8, 1>>( //
        so_name, "hess_psi_prod", dims(n, p, m, m, 1, m, m, n), dims(n));
    impl->hess_ψ      = try_load<CasADiFunctionEvaluator<Conf, 7, 1>>( //
        so_name, "hess_psi", dims(n, p, m, m, 1, m, m), dims(dim(n, n)));
    impl->jac_g       = try_load<CasADiFunctionEvaluator<Conf, 2, 1>>( //
        so_name, "jacobian_g", dims(n, p), dims(dim(m, n)));

    auto data_filepath = fs::path{so_name}.replace_extension("csv");
    if (fs::exists(data_filepath))
        load_numerical_data(data_filepath);
}

template <Config Conf>
void CasADiProblem<Conf>::load_numerical_data(
    const std::filesystem::path &filepath, char sep) {
    // Open data file
    std::ifstream data_file{filepath};
    if (!data_file)
        throw std::runtime_error("Unable to open data file \"" +
                                 filepath.string() + '"');

    // Helper function for reading single line of (float) data
    index_t line        = 0;
    auto wrap_data_load = [&](std::string_view name, auto &v, bool fixed_size) {
        try {
            ++line;
            if (data_file.peek() == '\n') // Ignore empty lines
                return static_cast<void>(data_file.get());
            if (fixed_size) {
                csv::read_row(data_file, v, sep);
            } else { // Dynamic size
                auto s = csv::read_row_std_vector<real_t>(data_file, sep);
                v      = cmvec{s.data(), static_cast<index_t>(s.size())};
            }
        } catch (csv::read_error &e) {
            // Transform any errors in something more readable
            throw std::runtime_error("Unable to read " + std::string(name) +
                                     " from data file \"" + filepath.string() +
                                     ':' + std::to_string(line) +
                                     "\": " + e.what());
        }
    };
    // Helper function for reading a single value
    auto read_single = [&](std::string_view name, auto &v) {
        data_file >> v;
        if (!data_file)
            throw std::runtime_error("Unable to read " + std::string(name) +
                                     " from data file \"" + filepath.string() +
                                     ':' + std::to_string(line) + '"');
    };
    // Read the bounds, parameter value, and regularization
    wrap_data_load("C.lowerbound", this->C.lowerbound, true);
    wrap_data_load("C.upperbound", this->C.upperbound, true);
    wrap_data_load("D.lowerbound", this->D.lowerbound, true);
    wrap_data_load("D.upperbound", this->D.upperbound, true);
    wrap_data_load("param", this->param, true);
    wrap_data_load("l1_reg", this->l1_reg, false);
    // Penalty/ALM split is a single integer
    read_single("penalty_alm_split", this->penalty_alm_split);
}

template <Config Conf>
CasADiProblem<Conf>::CasADiProblem(const CasADiProblem &) = default;
template <Config Conf>
CasADiProblem<Conf> &
CasADiProblem<Conf>::operator=(const CasADiProblem &) = default;
template <Config Conf>
CasADiProblem<Conf>::CasADiProblem(CasADiProblem &&) noexcept = default;
template <Config Conf>
CasADiProblem<Conf> &
CasADiProblem<Conf>::operator=(CasADiProblem &&) noexcept = default;

template <Config Conf>
CasADiProblem<Conf>::~CasADiProblem() = default;

template <Config Conf>
auto CasADiProblem<Conf>::eval_f(crvec x) const -> real_t {
    real_t f;
    impl->f({x.data(), param.data()}, {&f});
    return f;
}

template <Config Conf>
void CasADiProblem<Conf>::eval_grad_f(crvec x, rvec grad_fx) const {
    real_t f;
    impl->f_grad_f({x.data(), param.data()}, {&f, grad_fx.data()});
}

template <Config Conf>
auto CasADiProblem<Conf>::eval_f_grad_f(crvec x, rvec grad_fx) const -> real_t {
    real_t f;
    impl->f_grad_f({x.data(), param.data()}, {&f, grad_fx.data()});
    return f;
}

template <Config Conf>
void CasADiProblem<Conf>::eval_g(crvec x, rvec g) const {
    if (impl->constr)
        impl->constr->g({x.data(), param.data()}, {g.data()});
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
        eval_f_grad_f(x, grad_L);
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
void CasADiProblem<Conf>::eval_grad_gi(crvec, index_t, rvec) const {
    throw not_implemented_error("CasADiProblem::eval_grad_gi"); // TODO
}

template <Config Conf>
auto CasADiProblem<Conf>::get_jac_g_num_nonzeros() const -> length_t {
    if (!impl->jac_g.has_value())
        return -1;
    auto &&sparsity = impl->jac_g->fun.sparsity_out(0);
    return sparsity.is_dense() ? -1 : static_cast<length_t>(sparsity.nnz());
}

template <Config Conf>
void CasADiProblem<Conf>::eval_jac_g(crvec x, rindexvec inner_idx,
                                     rindexvec outer_ptr, rvec J_values) const {
    assert(impl->jac_g.has_value());
    if (J_values.size() > 0) {
        (*impl->jac_g)({x.data(), param.data()}, {J_values.data()});
    } else {
        auto &&sparsity = impl->jac_g->fun.sparsity_out(0);
        using detail::casadi_to_index;
        if (!sparsity.is_dense()) {
            std::transform(sparsity.row(), sparsity.row() + sparsity.nnz(),
                           inner_idx.begin(), casadi_to_index<config_t>);
            std::transform(sparsity.colind(),
                           sparsity.colind() + this->get_n() + 1,
                           outer_ptr.begin(), casadi_to_index<config_t>);
        }
    }
}

template <Config Conf>
void CasADiProblem<Conf>::eval_hess_L_prod(crvec x, crvec y, real_t scale,
                                           crvec v, rvec Hv) const {
    assert(impl->hess_L_prod.has_value());
    (*impl->hess_L_prod)({x.data(), param.data(), y.data(), &scale, v.data()},
                         {Hv.data()});
}

template <Config Conf>
auto CasADiProblem<Conf>::get_hess_L_num_nonzeros() const -> length_t {
    if (!impl->hess_L.has_value())
        return -1;
    auto &&sparsity = impl->hess_L->fun.sparsity_out(0);
    return sparsity.is_dense() ? -1 : static_cast<length_t>(sparsity.nnz());
}

template <Config Conf>
void CasADiProblem<Conf>::eval_hess_L(crvec x, crvec y, real_t scale,
                                      rindexvec inner_idx, rindexvec outer_ptr,
                                      rvec H_values) const {
    assert(impl->hess_L.has_value());
    if (H_values.size() > 0) {
        (*impl->hess_L)({x.data(), param.data(), y.data(), &scale},
                        {H_values.data()});
    } else {
        auto &&sparsity = impl->hess_L->fun.sparsity_out(0);
        using detail::casadi_to_index;
        if (!sparsity.is_dense()) {
            std::transform(sparsity.row(), sparsity.row() + sparsity.nnz(),
                           inner_idx.begin(), casadi_to_index<config_t>);
            std::transform(sparsity.colind(),
                           sparsity.colind() + this->get_n() + 1,
                           outer_ptr.begin(), casadi_to_index<config_t>);
        }
    }
}

template <Config Conf>
void CasADiProblem<Conf>::eval_hess_ψ_prod(crvec x, crvec y, crvec Σ,
                                           real_t scale, crvec v,
                                           rvec Hv) const {
    assert(impl->hess_ψ_prod.has_value());
    (*impl->hess_ψ_prod)({x.data(), param.data(), y.data(), Σ.data(), &scale,
                          this->D.lowerbound.data(), this->D.upperbound.data(),
                          v.data()},
                         {Hv.data()});
}

template <Config Conf>
auto CasADiProblem<Conf>::get_hess_ψ_num_nonzeros() const -> length_t {
    if (!impl->hess_ψ.has_value())
        return 0;
    auto &&sparsity = impl->hess_ψ->fun.sparsity_out(0);
    return sparsity.is_dense() ? 0 : static_cast<length_t>(sparsity.nnz());
}

template <Config Conf>
void CasADiProblem<Conf>::eval_hess_ψ(crvec x, crvec y, crvec Σ, real_t scale,
                                      rindexvec inner_idx, rindexvec outer_ptr,
                                      rvec H_values) const {
    assert(impl->hess_ψ.has_value());
    if (H_values.size() > 0) {
        (*impl->hess_ψ)({x.data(), param.data(), y.data(), Σ.data(), &scale,
                         this->D.lowerbound.data(), this->D.upperbound.data()},
                        {H_values.data()});
    } else {
        auto &&sparsity = impl->hess_ψ->fun.sparsity_out(0);
        using detail::casadi_to_index;
        if (!sparsity.is_dense()) {
            std::transform(sparsity.row(), sparsity.row() + sparsity.nnz(),
                           inner_idx.begin(), casadi_to_index<config_t>);
            std::transform(sparsity.colind(),
                           sparsity.colind() + this->get_n() + 1,
                           outer_ptr.begin(), casadi_to_index<config_t>);
        }
    }
}

template <Config Conf>
bool CasADiProblem<Conf>::provides_eval_grad_gi() const {
    return false; // TODO
}
template <Config Conf>
bool CasADiProblem<Conf>::provides_eval_jac_g() const {
    return impl->jac_g.has_value();
}
template <Config Conf>
bool CasADiProblem<Conf>::provides_eval_hess_L_prod() const {
    return impl->hess_L_prod.has_value();
}
template <Config Conf>
bool CasADiProblem<Conf>::provides_eval_hess_L() const {
    return impl->hess_L.has_value();
}
template <Config Conf>
bool CasADiProblem<Conf>::provides_eval_hess_ψ_prod() const {
    return impl->hess_ψ_prod.has_value();
}
template <Config Conf>
bool CasADiProblem<Conf>::provides_eval_hess_ψ() const {
    return impl->hess_ψ.has_value();
}

} // namespace alpaqa