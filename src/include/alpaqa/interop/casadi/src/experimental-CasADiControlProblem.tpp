#pragma once

#include <alpaqa/config/config.hpp>
#include <alpaqa/interop/casadi/CasADiFunctionWrapper.hpp>
#include <alpaqa/interop/casadi/experimental-CasADiControlProblem.hpp>
#include <alpaqa/util/not-implemented.hpp>
#include "CasADiLoader-util.hpp"
#include "alpaqa/util/iter-adapter.hpp"
#include "alpaqa/util/set-intersection.hpp"

#include <Eigen/Sparse>
#include <casadi/core/external.hpp>

#include <memory>
#include <optional>
#include <ranges>
#include <span>
#include <stdexcept>
#include <type_traits>

namespace alpaqa::experimental {

namespace casadi_loader {

using namespace alpaqa::casadi_loader;

template <Config Conf>
struct CasADiControlFunctionsWithParam {
    USING_ALPAQA_CONFIG(Conf);

    static constexpr bool WithParam = true;
    CasADiFunctionEvaluator<Conf, 2 + WithParam, 1> f;
    CasADiFunctionEvaluator<Conf, 2 + WithParam, 1> jac_f;
    CasADiFunctionEvaluator<Conf, 3 + WithParam, 1> grad_f_prod;
    CasADiFunctionEvaluator<Conf, 2 + WithParam, 1> h;
    CasADiFunctionEvaluator<Conf, 1 + WithParam, 1> h_N;
    CasADiFunctionEvaluator<Conf, 1 + WithParam, 1> l;
    CasADiFunctionEvaluator<Conf, 1 + WithParam, 1> l_N;

    CasADiFunctionEvaluator<Conf, 2 + WithParam, 1> qr;
    CasADiFunctionEvaluator<Conf, 2 + WithParam, 1> q_N;
    CasADiFunctionEvaluator<Conf, 2 + WithParam, 1> Q;
    CasADiFunctionEvaluator<Conf, 2 + WithParam, 1> Q_N;
    CasADiFunctionEvaluator<Conf, 2 + WithParam, 1> R;
    CasADiFunctionEvaluator<Conf, 2 + WithParam, 1> S;

    CasADiFunctionEvaluator<Conf, 1 + WithParam, 1> c;
    CasADiFunctionEvaluator<Conf, 2 + WithParam, 1> grad_c_prod;
    CasADiFunctionEvaluator<Conf, 2 + WithParam, 1> gn_hess_c;
};

} // namespace casadi_loader

template <Config Conf>
CasADiControlProblem<Conf>::CasADiControlProblem(const std::string &so_name,
                                                 length_t N)
    : N{N} {
    length_t p;
    using namespace casadi_loader;
    auto load_f = [&]() -> CasADiFunctionEvaluator<Conf, 3, 1> {
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
        // if (ffun.size2_in(0) != 1)
        //     throw std::invalid_argument(
        //         "First input argument should be a column vector.");
        // if (ffun.size2_in(1) != 1)
        //     throw std::invalid_argument(
        //         "Second input argument should be a column vector.");
        // if (ffun.size2_in(2) != 1)
        //     throw std::invalid_argument(
        //         "Third input argument should be a column vector.");
        // if (ffun.size2_out(0) != 1)
        //     throw std::invalid_argument(
        //         "First output argument should be a column vector.");
        nx = static_cast<length_t>(ffun.size1_in(0));
        nu = static_cast<length_t>(ffun.size1_in(1));
        p  = static_cast<length_t>(ffun.size1_in(2));
        CasADiFunctionEvaluator<Conf, 3, 1> f{std::move(ffun)};
        f.validate_dimensions({dim(nx, 1), dim(nu, 1), dim(p, 1)},
                              {dim(nx, 1)});
        return f;
    };
    auto load_h = [&]() -> CasADiFunctionEvaluator<Conf, 3, 1> {
        casadi::Function hfun = casadi::external("h", so_name);
        using namespace std::literals::string_literals;
        if (hfun.n_in() != 3)
            throw std::invalid_argument(
                "Invalid number of input arguments: got "s +
                std::to_string(hfun.n_in()) + ", should be 3.");
        if (hfun.n_out() != 1)
            throw std::invalid_argument(
                "Invalid number of output arguments: got "s +
                std::to_string(hfun.n_in()) + ", should be 1.");
        // if (hfun.n_out() == 1 && hfun.size2_out(0) != 1)
        //     throw std::invalid_argument(
        //         "First output argument should be a column vector.");
        nh = static_cast<length_t>(hfun.size1_out(0));
        CasADiFunctionEvaluator<Conf, 3, 1> h{std::move(hfun)};
        h.validate_dimensions({dim(nx, 1), dim(nu, 1), dim(p, 1)},
                              {dim(nh, 1)});
        return h;
    };
    auto load_h_N = [&]() -> CasADiFunctionEvaluator<Conf, 2, 1> {
        casadi::Function hfun = casadi::external("h_N", so_name);
        using namespace std::literals::string_literals;
        if (hfun.n_in() != 2)
            throw std::invalid_argument(
                "Invalid number of input arguments: got "s +
                std::to_string(hfun.n_in()) + ", should be 2.");
        if (hfun.n_out() != 1)
            throw std::invalid_argument(
                "Invalid number of output arguments: got "s +
                std::to_string(hfun.n_in()) + ", should be 1.");
        nh_N = static_cast<length_t>(hfun.size1_out(0));
        CasADiFunctionEvaluator<Conf, 2, 1> h{std::move(hfun)};
        h.validate_dimensions({dim(nx, 1), dim(p, 1)}, {dim(nh_N, 1)});
        return h;
    };
    auto load_c = [&]() -> CasADiFunctionEvaluator<Conf, 2, 1> {
        casadi::Function cfun = casadi::external("c", so_name);
        using namespace std::literals::string_literals;
        if (cfun.n_in() != 2)
            throw std::invalid_argument(
                "Invalid number of input arguments: got "s +
                std::to_string(cfun.n_in()) + ", should be 2.");
        if (cfun.n_out() != 1)
            throw std::invalid_argument(
                "Invalid number of output arguments: got "s +
                std::to_string(cfun.n_in()) + ", should be 1.");
        // if (cfun.n_out() == 1 && cfun.size2_out(0) != 1)
        //     throw std::invalid_argument(
        //         "First output argument should be a column vector.");
        nc = static_cast<length_t>(cfun.size1_out(0));
        CasADiFunctionEvaluator<Conf, 2, 1> c{std::move(cfun)};
        c.validate_dimensions({dim(nx, 1), dim(p, 1)}, {dim(nc, 1)});
        return c;
    };
    // Load the functions "f", "h", and "c" to determine the unknown dimensions.
    auto f   = wrap_load(so_name, "f", load_f);
    auto h   = wrap_load(so_name, "h", load_h);
    auto h_N = wrap_load(so_name, "h_N", load_h_N);
    auto c   = wrap_load(so_name, "c", load_c);

    this->x_init = vec::Constant(nx, alpaqa::NaN<Conf>);
    this->param  = vec::Constant(p, alpaqa::NaN<Conf>);
    this->U      = Box{nu};
    this->D      = Box{nc};
    this->D_N    = Box{nc};

    impl = std::make_unique<CasADiControlFunctionsWithParam<Conf>>(
        CasADiControlFunctionsWithParam<Conf>{
            .f     = std::move(f),
            .jac_f = wrapped_load<CasADiFunctionEvaluator<Conf, 3, 1>>(
                so_name, "jacobian_f", dims(nx, nu, p), dims(dim(nx, nx + nu))),
            .grad_f_prod = wrapped_load<CasADiFunctionEvaluator<Conf, 4, 1>>(
                so_name, "grad_f_prod", dims(nx, nu, p, nx), dims(nx + nu)),
            .h   = std::move(h),
            .h_N = std::move(h_N),
            .l   = wrapped_load<CasADiFunctionEvaluator<Conf, 2, 1>>(
                so_name, "l", dims(nh, p), dims(1)),
            .l_N = wrapped_load<CasADiFunctionEvaluator<Conf, 2, 1>>(
                so_name, "l_N", dims(nh_N, p), dims(1)),
            .qr = wrapped_load<CasADiFunctionEvaluator<Conf, 3, 1>>(
                so_name, "qr", dims(nx + nu, nh, p), dims(nx + nu)),
            .q_N = wrapped_load<CasADiFunctionEvaluator<Conf, 3, 1>>(
                so_name, "q_N", dims(nx, nh_N, p), dims(nx)),
            .Q = wrapped_load<CasADiFunctionEvaluator<Conf, 3, 1>>(
                so_name, "Q", dims(nx + nu, nh, p), dims(dim{nx, nx})),
            .Q_N = wrapped_load<CasADiFunctionEvaluator<Conf, 3, 1>>(
                so_name, "Q_N", dims(nx, nh_N, p), dims(dim{nx, nx})),
            .R = wrapped_load<CasADiFunctionEvaluator<Conf, 3, 1>>(
                so_name, "R", dims(nx + nu, nh, p), dims(dim{nu, nu})),
            .S = wrapped_load<CasADiFunctionEvaluator<Conf, 3, 1>>(
                so_name, "S", dims(nx + nu, nh, p), dims(dim{nu, nx})),
            .c           = std::move(c),
            .grad_c_prod = wrapped_load<CasADiFunctionEvaluator<Conf, 3, 1>>(
                so_name, "grad_c_prod", dims(nx, p, nc), dims(nx)),
            .gn_hess_c = wrapped_load<CasADiFunctionEvaluator<Conf, 3, 1>>(
                so_name, "gn_hess_c", dims(nx, p, nc), dims(dim{nx, nx})),
        });

    auto n_work = std::max({
        impl->Q.fun.sparsity_out(0).nnz(),
        impl->Q_N.fun.sparsity_out(0).nnz(),
        impl->gn_hess_c.fun.sparsity_out(0).nnz(),
    });
    this->work  = vec::Constant(static_cast<length_t>(n_work), NaN<Conf>);
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
    assert(x.size() == nx);
    assert(u.size() == nu);
    assert(fxu.size() == nx);
    impl->f({x.data(), u.data(), param.data()}, {fxu.data()});
}
template <Config Conf>
void CasADiControlProblem<Conf>::eval_jac_f(index_t, crvec x, crvec u,
                                            rmat J_fxu) const {
    assert(x.size() == nx);
    assert(u.size() == nu);
    assert(J_fxu.rows() == nx);
    assert(J_fxu.cols() == nx + nu);
    impl->jac_f({x.data(), u.data(), param.data()}, {J_fxu.data()});
}
template <Config Conf>
void CasADiControlProblem<Conf>::eval_grad_f_prod(index_t, crvec x, crvec u,
                                                  crvec p,
                                                  rvec grad_fxu_p) const {
    assert(x.size() == nx);
    assert(u.size() == nu);
    assert(p.size() == nx);
    assert(grad_fxu_p.size() == nx + nu);
    impl->grad_f_prod({x.data(), u.data(), param.data(), p.data()},
                      {grad_fxu_p.data()});
}
template <Config Conf>
void CasADiControlProblem<Conf>::eval_h(index_t, crvec x, crvec u,
                                        rvec h) const {
    assert(x.size() == nx);
    assert(u.size() == nu);
    assert(h.size() == nh);
    impl->h({x.data(), u.data(), param.data()}, {h.data()});
}
template <Config Conf>
void CasADiControlProblem<Conf>::eval_h_N(crvec x, rvec h) const {
    assert(x.size() == nx);
    assert(h.size() == nh);
    impl->h({x.data(), param.data()}, {h.data()});
}
template <Config Conf>
auto CasADiControlProblem<Conf>::eval_l(index_t, crvec h) const -> real_t {
    assert(h.size() == nh);
    real_t l;
    impl->l({h.data(), param.data()}, {&l});
    return l;
}
template <Config Conf>
auto CasADiControlProblem<Conf>::eval_l_N(crvec h) const -> real_t {
    assert(h.size() == nh);
    real_t l;
    impl->l_N({h.data(), param.data()}, {&l});
    return l;
}
template <Config Conf>
void CasADiControlProblem<Conf>::eval_qr(index_t, crvec xu, crvec h,
                                         rvec qr) const {
    assert(xu.size() == nx + nu);
    assert(h.size() == nh);
    assert(qr.size() == nx + nu);
    impl->qr({xu.data(), h.data(), param.data()}, {qr.data()});
}
template <Config Conf>
void CasADiControlProblem<Conf>::eval_q_N(crvec x, crvec h, rvec q) const {
    assert(x.size() == nx);
    assert(h.size() == nh);
    assert(q.size() == nx);
    impl->q_N({x.data(), h.data(), param.data()}, {q.data()});
}
template <Config Conf>
void CasADiControlProblem<Conf>::eval_add_Q(index_t, crvec xu, crvec h,
                                            rmat Q) const {
    assert(xu.size() == nx + nu);
    assert(h.size() == nh);
    assert(Q.rows() == nx);
    assert(Q.cols() == nx);
    impl->Q({xu.data(), h.data(), param.data()}, {work.data()});
    using spmat   = Eigen::SparseMatrix<real_t, Eigen::ColMajor, casadi_int>;
    using cmspmat = Eigen::Map<const spmat>;
    auto &&sparse = impl->Q.fun.sparsity_out(0);
    if (sparse.is_dense())
        Q += cmmat{work.data(), nx, nx};
    else
        Q += cmspmat{
            nx, nx, sparse.nnz(), sparse.colind(), sparse.row(), work.data(),
        };
}
template <Config Conf>
void CasADiControlProblem<Conf>::eval_add_Q_N(crvec x, crvec h, rmat Q) const {
    assert(x.size() == nx);
    assert(h.size() == nh);
    assert(Q.rows() == nx);
    assert(Q.cols() == nx);
    impl->Q_N({x.data(), h.data(), param.data()}, {work.data()});
    auto &&sparse = impl->Q_N.fun.sparsity_out(0);
    using spmat   = Eigen::SparseMatrix<real_t, Eigen::ColMajor, casadi_int>;
    using cmspmat = Eigen::Map<const spmat>;
    if (sparse.is_dense())
        Q += cmmat{work.data(), nx, nx};
    else
        Q += cmspmat{
            nx, nx, sparse.nnz(), sparse.colind(), sparse.row(), work.data(),
        };
}

namespace detail {
template <class SpMat, class MaskVec>
auto select_rows_in_col(const SpMat &sp_mat, MaskVec mask, auto column) {
    using row_iter_t = typename SpMat::InnerIterator;
    util::iter_range_adapter<row_iter_t> col_range{{sp_mat, column}};
    static constexpr auto proj_row = [](const row_iter_t &it) {
        return static_cast<typename MaskVec::value_type>(it.row());
    };
    std::span mask_span{mask.data(), static_cast<size_t>(mask.size())};
    auto intersection = util::iter_set_intersection(
        std::move(col_range), std::move(mask_span), std::less{}, proj_row);
    return std::views::transform(std::move(intersection),
                                 []<class T>(T &&tup) -> decltype(auto) {
                                     return std::get<0>(std::forward<T>(tup));
                                 });
}
template <class SpMat, class MaskVec>
auto select_rows_in_col_iota(const SpMat &sp_mat, MaskVec mask, auto column) {
    using row_iter_t = typename SpMat::InnerIterator;
    util::iter_range_adapter<row_iter_t> col_range{{sp_mat, column}};
    std::span mask_span{mask.data(), static_cast<size_t>(mask.size())};
    static constexpr auto proj_row = [](const row_iter_t &it) {
        return static_cast<typename MaskVec::value_type>(it.row());
    };
    static constexpr auto proj_mask = [](const auto &tup) -> decltype(auto) {
        return std::get<1>(tup);
    };
    return util::iter_set_intersection(std::move(col_range),
                                       util::enumerate(std::move(mask_span)),
                                       std::less{}, proj_row, proj_mask);
}
} // namespace detail

template <Config Conf>
void CasADiControlProblem<Conf>::eval_add_R_masked(index_t, crvec xu, crvec h,
                                                   crindexvec mask, rmat R,
                                                   rvec work) const {
    auto &&sparse = impl->R.fun.sparsity_out(0);
    assert(xu.size() == nx + nu);
    assert(h.size() == nh);
    assert(R.rows() == nu);
    assert(R.cols() == nu);
    assert(work.size() >= sparse.nnz());
    impl->R({xu.data(), h.data(), param.data()}, {work.data()});
    using spmat   = Eigen::SparseMatrix<real_t, Eigen::ColMajor, casadi_int>;
    using cmspmat = Eigen::Map<const spmat>;
    if (sparse.is_dense())
        R += cmmat{work.data(), nu, nu}(mask, mask);
    else {
        cmspmat R_full{
            nu, nu, sparse.nnz(), sparse.colind(), sparse.row(), work.data(),
        };
        // Iterate over all columns in the mask
        for (index_t c : mask)
            // Iterate over rows in intersection of mask and sparse column
            for (auto r : detail::select_rows_in_col(R_full, mask, c))
                R(r.row(), c) += r.value();
    }
}

template <Config Conf>
void CasADiControlProblem<Conf>::eval_add_S_masked(index_t, crvec xu, crvec h,
                                                   crindexvec mask, rmat S,
                                                   rvec work) const {
    auto &&sparse = impl->S.fun.sparsity_out(0);
    assert(xu.size() == nx + nu);
    assert(h.size() == nh);
    assert(S.rows() == nu);
    assert(S.cols() == nx);
    assert(work.size() >= sparse.nnz());
    impl->S({xu.data(), h.data(), param.data()}, {work.data()});
    using spmat   = Eigen::SparseMatrix<real_t, Eigen::ColMajor, casadi_int>;
    using cmspmat = Eigen::Map<const spmat>;
    if (sparse.is_dense())
        S += cmmat{work.data(), nu, nx}(mask, Eigen::all);
    else {
        cmspmat S_full{
            nu, nx, sparse.nnz(), sparse.colind(), sparse.row(), work.data(),
        };
        // Iterate over all columns
        for (index_t c = 0; c < S_full.cols(); ++c)
            // Iterate over rows in intersection of mask and sparse column
            for (auto r : detail::select_rows_in_col(S_full, mask, c))
                S(r.row(), c) += r.value();
    }
}

template <Config Conf>
void CasADiControlProblem<Conf>::eval_add_R_prod_masked(index_t, crvec, crvec,
                                                        crindexvec mask_J,
                                                        crindexvec mask_K,
                                                        crvec v, rvec out,
                                                        rvec work) const {
    auto &&sparse = impl->R.fun.sparsity_out(0);
    assert(v.size() == nu);
    assert(out.size() == mask_J.size());
    assert(work.size() >= sparse.nnz());
    using spmat   = Eigen::SparseMatrix<real_t, Eigen::ColMajor, casadi_int>;
    using cmspmat = Eigen::Map<const spmat>;
    if (sparse.is_dense()) {
        auto R_JK = cmmat{work.data(), nu, nu}(mask_J, mask_K);
        out.noalias() += R_JK * v(mask_K);
    } else {
        cmspmat R_full{
            nu, nu, sparse.nnz(), sparse.colind(), sparse.row(), work.data(),
        };
        using detail::select_rows_in_col_iota;
        // out += R_JK * v_K;
        // Iterate over all columns in the mask K
        for (index_t c : mask_K)
            // Iterate over rows in intersection of mask J and sparse column
            for (auto &&[r, rJ] : select_rows_in_col_iota(R_full, mask_J, c))
                out(std::get<0>(rJ)) += r.value() * v(c);
    }
}

template <Config Conf>
void CasADiControlProblem<Conf>::eval_add_S_prod_masked(index_t, crvec, crvec,
                                                        crindexvec mask_K,
                                                        crvec v, rvec out,
                                                        rvec work) const {
    auto &&sparse = impl->R.fun.sparsity_out(0);
    assert(v.size() == nu);
    assert(out.size() == nx);
    assert(work.size() >= sparse.nnz());
    using spmat   = Eigen::SparseMatrix<real_t, Eigen::ColMajor, casadi_int>;
    using cmspmat = Eigen::Map<const spmat>;
    if (sparse.is_dense()) {
        auto Sᵀ_K = cmmat{work.data(), nu, nu}.transpose()(Eigen::all, mask_K);
        out.noalias() += Sᵀ_K * v(mask_K);
    } else {
        cmspmat S_full{
            nu, nu, sparse.nnz(), sparse.colind(), sparse.row(), work.data(),
        };
        // out += Sᵀ_K * v(mask_K);
        // Iterate over all rows of Sᵀ
        for (index_t c = 0; c < S_full.cols(); ++c)
            // Iterate over columns in intersection of mask K and sparse row
            for (auto r : detail::select_rows_in_col(S_full, mask_K, c))
                out(c) += r.value() * v(r.row());
    }
}

template <Config Conf>
auto CasADiControlProblem<Conf>::get_R_work_size() const -> length_t {
    auto &&sparse = impl->R.fun.sparsity_out(0);
    return sparse.nnz();
}

template <Config Conf>
auto CasADiControlProblem<Conf>::get_S_work_size() const -> length_t {
    auto &&sparse = impl->S.fun.sparsity_out(0);
    return sparse.nnz();
}

template <Config Conf>
void CasADiControlProblem<Conf>::eval_constr(index_t, crvec x, rvec c) const {
    if (nc == 0)
        return;
    assert(x.size() == nx);
    assert(c.size() == nc);
    impl->c({x.data(), param.data()}, {c.data()});
}

template <Config Conf>
void CasADiControlProblem<Conf>::eval_jac_constr(index_t, crvec, rmat) const {
    throw std::logic_error("'eval_jac_constr' not implemented");
}

template <Config Conf>
void CasADiControlProblem<Conf>::eval_grad_constr_prod(index_t, crvec x,
                                                       crvec p,
                                                       rvec grad_cx_p) const {
    assert(x.size() == nx);
    assert(p.size() == nc);
    assert(grad_cx_p.size() == nx);
    impl->grad_c_prod({x.data(), param.data(), p.data()}, {grad_cx_p.data()});
}

template <Config Conf>
void CasADiControlProblem<Conf>::eval_add_gn_hess_constr(index_t, crvec x,
                                                         crvec M,
                                                         rmat out) const {
    auto &&sparse = impl->gn_hess_c.fun.sparsity_out(0);
    assert(x.size() == nx);
    assert(M.size() == nc);
    assert(out.rows() == nx);
    assert(out.cols() == nx);
    assert(work.size() >= sparse.nnz());
    impl->gn_hess_c({x.data(), param.data(), M.data()}, {work.data()});
    using spmat   = Eigen::SparseMatrix<real_t, Eigen::ColMajor, casadi_int>;
    using cmspmat = Eigen::Map<const spmat>;
    if (sparse.is_dense())
        out += cmmat{work.data(), nx, nx};
    else
        out += cmspmat{
            nx, nx, sparse.nnz(), sparse.colind(), sparse.row(), work.data(),
        };
}
} // namespace alpaqa::experimental