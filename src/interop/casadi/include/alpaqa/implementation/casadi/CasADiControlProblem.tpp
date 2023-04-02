#pragma once

#include <alpaqa/casadi/CasADiControlProblem.hpp>
#include <alpaqa/casadi/CasADiFunctionWrapper.hpp>
#include <alpaqa/config/config.hpp>
#include <alpaqa/util/io/csv.hpp>
#include <alpaqa/util/not-implemented.hpp>
#include <alpaqa/util/sparse-ops.hpp>
#include "CasADiLoader-util.hpp"

#include <Eigen/Sparse>
#include <casadi/core/external.hpp>

#include <memory>
#include <optional>
#include <stdexcept>
#include <type_traits>

namespace alpaqa {

namespace fs = std::filesystem;

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

    CasADiFunctionEvaluator<Conf, 1 + WithParam, 1> c_N;
    CasADiFunctionEvaluator<Conf, 2 + WithParam, 1> grad_c_prod_N;
    CasADiFunctionEvaluator<Conf, 2 + WithParam, 1> gn_hess_c_N;
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
        nc = static_cast<length_t>(cfun.size1_out(0));
        CasADiFunctionEvaluator<Conf, 2, 1> c{std::move(cfun)};
        c.validate_dimensions({dim(nx, 1), dim(p, 1)}, {dim(nc, 1)});
        return c;
    };
    auto load_c_N = [&]() -> CasADiFunctionEvaluator<Conf, 2, 1> {
        casadi::Function cfun = casadi::external("c_N", so_name);
        using namespace std::literals::string_literals;
        if (cfun.n_in() != 2)
            throw std::invalid_argument(
                "Invalid number of input arguments: got "s +
                std::to_string(cfun.n_in()) + ", should be 2.");
        if (cfun.n_out() != 1)
            throw std::invalid_argument(
                "Invalid number of output arguments: got "s +
                std::to_string(cfun.n_in()) + ", should be 1.");
        nc_N = static_cast<length_t>(cfun.size1_out(0));
        CasADiFunctionEvaluator<Conf, 2, 1> c{std::move(cfun)};
        c.validate_dimensions({dim(nx, 1), dim(p, 1)}, {dim(nc_N, 1)});
        return c;
    };
    // Load the functions "f", "h", and "c" to determine the unknown dimensions.
    auto f   = wrap_load(so_name, "f", load_f);
    auto h   = wrap_load(so_name, "h", load_h);
    auto h_N = wrap_load(so_name, "h_N", load_h_N);
    auto c   = wrap_load(so_name, "c", load_c);
    auto c_N = wrap_load(so_name, "c_N", load_c_N);

    this->x_init = vec::Constant(nx, alpaqa::NaN<Conf>);
    this->param  = vec::Constant(p, alpaqa::NaN<Conf>);
    this->U      = Box{nu};
    this->D      = Box{nc};
    this->D_N    = Box{nc_N};

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
            .c_N           = std::move(c_N),
            .grad_c_prod_N = wrapped_load<CasADiFunctionEvaluator<Conf, 3, 1>>(
                so_name, "grad_c_prod_N", dims(nx, p, nc_N), dims(nx)),
            .gn_hess_c_N = wrapped_load<CasADiFunctionEvaluator<Conf, 3, 1>>(
                so_name, "gn_hess_c_N", dims(nx, p, nc_N), dims(dim{nx, nx})),
        });

    auto n_work = std::max({
        impl->Q.fun.sparsity_out(0).nnz(),
        impl->Q_N.fun.sparsity_out(0).nnz(),
        impl->gn_hess_c.fun.sparsity_out(0).nnz(),
        impl->gn_hess_c_N.fun.sparsity_out(0).nnz(),
    });
    this->work  = vec::Constant(static_cast<length_t>(n_work), NaN<Conf>);

    auto bounds_filepath = fs::path{so_name}.replace_extension("csv");
    if (fs::exists(bounds_filepath))
        load_numerical_data(bounds_filepath);
}

template <Config Conf>
void CasADiControlProblem<Conf>::load_numerical_data(
    const std::filesystem::path &filepath, char sep) {
    std::ifstream data_file{filepath};
    if (!data_file)
        throw std::runtime_error("Unable to open bounds file \"" +
                                 filepath.string() + '"');
    index_t line          = 0;
    auto wrap_bounds_load = [&](std::string_view name, auto &v) {
        try {
            ++line;
            csv::read_row(data_file, v, sep);
        } catch (csv::read_error &e) {
            throw std::runtime_error("Unable to read " + std::string(name) +
                                     " from bounds file \"" +
                                     filepath.string() + ':' +
                                     std::to_string(line) + "\": " + e.what());
        }
    };
    wrap_bounds_load("U.lowerbound", this->U.lowerbound);
    wrap_bounds_load("U.upperbound", this->U.upperbound);
    wrap_bounds_load("D.lowerbound", this->D.lowerbound);
    wrap_bounds_load("D.upperbound", this->D.upperbound);
    wrap_bounds_load("D_N.lowerbound", this->D_N.lowerbound);
    wrap_bounds_load("D_N.upperbound", this->D_N.upperbound);
    wrap_bounds_load("x_init", this->x_init);
    wrap_bounds_load("param", this->param);
}

template <Config Conf>
CasADiControlProblem<Conf>::CasADiControlProblem(const CasADiControlProblem &) =
    default;
template <Config Conf>
CasADiControlProblem<Conf> &
CasADiControlProblem<Conf>::operator=(const CasADiControlProblem &) = default;

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
    assert(h.size() == nh_N);
    impl->h_N({x.data(), param.data()}, {h.data()});
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
    assert(h.size() == nh_N);
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
    assert(h.size() == nh_N);
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
            nx,
            nx,
            static_cast<length_t>(sparse.nnz()),
            sparse.colind(),
            sparse.row(),
            work.data(),
        };
}
template <Config Conf>
void CasADiControlProblem<Conf>::eval_add_Q_N(crvec x, crvec h, rmat Q) const {
    assert(x.size() == nx);
    assert(h.size() == nh_N);
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
            nx,
            nx,
            static_cast<length_t>(sparse.nnz()),
            sparse.colind(),
            sparse.row(),
            work.data(),
        };
}

template <Config Conf>
void CasADiControlProblem<Conf>::eval_add_R_masked(index_t, crvec xu, crvec h,
                                                   crindexvec mask, rmat R,
                                                   rvec work) const {
    auto &&sparse = impl->R.fun.sparsity_out(0);
    assert(xu.size() == nx + nu);
    assert(h.size() == nh);
    assert(R.rows() <= nu);
    assert(R.cols() <= nu);
    assert(R.rows() == mask.size());
    assert(R.cols() == mask.size());
    assert(work.size() >= static_cast<length_t>(sparse.nnz()));
    impl->R({xu.data(), h.data(), param.data()}, {work.data()});
    using spmat   = Eigen::SparseMatrix<real_t, Eigen::ColMajor, casadi_int>;
    using cmspmat = Eigen::Map<const spmat>;
    if (sparse.is_dense()) {
        cmmat R_full{work.data(), nu, nu};
        R += R_full(mask, mask);
    } else {
        cmspmat R_full{
            nu,
            nu,
            static_cast<length_t>(sparse.nnz()),
            sparse.colind(),
            sparse.row(),
            work.data(),
        };
        util::sparse_add_masked(R_full, R, mask);
    }
}

template <Config Conf>
void CasADiControlProblem<Conf>::eval_add_S_masked(index_t, crvec xu, crvec h,
                                                   crindexvec mask, rmat S,
                                                   rvec work) const {
    auto &&sparse = impl->S.fun.sparsity_out(0);
    assert(xu.size() == nx + nu);
    assert(h.size() == nh);
    assert(S.rows() <= nu);
    assert(S.rows() == mask.size());
    assert(S.cols() == nx);
    assert(work.size() >= static_cast<length_t>(sparse.nnz()));
    impl->S({xu.data(), h.data(), param.data()}, {work.data()});
    using spmat   = Eigen::SparseMatrix<real_t, Eigen::ColMajor, casadi_int>;
    using cmspmat = Eigen::Map<const spmat>;
    using Eigen::indexing::all;
    if (sparse.is_dense()) {
        cmmat S_full{work.data(), nu, nx};
        S += S_full(mask, all);
    } else {
        cmspmat S_full{
            nu,
            nx,
            static_cast<length_t>(sparse.nnz()),
            sparse.colind(),
            sparse.row(),
            work.data(),
        };
        util::sparse_add_masked_rows(S_full, S, mask);
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
    assert(work.size() >= static_cast<length_t>(sparse.nnz()));
    using spmat   = Eigen::SparseMatrix<real_t, Eigen::ColMajor, casadi_int>;
    using cmspmat = Eigen::Map<const spmat>;
    if (sparse.is_dense()) {
        auto R = cmmat{work.data(), nu, nu};
        out.noalias() += R(mask_J, mask_K) * v(mask_K);
    } else {
        cmspmat R{
            nu,
            nu,
            static_cast<length_t>(sparse.nnz()),
            sparse.colind(),
            sparse.row(),
            work.data(),
        };
        // out += R_full(mask_J,mask_K) * v(mask_K);
        util::sparse_matvec_add_masked_rows_cols(R, v, out, mask_J, mask_K);
    }
}

template <Config Conf>
void CasADiControlProblem<Conf>::eval_add_S_prod_masked(index_t, crvec, crvec,
                                                        crindexvec mask_K,
                                                        crvec v, rvec out,
                                                        rvec work) const {
    auto &&sparse = impl->S.fun.sparsity_out(0);
    assert(v.size() == nu);
    assert(out.size() == nx);
    assert(work.size() >= static_cast<length_t>(sparse.nnz()));
    using spmat   = Eigen::SparseMatrix<real_t, Eigen::ColMajor, casadi_int>;
    using cmspmat = Eigen::Map<const spmat>;
    using Eigen::indexing::all;
    if (sparse.is_dense()) {
        auto Sᵀ = cmmat{work.data(), nu, nx}.transpose();
        out.noalias() += Sᵀ(all, mask_K) * v(mask_K);
    } else {
        cmspmat S{
            nu,
            nx,
            static_cast<length_t>(sparse.nnz()),
            sparse.colind(),
            sparse.row(),
            work.data(),
        };
        // out += S(mask_K,:)ᵀ * v(mask_K);
        util::sparse_matvec_add_transpose_masked_rows(S, v, out, mask_K);
    }
}

template <Config Conf>
auto CasADiControlProblem<Conf>::get_R_work_size() const -> length_t {
    auto &&sparse = impl->R.fun.sparsity_out(0);
    return static_cast<length_t>(sparse.nnz());
}

template <Config Conf>
auto CasADiControlProblem<Conf>::get_S_work_size() const -> length_t {
    auto &&sparse = impl->S.fun.sparsity_out(0);
    return static_cast<length_t>(sparse.nnz());
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
    assert(work.size() >= static_cast<length_t>(sparse.nnz()));
    impl->gn_hess_c({x.data(), param.data(), M.data()}, {work.data()});
    using spmat   = Eigen::SparseMatrix<real_t, Eigen::ColMajor, casadi_int>;
    using cmspmat = Eigen::Map<const spmat>;
    if (sparse.is_dense())
        out += cmmat{work.data(), nx, nx};
    else
        out += cmspmat{
            nx,
            nx,
            static_cast<length_t>(sparse.nnz()),
            sparse.colind(),
            sparse.row(),
            work.data(),
        };
}

template <Config Conf>
void CasADiControlProblem<Conf>::eval_constr_N(crvec x, rvec c) const {
    if (nc_N == 0)
        return;
    assert(x.size() == nx);
    assert(c.size() == nc_N);
    impl->c_N({x.data(), param.data()}, {c.data()});
}

template <Config Conf>
void CasADiControlProblem<Conf>::eval_grad_constr_prod_N(crvec x, crvec p,
                                                         rvec grad_cx_p) const {
    assert(x.size() == nx);
    assert(p.size() == nc_N);
    assert(grad_cx_p.size() == nx);
    impl->grad_c_prod_N({x.data(), param.data(), p.data()}, {grad_cx_p.data()});
}

template <Config Conf>
void CasADiControlProblem<Conf>::eval_add_gn_hess_constr_N(crvec x, crvec M,
                                                           rmat out) const {
    auto &&sparse = impl->gn_hess_c.fun.sparsity_out(0);
    assert(x.size() == nx);
    assert(M.size() == nc_N);
    assert(out.rows() == nx);
    assert(out.cols() == nx);
    assert(work.size() >= static_cast<length_t>(sparse.nnz()));
    impl->gn_hess_c_N({x.data(), param.data(), M.data()}, {work.data()});
    using spmat   = Eigen::SparseMatrix<real_t, Eigen::ColMajor, casadi_int>;
    using cmspmat = Eigen::Map<const spmat>;
    if (sparse.is_dense())
        out += cmmat{work.data(), nx, nx};
    else
        out += cmspmat{
            nx,
            nx,
            static_cast<length_t>(sparse.nnz()),
            sparse.colind(),
            sparse.row(),
            work.data(),
        };
}

} // namespace alpaqa