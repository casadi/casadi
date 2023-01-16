#if !defined(__clang_major__) || __clang_major__ > 15 || defined(__clangd__)

#include <test-util/eigen-matchers.hpp>

#include <alpaqa/config/config.hpp>
#include <alpaqa/util/sparse-ops.hpp>

USING_ALPAQA_CONFIG(alpaqa::DefaultConfig);
using spmat = Eigen::SparseMatrix<real_t, Eigen::ColMajor, Eigen::Index>;

spmat random_sparse(length_t n, length_t m) {
    mat R = mat::Random(n, m);
    R     = R.unaryExpr([](real_t r) {
        return r > 0.5 ? r - 0.5 : r < -0.5 ? r + 0.5 : 0;
    });
    return R.sparseView();
}

indexvec random_mask(length_t n) {
    vec v = vec::Random(n);
    indexvec J(n);
    index_t i = 0;
    for (auto vv : v) {
        if (vv > 0) {
            J(i) = i;
            ++i;
        }
    }
    return J.topRows(i);
}

TEST(SparseOps, R) {
    length_t n = 1024;
    auto R     = random_sparse(n, n);
    auto J     = random_mask(n);
    mat O      = mat::Random(J.size(), J.size());

    auto R_dense = R.toDense();
    mat expected = O + R_dense(J, J);

    mat result = O;
    alpaqa::util::sparse_add_masked(R, rmat{result}, crindexvec{J});

    EXPECT_THAT(result, EigenEqual(expected));
}

TEST(SparseOps, Rvec) {
    length_t n = 1024;
    vec v      = vec::Random(n);
    auto R     = random_sparse(n, n);
    auto J     = random_mask(n);
    auto K     = random_mask(n);
    vec o      = vec::Random(J.size());

    auto R_dense = R.toDense();
    vec expected = o + R_dense(J, K) * v(K);

    vec result = o;
    alpaqa::util::sparse_matvec_add_masked_rows_cols(
        R, crvec{v}, rvec{result}, crindexvec{J}, crindexvec{K});

    EXPECT_THAT(result, EigenAlmostEqual(expected, 1e-13));
}

TEST(SparseOps, S) {
    length_t n = 1024;
    length_t m = 512;
    auto S     = random_sparse(n, m);
    auto J     = random_mask(n);
    mat O      = mat::Random(J.size(), m);

    auto S_dense = S.toDense();
    mat expected = O + S_dense(J, Eigen::placeholders::all);

    mat result = O;
    alpaqa::util::sparse_add_masked_rows(S, rmat{result}, crindexvec{J});

    EXPECT_THAT(result, EigenEqual(expected));
}

TEST(SparseOps, Svec) {
    length_t n = 1024;
    length_t m = 512;
    vec v      = vec::Random(n);
    auto S     = random_sparse(n, m);
    auto J     = random_mask(n);
    vec o      = vec::Random(m);

    auto S_dense = S.toDense();
    vec expected = o + S_dense(J, Eigen::placeholders::all).transpose() * v(J);

    vec result = o;
    alpaqa::util::sparse_matvec_add_transpose_masked_rows(
        S, crvec{v}, rvec{result}, crindexvec{J});

    EXPECT_THAT(result, EigenAlmostEqual(expected, 1e-13));
}

#endif
