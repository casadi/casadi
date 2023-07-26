#include <alpaqa/config/config.hpp>
#include <Eigen/Sparse>
USING_ALPAQA_CONFIG(alpaqa::DefaultConfig);
using spmat = Eigen::SparseMatrix<real_t, Eigen::ColMajor, Eigen::Index>;

#if ALPAQA_WITH_OCP

#include <test-util/eigen-matchers.hpp>

#include <alpaqa/util/sparse-ops.hpp>

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
    mat expected = O + S_dense(J, Eigen::indexing::all);

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
    vec expected = o + S_dense(J, Eigen::indexing::all).transpose() * v(J);

    vec result = o;
    alpaqa::util::sparse_matvec_add_transpose_masked_rows(
        S, crvec{v}, rvec{result}, crindexvec{J});

    EXPECT_THAT(result, EigenAlmostEqual(expected, 1e-13));
}

#endif

#if ALPAQA_WITH_CXX23_TESTS

#include <alpaqa/util/sparse-ops.hpp>

using namespace alpaqa::util;

#if __cpp_lib_ranges_zip >= 202110L && __cpp_lib_ranges_enumerate >= 202302L

TEST(SparseOps, convertTripletToCCS) {
    // Generate some row and column indices of a sparse matrix
    index_t nnz = 12;
    Eigen::MatrixX<int> rows_cols(nnz, 2);
    rows_cols << 0, 0, //
        2, 0,          //
        4, 0,          //
        1, 1,          //
        3, 1,          //
        5, 1,          //
        2, 2,          //
        3, 3,          //
        5, 3,          //
        0, 5,          //
        2, 5,          //
        5, 5;
    // Convert to compressed column storage format
    indexvec inner(nnz), outer(6 + 1);
    convert_triplets_to_ccs<config_t>(rows_cols.col(0), rows_cols.col(1), inner,
                                      outer);
    std::cout << "Inner: " << inner.transpose() << std::endl;
    std::cout << "Outer: " << outer.transpose() << std::endl;

    // Compare to Eigen's implementation
    using triplet_t = Eigen::Triplet<double>;
    std::vector<triplet_t> triplets(static_cast<size_t>(nnz));
    for (auto &&[i, t] : std::views::enumerate(triplets))
        t = {rows_cols(i, 0), rows_cols(i, 1), static_cast<double>(i + 1)};
    spmat A(6, 6);
    A.setFromTriplets(std::begin(triplets), std::end(triplets));
    std::cout << "\n" << A << std::endl;

    EXPECT_TRUE(std::equal(outer.begin(), outer.end(), A.outerIndexPtr()));
    EXPECT_TRUE(std::equal(inner.begin(), inner.end(), A.innerIndexPtr()));
}

#endif

#endif