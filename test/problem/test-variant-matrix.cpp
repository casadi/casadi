#include <alpaqa/config/config.hpp>
#include <alpaqa/problem/variant-matrix.hpp>
#include <gtest/gtest.h>
#include <test-util/eigen-matchers.hpp>

template <class Mat>
auto mat_to_span(Mat &m) {
    return std::span{m.data(), static_cast<size_t>(m.size())};
}

template <class Mat>
    requires requires(Mat m) { m.valuePtr(); }
auto mat_to_span(Mat &m) {
    return std::span{m.valuePtr(), static_cast<size_t>(m.nonZeros())};
}

TEST(VariantMatrix, Dense) {
    USING_ALPAQA_CONFIG(alpaqa::DefaultConfig);
    using varmat = alpaqa::VariantMatrix<config_t>;
    length_t m = 11, n = 17, p = 13;
    mat A = mat::Random(m, n);
    varmat A_struc{alpaqa::DenseMatrixStructure<config_t>{m, n}};
    {
        mat B = mat::Random(n, p);
        mat C = A * B;
        mat C_result{m, p};
        C_result.setZero();
        alpaqa::multiply_add(mat_to_span(A), A_struc, B, C_result);
        EXPECT_DOUBLE_EQ(C(0, 0), C_result(0, 0));
        EXPECT_THAT(C, EigenAlmostEqual(C_result, 1e-14));
    }
    {
        mat B = mat::Random(m, p);
        mat C = A.transpose() * B;
        mat C_result{n, p};
        C_result.setZero();
        alpaqa::multiply_add(mat_to_span(A), A_struc.transpose(), B, C_result);
        EXPECT_DOUBLE_EQ(C(0, 0), C_result(0, 0));
        EXPECT_THAT(C, EigenAlmostEqual(C_result, 1e-14));
    }
    {
        mat B = mat::Random(p, m);
        mat C = B * A;
        mat C_result{p, n};
        C_result.setZero();
        alpaqa::multiply_add(mat_to_span(A), B, A_struc, C_result);
        EXPECT_DOUBLE_EQ(C(0, 0), C_result(0, 0));
        EXPECT_THAT(C, EigenAlmostEqual(C_result, 1e-14));
    }
    {
        mat B = mat::Random(p, n);
        mat C = B * A.transpose();
        mat C_result{p, m};
        C_result.setZero();
        alpaqa::multiply_add(mat_to_span(A), B, A_struc.transpose(), C_result);
        EXPECT_DOUBLE_EQ(C(0, 0), C_result(0, 0));
        EXPECT_THAT(C, EigenAlmostEqual(C_result, 1e-14));
    }
}

TEST(VariantMatrix, Diagonal) {
    USING_ALPAQA_CONFIG(alpaqa::DefaultConfig);
    using varmat = alpaqa::VariantMatrix<config_t>;
    length_t m = 11, p = 13;
    vec A_sto = vec::Random(m);
    mat A     = A_sto.asDiagonal();
    varmat A_struc{alpaqa::DiagonalMatrixStructure<config_t>{m, m}};
    {
        mat B = mat::Random(m, p);
        mat C = A * B;
        mat C_result{m, p};
        C_result.setZero();
        alpaqa::multiply_add(mat_to_span(A_sto), A_struc, B, C_result);
        EXPECT_DOUBLE_EQ(C(0, 0), C_result(0, 0));
        EXPECT_THAT(C, EigenAlmostEqual(C_result, 1e-14));
    }
    {
        mat B = mat::Random(m, p);
        mat C = A * B;
        mat C_result{m, p};
        C_result.setZero();
        alpaqa::multiply_add(mat_to_span(A_sto), A_struc.transpose(), B,
                             C_result);
        EXPECT_DOUBLE_EQ(C(0, 0), C_result(0, 0));
        EXPECT_THAT(C, EigenAlmostEqual(C_result, 1e-14));
    }
    {
        mat B = mat::Random(p, m);
        mat C = B * A;
        mat C_result{p, m};
        C_result.setZero();
        alpaqa::multiply_add(mat_to_span(A_sto), B, A_struc, C_result);
        EXPECT_DOUBLE_EQ(C(0, 0), C_result(0, 0));
        EXPECT_THAT(C, EigenAlmostEqual(C_result, 1e-14));
    }
    {
        mat B = mat::Random(p, m);
        mat C = B * A;
        mat C_result{p, m};
        C_result.setZero();
        alpaqa::multiply_add(mat_to_span(A_sto), B, A_struc.transpose(),
                             C_result);
        EXPECT_DOUBLE_EQ(C(0, 0), C_result(0, 0));
        EXPECT_THAT(C, EigenAlmostEqual(C_result, 1e-14));
    }
}

TEST(VariantMatrix, Symmetric) {
    USING_ALPAQA_CONFIG(alpaqa::DefaultConfig);
    using varmat = alpaqa::VariantMatrix<config_t>;
    length_t m = 11, p = 13;
    mat A_sto = mat::Random(m, m);
    mat A     = A_sto.selfadjointView<Eigen::Upper>();
    varmat A_struc{
        alpaqa::SymmetricMatrixStructure<config_t, Eigen::Upper>{m, m}};
    {
        mat B = mat::Random(m, p);
        mat C = A * B;
        mat C_result{m, p};
        C_result.setZero();
        alpaqa::multiply_add(mat_to_span(A_sto), A_struc, B, C_result);
        EXPECT_DOUBLE_EQ(C(0, 0), C_result(0, 0));
        EXPECT_THAT(C, EigenAlmostEqual(C_result, 1e-14));
    }
    {
        mat B = mat::Random(m, p);
        mat C = A * B;
        mat C_result{m, p};
        C_result.setZero();
        alpaqa::multiply_add(mat_to_span(A_sto), A_struc.transpose(), B,
                             C_result);
        EXPECT_DOUBLE_EQ(C(0, 0), C_result(0, 0));
        EXPECT_THAT(C, EigenAlmostEqual(C_result, 1e-14));
    }
    {
        mat B = mat::Random(p, m);
        mat C = B * A;
        mat C_result{p, m};
        C_result.setZero();
        alpaqa::multiply_add(mat_to_span(A_sto), B, A_struc, C_result);
        EXPECT_DOUBLE_EQ(C(0, 0), C_result(0, 0));
        EXPECT_THAT(C, EigenAlmostEqual(C_result, 1e-14));
    }
    {
        mat B = mat::Random(p, m);
        mat C = B * A;
        mat C_result{p, m};
        C_result.setZero();
        alpaqa::multiply_add(mat_to_span(A_sto), B, A_struc.transpose(),
                             C_result);
        EXPECT_DOUBLE_EQ(C(0, 0), C_result(0, 0));
        EXPECT_THAT(C, EigenAlmostEqual(C_result, 1e-14));
    }
}

TEST(VariantMatrix, Sparse) {
    USING_ALPAQA_CONFIG(alpaqa::DefaultConfig);
    using varmat = alpaqa::VariantMatrix<config_t>;
    length_t m = 11, n = 17, p = 13;
    using spmat    = Eigen::SparseMatrix<real_t, Eigen::ColMajor, int>;
    spmat A        = mat::Random(m, n).sparseView();
    using sp_struc = alpaqa::SparseMatrixStructure<config_t, int>;
    varmat A_struc{
        sp_struc{std::make_shared<sp_struc::Sparsity>(sp_struc::Sparsity{
            .rows                = m,
            .cols                = n,
            .nnz                 = A.nonZeros(),
            .outer_index_ptr     = A.outerIndexPtr(),
            .inner_index_ptr     = A.innerIndexPtr(),
            .inner_non_zeros_ptr = A.innerNonZeroPtr(),
        })},
    };
    {
        mat B = mat::Random(n, p);
        mat C = A * B;
        mat C_result{m, p};
        C_result.setZero();
        alpaqa::multiply_add(mat_to_span(A), A_struc, B, C_result);
        EXPECT_DOUBLE_EQ(C(0, 0), C_result(0, 0));
        EXPECT_THAT(C, EigenAlmostEqual(C_result, 1e-14));
    }
    {
        mat B = mat::Random(m, p);
        mat C = A.transpose() * B;
        mat C_result{n, p};
        C_result.setZero();
        alpaqa::multiply_add(mat_to_span(A), A_struc.transpose(), B, C_result);
        EXPECT_DOUBLE_EQ(C(0, 0), C_result(0, 0));
        EXPECT_THAT(C, EigenAlmostEqual(C_result, 1e-14));
    }
    {
        mat B = mat::Random(p, m);
        mat C = B * A;
        mat C_result{p, n};
        C_result.setZero();
        alpaqa::multiply_add(mat_to_span(A), B, A_struc, C_result);
        EXPECT_DOUBLE_EQ(C(0, 0), C_result(0, 0));
        EXPECT_THAT(C, EigenAlmostEqual(C_result, 1e-14));
    }
    {
        mat B = mat::Random(p, n);
        mat C = B * A.transpose();
        mat C_result{p, m};
        C_result.setZero();
        alpaqa::multiply_add(mat_to_span(A), B, A_struc.transpose(), C_result);
        EXPECT_DOUBLE_EQ(C(0, 0), C_result(0, 0));
        EXPECT_THAT(C, EigenAlmostEqual(C_result, 1e-14));
    }
}