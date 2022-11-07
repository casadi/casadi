#include <gtest/gtest.h>

#include <alpaqa/config/config.hpp>
#include <alpaqa/util/index-set.hpp>
#include <algorithm>

TEST(IndexSet, simple) {
    USING_ALPAQA_CONFIG(alpaqa::DefaultConfig);
    constexpr int n = 10;
    constexpr int N = 5;
    Eigen::Matrix<bool, N, n> values_T;
    values_T << 0, 0, 1, 1, 1, 1, 0, 1, 1, 0, //
        1, 1, 0, 1, 0, 0, 0, 1, 0, 0,         //
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0,         //
        1, 1, 1, 1, 1, 1, 1, 1, 1, 1,         //
        1, 1, 0, 0, 0, 1, 0, 1, 0, 1;         //
    auto values = values_T.transpose();
    alpaqa::detail::IndexSet<config_t> set{N, n};
    set.update([&](index_t t, index_t i) { return values(i, t); });
    for (index_t i = 0; i < N; ++i) {
        auto J = set.indices(i);
        auto K = set.compl_indices(i);
        std::cout << '[' << values.col(i).transpose() << "]\n";
        std::cout << '[' << J.transpose() << "]\n";
        std::cout << '[' << K.transpose() << "]\n";
        ASSERT_EQ(J.size() + K.size(), n);
        ASSERT_LT(J.end()[-1], n);
        ASSERT_TRUE(std::is_sorted(J.begin(), J.end()));
        EXPECT_EQ(std::adjacent_find(J.begin(), J.end()), J.end());
        for (index_t j : J) {
            EXPECT_TRUE(values(j, i));
        }
        ASSERT_LT(K.end()[-1], n);
        ASSERT_TRUE(std::is_sorted(K.begin(), K.end()));
        EXPECT_EQ(std::adjacent_find(K.begin(), K.end()), K.end());
        for (index_t k : K) {
            EXPECT_FALSE(values(k, i));
        }
        std::cout << "---\n";
    }
}