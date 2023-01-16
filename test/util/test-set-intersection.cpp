#if !defined(__clang_major__) || __clang_major__ > 15 || defined(__clangd__)

#include <gtest/gtest.h>

#include <alpaqa/util/set-intersection.hpp>
#include <algorithm>
#include <iterator>

TEST(SetIntersection, sameFirst) {
    std::array set_A{1, 3, 10, 22, 23, 24, 25, 30, 32};
    std::array set_B{1, 2, 3, 24, 25, 26, 30, 31, 33};
    std::vector<int> result_A;
    std::vector<int> result_B;
    std::vector<int> expected;
    std::ranges::set_intersection(set_A, set_B, std::back_inserter(expected));
    for (auto [a, b] : alpaqa::util::iter_set_intersection(set_A, set_B)) {
        result_A.push_back(a);
        result_B.push_back(b);
    }
    EXPECT_EQ(result_A, expected);
    EXPECT_EQ(result_B, expected);
}

TEST(SetIntersection, sameLast) {
    std::array set_A{3, 10, 22, 23, 24, 25, 30, 32, 34};
    std::array set_B{2, 3, 24, 25, 26, 30, 31, 33, 34};
    std::vector<int> result_A;
    std::vector<int> result_B;
    std::vector<int> expected;
    std::ranges::set_intersection(set_A, set_B, std::back_inserter(expected));
    for (auto [a, b] : alpaqa::util::iter_set_intersection(set_A, set_B)) {
        result_A.push_back(a);
        result_B.push_back(b);
    }
    EXPECT_EQ(result_A, expected);
    EXPECT_EQ(result_B, expected);
}

TEST(SetIntersection, differentFirstLast) {
    std::array set_A{1, 3, 10, 22, 23, 24, 25, 30, 32};
    std::array set_B{0, 2, 3, 24, 25, 26, 30, 31, 33};
    std::vector<int> result_A;
    std::vector<int> result_B;
    std::vector<int> expected;
    std::ranges::set_intersection(set_A, set_B, std::back_inserter(expected));
    for (auto [a, b] : alpaqa::util::iter_set_intersection(set_A, set_B)) {
        result_A.push_back(a);
        result_B.push_back(b);
    }
    EXPECT_EQ(result_A, expected);
    EXPECT_EQ(result_B, expected);
}

#endif
