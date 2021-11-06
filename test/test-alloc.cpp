#include <new>
#include <alpaqa/util/alloc.hpp>

#include "eigen-matchers.hpp"

TEST(Alloc, allocrepeat) {
    alpaqa::vec_allocator alloc{3, 2};
    auto v0 = alloc.alloc();
    EXPECT_EQ(v0.size(), 2);
    v0(0) = 1.2;
    v0(1) = 1.3;
    EXPECT_EQ(v0(0), 1.2);
    EXPECT_EQ(v0(1), 1.3);
    auto p0 = v0.data();
    alloc.free(v0);
    auto v1 = alloc.alloc();
    EXPECT_EQ(v1.size(), 2);
    v1(0) = 1.2;
    v1(1) = 1.3;
    EXPECT_EQ(v1(0), 1.2);
    EXPECT_EQ(v1(1), 1.3);
    auto p1 = v1.data();
    EXPECT_EQ(p1, p0);
    alloc.free(v1);
    EXPECT_EQ(alloc.highwater(), 1);
}

TEST(Alloc, badalloc) {
    alpaqa::vec_allocator alloc{3, 2};
    auto v0 = alloc.alloc();
    auto v1 = alloc.alloc();
    auto v2 = alloc.alloc();
    try {
        auto v3 = alloc.alloc();
        alloc.free(v3);
        FAIL();
    } catch (std::bad_alloc &) {
        // expected
    }
    alloc.free(v2);
    alloc.free(v1);
    alloc.free(v0);
    EXPECT_EQ(alloc.highwater(), 3);
}

TEST(Alloc, raii) {
    alpaqa::vec_allocator alloc{3, 2};
    {
        auto v0 = alloc.alloc_raii();
        EXPECT_EQ(alloc.size(), 2);
    }
    EXPECT_EQ(alloc.size(), 3);
}