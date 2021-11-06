#include <alpaqa/util/box.hpp>

#include "eigen-matchers.hpp"

TEST(Box, project) {
    pa::Box b{pa::vec(3), pa::vec(3)};
    b.upperbound << 10, 11, 12;
    b.lowerbound << -12, -11, -10;
    pa::vec v(3);
    v << 11, -12, -5;

    pa::vec result(3);
    pa::vec expected(3);
    expected << 10, -11, -5;

    result = pa::project(v, b);
    EXPECT_THAT(print_wrap(result), EigenEqual(print_wrap(expected)));
}
