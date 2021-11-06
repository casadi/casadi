#include <alpaqa/util/box.hpp>

#include "eigen-matchers.hpp"

TEST(Box, project) {
    alpaqa::Box b{alpaqa::vec(3), alpaqa::vec(3)};
    b.upperbound << 10, 11, 12;
    b.lowerbound << -12, -11, -10;
    alpaqa::vec v(3);
    v << 11, -12, -5;

    alpaqa::vec result(3);
    alpaqa::vec expected(3);
    expected << 10, -11, -5;

    result = alpaqa::project(v, b);
    EXPECT_THAT(print_wrap(result), EigenEqual(print_wrap(expected)));
}
