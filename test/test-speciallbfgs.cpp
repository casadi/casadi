#include "eigen-matchers.hpp"

#include <alpaqa/inner/detail/panoc-helpers.hpp>
#include <alpaqa/inner/directions/specialized-lbfgs.hpp>

[[maybe_unused]] const static auto printvec = [](const char *name,
                                                 const auto &vec) {
    std::cout << name << ": " << vec.transpose() << std::endl;
};

TEST(SpecializedLBFGS, constantgamma) {
    using alpaqa::inf;
    using alpaqa::real_t;
    using alpaqa::vec;
    alpaqa::SpecializedLBFGS l(alpaqa::LBFGSParams{}, 2, 3);

    vec x0 = vec::Zero(2);
    vec g0(2);
    g0 << 1, 1;
    real_t γ = 1;
    alpaqa::Box C;
    C.lowerbound = vec::Constant(2, 0.);
    C.upperbound = vec::Constant(2, inf);
    vec p0       = alpaqa::detail::projected_gradient_step(C, γ, x0, g0);

    l.initialize(x0, g0);

    EXPECT_EQ(l.x(0), x0);
    EXPECT_EQ(l.g(0), g0);

    // ------

    vec x1(2), g1(2);
    x1 << 2, -1;
    g1 << 2, 3;
    vec p1 = alpaqa::detail::projected_gradient_step(C, γ, x1, g1);

    ASSERT_TRUE(l.update(x0, x1, p0, p1, g1, C, γ));
    EXPECT_EQ(l.x(0), x0);
    EXPECT_EQ(l.g(0), g0);

    EXPECT_EQ(l.x(1), x1);
    EXPECT_EQ(l.g(1), g1);

    EXPECT_EQ(l.s(0), x1 - x0);
    EXPECT_EQ(l.y(0), p0 - p1);

    // ------

    vec x2(2), g2(2);
    x2 << 5, -3;
    g2 << 4, -3;
    vec p2 = alpaqa::detail::projected_gradient_step(C, γ, x2, g2);

    ASSERT_TRUE(l.update(x1, x2, p1, p2, g2, C, γ));
    EXPECT_EQ(l.x(0), x0);
    EXPECT_EQ(l.g(0), g0);

    EXPECT_EQ(l.x(1), x1);
    EXPECT_EQ(l.g(1), g1);

    EXPECT_EQ(l.s(0), x1 - x0);
    EXPECT_EQ(l.y(0), p0 - p1);

    EXPECT_EQ(l.x(2), x2);
    EXPECT_EQ(l.g(2), g2);

    EXPECT_EQ(l.s(1), x2 - x1);
    EXPECT_EQ(l.y(1), p1 - p2);

    // ------

    vec x3(2), g3(2);
    x3 << 4, -3;
    g3 << -9, -1;
    vec p3 = alpaqa::detail::projected_gradient_step(C, γ, x3, g3);

    ASSERT_TRUE(l.update(x2, x3, p2, p3, g3, C, γ));
    EXPECT_EQ(l.x(1), x1);
    EXPECT_EQ(l.g(1), g1);

    EXPECT_EQ(l.s(0), x1 - x0);
    EXPECT_EQ(l.y(0), p0 - p1);

    EXPECT_EQ(l.x(2), x2);
    EXPECT_EQ(l.g(2), g2);

    EXPECT_EQ(l.s(1), x2 - x1);
    EXPECT_EQ(l.y(1), p1 - p2);

    EXPECT_EQ(l.x(0), x3);
    EXPECT_EQ(l.g(0), g3);

    EXPECT_EQ(l.s(2), x3 - x2);
    EXPECT_EQ(l.y(2), p2 - p3);

    // ------

    vec x4(2), g4(2);
    x4 << -1, -8;
    g4 << 3, -15;
    vec p4 = alpaqa::detail::projected_gradient_step(C, γ, x4, g4);

    ASSERT_TRUE(l.update(x3, x4, p3, p4, g4, C, γ));
    EXPECT_EQ(l.x(2), x2);
    EXPECT_EQ(l.g(2), g2);

    EXPECT_EQ(l.s(1), x2 - x1);
    EXPECT_EQ(l.y(1), p1 - p2);

    EXPECT_EQ(l.x(0), x3);
    EXPECT_EQ(l.g(0), g3);

    EXPECT_EQ(l.s(2), x3 - x2);
    EXPECT_EQ(l.y(2), p2 - p3);

    EXPECT_EQ(l.x(1), x4);
    EXPECT_EQ(l.g(1), g4);

    EXPECT_EQ(l.s(0), x4 - x3);
    EXPECT_EQ(l.y(0), p3 - p4);
}

TEST(SpecializedLBFGS, updategamma1) {
    using alpaqa::inf;
    using alpaqa::real_t;
    using alpaqa::vec;
    alpaqa::SpecializedLBFGS l(alpaqa::LBFGSParams{}, 2, 3);

    vec x0 = vec::Zero(2);
    vec g0(2);
    g0 << 1, 1;
    real_t γ = 1;
    alpaqa::Box C;
    C.lowerbound = vec::Constant(2, 0.);
    C.upperbound = vec::Constant(2, inf);
    vec p0       = alpaqa::detail::projected_gradient_step(C, γ, x0, g0);

    l.initialize(x0, g0);

    EXPECT_EQ(l.x(0), x0);
    EXPECT_EQ(l.g(0), g0);

    // ------

    vec x1(2), g1(2);
    x1 << 2, -1;
    g1 << 2, 3;
    vec p1 = alpaqa::detail::projected_gradient_step(C, γ, x1, g1);

    ASSERT_TRUE(l.update(x0, x1, p0, p1, g1, C, γ));
    EXPECT_EQ(l.x(0), x0);
    EXPECT_EQ(l.g(0), g0);

    EXPECT_EQ(l.x(1), x1);
    EXPECT_EQ(l.g(1), g1);

    EXPECT_EQ(l.s(0), x1 - x0);
    EXPECT_EQ(l.y(0), p0 - p1);

    // ------ GAMMA CHANGE

    vec x2(2), g2(2);
    x2 << 5, -3;
    g2 << 4, -1;

    γ      = 0.5;
    vec p2 = alpaqa::detail::projected_gradient_step(C, γ, x2, g2);

    ASSERT_TRUE(l.update(x1, x2, p1, p2, g2, C, γ));
    p0 = alpaqa::detail::projected_gradient_step(C, γ, x0, g0);
    p1 = alpaqa::detail::projected_gradient_step(C, γ, x1, g1);
    EXPECT_EQ(l.x(0), x0);
    EXPECT_EQ(l.g(0), g0);

    EXPECT_EQ(l.x(1), x1);
    EXPECT_EQ(l.g(1), g1);

    EXPECT_EQ(l.s(0), x1 - x0);
    EXPECT_EQ(l.y(0), p0 - p1);

    EXPECT_EQ(l.x(2), x2);
    EXPECT_EQ(l.g(2), g2);

    EXPECT_EQ(l.s(1), x2 - x1);
    EXPECT_EQ(l.y(1), p1 - p2);

    // ------

    vec x3(2), g3(2);
    x3 << 4, -3;
    g3 << -9, -1;
    vec p3 = alpaqa::detail::projected_gradient_step(C, γ, x3, g3);

    ASSERT_TRUE(l.update(x2, x3, p2, p3, g3, C, γ));
    EXPECT_EQ(l.x(1), x1);
    EXPECT_EQ(l.g(1), g1);

    EXPECT_EQ(l.s(0), x1 - x0);
    EXPECT_EQ(l.y(0), p0 - p1);

    EXPECT_EQ(l.x(2), x2);
    EXPECT_EQ(l.g(2), g2);

    EXPECT_EQ(l.s(1), x2 - x1);
    EXPECT_EQ(l.y(1), p1 - p2);

    EXPECT_EQ(l.x(0), x3);
    EXPECT_EQ(l.g(0), g3);

    EXPECT_EQ(l.s(2), x3 - x2);
    EXPECT_EQ(l.y(2), p2 - p3);

    // ------

    vec x4(2), g4(2);
    x4 << -1, -8;
    g4 << 3, -15;
    vec p4 = alpaqa::detail::projected_gradient_step(C, γ, x4, g4);

    ASSERT_TRUE(l.update(x3, x4, p3, p4, g4, C, γ));
    EXPECT_EQ(l.x(2), x2);
    EXPECT_EQ(l.g(2), g2);

    EXPECT_EQ(l.s(1), x2 - x1);
    EXPECT_EQ(l.y(1), p1 - p2);

    EXPECT_EQ(l.x(0), x3);
    EXPECT_EQ(l.g(0), g3);

    EXPECT_EQ(l.s(2), x3 - x2);
    EXPECT_EQ(l.y(2), p2 - p3);

    EXPECT_EQ(l.x(1), x4);
    EXPECT_EQ(l.g(1), g4);

    EXPECT_EQ(l.s(0), x4 - x3);
    EXPECT_EQ(l.y(0), p3 - p4);
}

TEST(SpecializedLBFGS, updategamma2) {
    using alpaqa::inf;
    using alpaqa::real_t;
    using alpaqa::vec;
    alpaqa::SpecializedLBFGS l(alpaqa::LBFGSParams{}, 2, 3);

    vec x0 = vec::Zero(2);
    vec g0(2);
    g0 << 1, 1;
    real_t γ = 1;
    alpaqa::Box C;
    C.lowerbound = vec::Constant(2, 0.);
    C.upperbound = vec::Constant(2, inf);
    vec p0       = alpaqa::detail::projected_gradient_step(C, γ, x0, g0);

    l.initialize(x0, g0);

    EXPECT_EQ(l.x(0), x0);
    EXPECT_EQ(l.g(0), g0);

    // ------

    vec x1(2), g1(2);
    x1 << 2, -1;
    g1 << 2, 3;
    vec p1 = alpaqa::detail::projected_gradient_step(C, γ, x1, g1);

    ASSERT_TRUE(l.update(x0, x1, p0, p1, g1, C, γ));
    EXPECT_EQ(l.x(0), x0);
    EXPECT_EQ(l.g(0), g0);

    EXPECT_EQ(l.x(1), x1);
    EXPECT_EQ(l.g(1), g1);

    EXPECT_EQ(l.s(0), x1 - x0);
    EXPECT_EQ(l.y(0), p0 - p1);

    // ------

    vec x2(2), g2(2);
    x2 << 5, -3;
    g2 << 4, -3;
    vec p2 = alpaqa::detail::projected_gradient_step(C, γ, x2, g2);

    ASSERT_TRUE(l.update(x1, x2, p1, p2, g2, C, γ));
    EXPECT_EQ(l.x(0), x0);
    EXPECT_EQ(l.g(0), g0);

    EXPECT_EQ(l.x(1), x1);
    EXPECT_EQ(l.g(1), g1);

    EXPECT_EQ(l.s(0), x1 - x0);
    EXPECT_EQ(l.y(0), p0 - p1);

    EXPECT_EQ(l.x(2), x2);
    EXPECT_EQ(l.g(2), g2);

    EXPECT_EQ(l.s(1), x2 - x1);
    EXPECT_EQ(l.y(1), p1 - p2);

    // ------ GAMMA CHANGE

    vec x3(2), g3(2);
    x3 << 4, -3;
    g3 << -9, -1;

    γ      = 0.5;
    vec p3 = alpaqa::detail::projected_gradient_step(C, γ, x3, g3);

    ASSERT_TRUE(l.update(x2, x3, p2, p3, g3, C, γ));
    p0 = alpaqa::detail::projected_gradient_step(C, γ, x0, g0);
    p1 = alpaqa::detail::projected_gradient_step(C, γ, x1, g1);
    p2 = alpaqa::detail::projected_gradient_step(C, γ, x2, g2);
    EXPECT_EQ(l.x(1), x1);
    EXPECT_EQ(l.g(1), g1);

    EXPECT_EQ(l.s(0), x1 - x0);
    EXPECT_EQ(l.y(0), p0 - p1);

    EXPECT_EQ(l.x(2), x2);
    EXPECT_EQ(l.g(2), g2);

    EXPECT_EQ(l.s(1), x2 - x1);
    EXPECT_EQ(l.y(1), p1 - p2);

    EXPECT_EQ(l.x(0), x3);
    EXPECT_EQ(l.g(0), g3);

    EXPECT_EQ(l.s(2), x3 - x2);
    EXPECT_EQ(l.y(2), p2 - p3);

    // ------

    vec x4(2), g4(2);
    x4 << -1, -8;
    g4 << 3, -15;
    vec p4 = alpaqa::detail::projected_gradient_step(C, γ, x4, g4);

    ASSERT_TRUE(l.update(x3, x4, p3, p4, g4, C, γ));
    EXPECT_EQ(l.x(2), x2);
    EXPECT_EQ(l.g(2), g2);

    EXPECT_EQ(l.s(1), x2 - x1);
    EXPECT_EQ(l.y(1), p1 - p2);

    EXPECT_EQ(l.x(0), x3);
    EXPECT_EQ(l.g(0), g3);

    EXPECT_EQ(l.s(2), x3 - x2);
    EXPECT_EQ(l.y(2), p2 - p3);

    EXPECT_EQ(l.x(1), x4);
    EXPECT_EQ(l.g(1), g4);

    EXPECT_EQ(l.s(0), x4 - x3);
    EXPECT_EQ(l.y(0), p3 - p4);
}

TEST(SpecializedLBFGS, updategamma3) {
    using alpaqa::inf;
    using alpaqa::real_t;
    using alpaqa::vec;
    alpaqa::SpecializedLBFGS l(alpaqa::LBFGSParams{}, 2, 3);

    vec x0 = vec::Zero(2);
    vec g0(2);
    g0 << 1, 1;
    real_t γ = 1;
    alpaqa::Box C;
    C.lowerbound = vec::Constant(2, 0.);
    C.upperbound = vec::Constant(2, inf);
    vec p0       = alpaqa::detail::projected_gradient_step(C, γ, x0, g0);

    l.initialize(x0, g0);

    EXPECT_EQ(l.x(0), x0);
    EXPECT_EQ(l.g(0), g0);

    // ------

    vec x1(2), g1(2);
    x1 << 2, -1;
    g1 << 2, 3;
    vec p1 = alpaqa::detail::projected_gradient_step(C, γ, x1, g1);

    ASSERT_TRUE(l.update(x0, x1, p0, p1, g1, C, γ));
    EXPECT_EQ(l.x(0), x0);
    EXPECT_EQ(l.g(0), g0);

    EXPECT_EQ(l.x(1), x1);
    EXPECT_EQ(l.g(1), g1);

    EXPECT_EQ(l.s(0), x1 - x0);
    EXPECT_EQ(l.y(0), p0 - p1);

    // ------

    vec x2(2), g2(2);
    x2 << 5, -3;
    g2 << 4, -3;
    vec p2 = alpaqa::detail::projected_gradient_step(C, γ, x2, g2);

    ASSERT_TRUE(l.update(x1, x2, p1, p2, g2, C, γ));
    EXPECT_EQ(l.x(0), x0);
    EXPECT_EQ(l.g(0), g0);

    EXPECT_EQ(l.x(1), x1);
    EXPECT_EQ(l.g(1), g1);

    EXPECT_EQ(l.s(0), x1 - x0);
    EXPECT_EQ(l.y(0), p0 - p1);

    EXPECT_EQ(l.x(2), x2);
    EXPECT_EQ(l.g(2), g2);

    EXPECT_EQ(l.s(1), x2 - x1);
    EXPECT_EQ(l.y(1), p1 - p2);

    // ------

    vec x3(2), g3(2);
    x3 << 4, -3;
    g3 << -9, -1;
    vec p3 = alpaqa::detail::projected_gradient_step(C, γ, x3, g3);

    ASSERT_TRUE(l.update(x2, x3, p2, p3, g3, C, γ));
    EXPECT_EQ(l.x(1), x1);
    EXPECT_EQ(l.g(1), g1);

    EXPECT_EQ(l.s(0), x1 - x0);
    EXPECT_EQ(l.y(0), p0 - p1);

    EXPECT_EQ(l.x(2), x2);
    EXPECT_EQ(l.g(2), g2);

    EXPECT_EQ(l.s(1), x2 - x1);
    EXPECT_EQ(l.y(1), p1 - p2);

    EXPECT_EQ(l.x(0), x3);
    EXPECT_EQ(l.g(0), g3);

    EXPECT_EQ(l.s(2), x3 - x2);
    EXPECT_EQ(l.y(2), p2 - p3);

    // ------ GAMMA CHANGE

    vec x4(2), g4(2);
    x4 << -1, -8;
    g4 << 3, -15;

    γ      = 0.5;
    vec p4 = alpaqa::detail::projected_gradient_step(C, γ, x4, g4);

    ASSERT_TRUE(l.update(x3, x4, p3, p4, g4, C, γ));
    p0 = alpaqa::detail::projected_gradient_step(C, γ, x0, g0);
    p1 = alpaqa::detail::projected_gradient_step(C, γ, x1, g1);
    p2 = alpaqa::detail::projected_gradient_step(C, γ, x2, g2);
    p3 = alpaqa::detail::projected_gradient_step(C, γ, x3, g3);
    EXPECT_EQ(l.x(2), x2);
    EXPECT_EQ(l.g(2), g2);

    EXPECT_EQ(l.s(1), x2 - x1);
    EXPECT_EQ(l.y(1), p1 - p2);

    EXPECT_EQ(l.x(0), x3);
    EXPECT_EQ(l.g(0), g3);

    EXPECT_EQ(l.s(2), x3 - x2);
    EXPECT_EQ(l.y(2), p2 - p3);

    EXPECT_EQ(l.x(1), x4);
    EXPECT_EQ(l.g(1), g4);

    EXPECT_EQ(l.s(0), x4 - x3);
    EXPECT_EQ(l.y(0), p3 - p4);
}