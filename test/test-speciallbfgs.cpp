#include "eigen-matchers.hpp"

#include <panoc-alm/inner/lbfgs.hpp>

[[maybe_unused]] const static auto printvec = [](const char *name,
                                                 const auto &vec) {
    std::cout << name << ": " << vec.transpose() << std::endl;
};

TEST(SpecializedLBFGS, constantgamma) {
    using pa::inf;
    using pa::real_t;
    using pa::vec;
    pa::SpecializedLBFGS l(2, 3);

    vec x0 = vec::Zero(2);
    vec g0(2);
    g0 << -1, 1;
    real_t γ = 1;
    pa::Box C;
    C.lowerbound = vec::Constant(2, 0.);
    C.upperbound = vec::Constant(2, inf);
    vec x̂0       = pa::project(x0 - γ * g0, C);

    l.initialize(x0, g0, x̂0, γ);

    EXPECT_EQ(l.x(0), x0);
    EXPECT_EQ(l.g(0), g0);
    EXPECT_EQ(l.x̂(), x̂0);

    // ------

    vec x1(2), g1(2);
    x1 << 2, -1;
    g1 << -2, 1;
    vec x̂1 = pa::project(x1 - γ * g1, C);

    l.update(x1, g1, x̂1, C, γ);
    EXPECT_EQ(l.x(0), x0);
    EXPECT_EQ(l.g(0), g0);

    EXPECT_EQ(l.x(1), x1);
    EXPECT_EQ(l.g(1), g1);
    EXPECT_EQ(l.x̂(), x̂1);

    EXPECT_EQ(l.s(0), x1 - x0);
    EXPECT_EQ(l.y(0), x1 - x0 + x̂0 - x̂1);

    // ------

    vec x2(2), g2(2);
    x2 << 5, -3;
    g2 << -3, -1;
    vec x̂2 = pa::project(x2 - γ * g2, C);

    l.update(x2, g2, x̂2, C, γ);
    EXPECT_EQ(l.x(0), x0);
    EXPECT_EQ(l.g(0), g0);

    EXPECT_EQ(l.x(1), x1);
    EXPECT_EQ(l.g(1), g1);

    EXPECT_EQ(l.s(0), x1 - x0);
    EXPECT_EQ(l.y(0), x1 - x0 + x̂0 - x̂1);

    EXPECT_EQ(l.x(2), x2);
    EXPECT_EQ(l.g(2), g2);
    EXPECT_EQ(l.x̂(), x̂2);

    EXPECT_EQ(l.s(1), x2 - x1);
    EXPECT_EQ(l.y(1), x2 - x1 + x̂1 - x̂2);

    // ------

    vec x3(2), g3(2);
    x3 << 4, -3;
    g3 << -9, -1;
    vec x̂3 = pa::project(x3 - γ * g3, C);

    l.update(x3, g3, x̂3, C, γ);
    EXPECT_EQ(l.x(1), x1);
    EXPECT_EQ(l.g(1), g1);

    EXPECT_EQ(l.s(0), x1 - x0);
    EXPECT_EQ(l.y(0), x1 - x0 + x̂0 - x̂1);

    EXPECT_EQ(l.x(2), x2);
    EXPECT_EQ(l.g(2), g2);

    EXPECT_EQ(l.s(1), x2 - x1);
    EXPECT_EQ(l.y(1), x2 - x1 + x̂1 - x̂2);

    EXPECT_EQ(l.x(0), x3);
    EXPECT_EQ(l.g(0), g3);
    EXPECT_EQ(l.x̂(), x̂3);

    EXPECT_EQ(l.s(2), x3 - x2);
    EXPECT_EQ(l.y(2), x3 - x2 + x̂2 - x̂3);

    // ------

    vec x4(2), g4(2);
    x4 << -1, -8;
    g4 << 3, -15;
    vec x̂4 = pa::project(x4 - γ * g4, C);

    l.update(x4, g4, x̂4, C, γ);
    EXPECT_EQ(l.x(2), x2);
    EXPECT_EQ(l.g(2), g2);

    EXPECT_EQ(l.s(1), x2 - x1);
    EXPECT_EQ(l.y(1), x2 - x1 + x̂1 - x̂2);

    EXPECT_EQ(l.x(0), x3);
    EXPECT_EQ(l.g(0), g3);

    EXPECT_EQ(l.s(2), x3 - x2);
    EXPECT_EQ(l.y(2), x3 - x2 + x̂2 - x̂3);

    EXPECT_EQ(l.x(1), x4);
    EXPECT_EQ(l.g(1), g4);
    EXPECT_EQ(l.x̂(), x̂4);

    EXPECT_EQ(l.s(0), x4 - x3);
    EXPECT_EQ(l.y(0), x4 - x3 + x̂3 - x̂4);
}

TEST(SpecializedLBFGS, updategamma1) {
    using pa::inf;
    using pa::real_t;
    using pa::vec;
    pa::SpecializedLBFGS l(2, 3);

    vec x0 = vec::Zero(2);
    vec g0(2);
    g0 << -1, 1;
    real_t γ = 1;
    pa::Box C;
    C.lowerbound = vec::Constant(2, 0.);
    C.upperbound = vec::Constant(2, inf);
    vec x̂0       = pa::project(x0 - γ * g0, C);

    l.initialize(x0, g0, x̂0, γ);

    EXPECT_EQ(l.x(0), x0);
    EXPECT_EQ(l.g(0), g0);
    EXPECT_EQ(l.x̂(), x̂0);

    // ------

    vec x1(2), g1(2);
    x1 << 2, -1;
    g1 << -2, 1;
    vec x̂1 = pa::project(x1 - γ * g1, C);

    l.update(x1, g1, x̂1, C, γ);
    EXPECT_EQ(l.x(0), x0);
    EXPECT_EQ(l.g(0), g0);

    EXPECT_EQ(l.x(1), x1);
    EXPECT_EQ(l.g(1), g1);
    EXPECT_EQ(l.x̂(), x̂1);

    EXPECT_EQ(l.s(0), x1 - x0);
    EXPECT_EQ(l.y(0), x1 - x0 + x̂0 - x̂1);

    // ------ GAMMA CHANGE

    vec x2(2), g2(2);
    x2 << 5, -3;
    g2 << -3, -1;

    γ      = 0.5;
    x̂0     = pa::project(x0 - γ * g0, C);
    x̂1     = pa::project(x1 - γ * g1, C);
    vec x̂2 = pa::project(x2 - γ * g2, C);

    l.update(x2, g2, x̂2, C, γ);
    EXPECT_EQ(l.x(0), x0);
    EXPECT_EQ(l.g(0), g0);

    EXPECT_EQ(l.x(1), x1);
    EXPECT_EQ(l.g(1), g1);

    EXPECT_EQ(l.s(0), x1 - x0);
    EXPECT_EQ(l.y(0), x1 - x0 + x̂0 - x̂1);

    EXPECT_EQ(l.x(2), x2);
    EXPECT_EQ(l.g(2), g2);
    EXPECT_EQ(l.x̂(), x̂2);

    EXPECT_EQ(l.s(1), x2 - x1);
    EXPECT_EQ(l.y(1), x2 - x1 + x̂1 - x̂2);

    // ------

    vec x3(2), g3(2);
    x3 << 4, -3;
    g3 << -9, -1;
    vec x̂3 = pa::project(x3 - γ * g3, C);

    l.update(x3, g3, x̂3, C, γ);
    EXPECT_EQ(l.x(1), x1);
    EXPECT_EQ(l.g(1), g1);

    EXPECT_EQ(l.s(0), x1 - x0);
    EXPECT_EQ(l.y(0), x1 - x0 + x̂0 - x̂1);

    EXPECT_EQ(l.x(2), x2);
    EXPECT_EQ(l.g(2), g2);

    EXPECT_EQ(l.s(1), x2 - x1);
    EXPECT_EQ(l.y(1), x2 - x1 + x̂1 - x̂2);

    EXPECT_EQ(l.x(0), x3);
    EXPECT_EQ(l.g(0), g3);
    EXPECT_EQ(l.x̂(), x̂3);

    EXPECT_EQ(l.s(2), x3 - x2);
    EXPECT_EQ(l.y(2), x3 - x2 + x̂2 - x̂3);

    // ------

    vec x4(2), g4(2);
    x4 << -1, -8;
    g4 << 3, -15;
    vec x̂4 = pa::project(x4 - γ * g4, C);

    l.update(x4, g4, x̂4, C, γ);
    EXPECT_EQ(l.x(2), x2);
    EXPECT_EQ(l.g(2), g2);

    EXPECT_EQ(l.s(1), x2 - x1);
    EXPECT_EQ(l.y(1), x2 - x1 + x̂1 - x̂2);

    EXPECT_EQ(l.x(0), x3);
    EXPECT_EQ(l.g(0), g3);

    EXPECT_EQ(l.s(2), x3 - x2);
    EXPECT_EQ(l.y(2), x3 - x2 + x̂2 - x̂3);

    EXPECT_EQ(l.x(1), x4);
    EXPECT_EQ(l.g(1), g4);
    EXPECT_EQ(l.x̂(), x̂4);

    EXPECT_EQ(l.s(0), x4 - x3);
    EXPECT_EQ(l.y(0), x4 - x3 + x̂3 - x̂4);
}

TEST(SpecializedLBFGS, updategamma2) {
    using pa::inf;
    using pa::real_t;
    using pa::vec;
    pa::SpecializedLBFGS l(2, 3);

    vec x0 = vec::Zero(2);
    vec g0(2);
    g0 << -1, 1;
    real_t γ = 1;
    pa::Box C;
    C.lowerbound = vec::Constant(2, 0.);
    C.upperbound = vec::Constant(2, inf);
    vec x̂0       = pa::project(x0 - γ * g0, C);

    l.initialize(x0, g0, x̂0, γ);

    EXPECT_EQ(l.x(0), x0);
    EXPECT_EQ(l.g(0), g0);
    EXPECT_EQ(l.x̂(), x̂0);

    // ------

    vec x1(2), g1(2);
    x1 << 2, -1;
    g1 << -2, 1;
    vec x̂1 = pa::project(x1 - γ * g1, C);

    l.update(x1, g1, x̂1, C, γ);
    EXPECT_EQ(l.x(0), x0);
    EXPECT_EQ(l.g(0), g0);

    EXPECT_EQ(l.x(1), x1);
    EXPECT_EQ(l.g(1), g1);
    EXPECT_EQ(l.x̂(), x̂1);

    EXPECT_EQ(l.s(0), x1 - x0);
    EXPECT_EQ(l.y(0), x1 - x0 + x̂0 - x̂1);

    // ------

    vec x2(2), g2(2);
    x2 << 5, -3;
    g2 << -3, -1;
    vec x̂2 = pa::project(x2 - γ * g2, C);

    l.update(x2, g2, x̂2, C, γ);
    EXPECT_EQ(l.x(0), x0);
    EXPECT_EQ(l.g(0), g0);

    EXPECT_EQ(l.x(1), x1);
    EXPECT_EQ(l.g(1), g1);

    EXPECT_EQ(l.s(0), x1 - x0);
    EXPECT_EQ(l.y(0), x1 - x0 + x̂0 - x̂1);

    EXPECT_EQ(l.x(2), x2);
    EXPECT_EQ(l.g(2), g2);
    EXPECT_EQ(l.x̂(), x̂2);

    EXPECT_EQ(l.s(1), x2 - x1);
    EXPECT_EQ(l.y(1), x2 - x1 + x̂1 - x̂2);

    // ------

    vec x3(2), g3(2);
    x3 << 4, -3;
    g3 << -9, -1;
    vec x̂3 = pa::project(x3 - γ * g3, C);

    l.update(x3, g3, x̂3, C, γ);
    EXPECT_EQ(l.x(1), x1);
    EXPECT_EQ(l.g(1), g1);

    EXPECT_EQ(l.s(0), x1 - x0);
    EXPECT_EQ(l.y(0), x1 - x0 + x̂0 - x̂1);

    EXPECT_EQ(l.x(2), x2);
    EXPECT_EQ(l.g(2), g2);

    EXPECT_EQ(l.s(1), x2 - x1);
    EXPECT_EQ(l.y(1), x2 - x1 + x̂1 - x̂2);

    EXPECT_EQ(l.x(0), x3);
    EXPECT_EQ(l.g(0), g3);
    EXPECT_EQ(l.x̂(), x̂3);

    EXPECT_EQ(l.s(2), x3 - x2);
    EXPECT_EQ(l.y(2), x3 - x2 + x̂2 - x̂3);

    // ------

    vec x4(2), g4(2);
    x4 << -1, -8;
    g4 << 3, -15;
    vec x̂4 = pa::project(x4 - γ * g4, C);

    l.update(x4, g4, x̂4, C, γ);
    EXPECT_EQ(l.x(2), x2);
    EXPECT_EQ(l.g(2), g2);

    EXPECT_EQ(l.s(1), x2 - x1);
    EXPECT_EQ(l.y(1), x2 - x1 + x̂1 - x̂2);

    EXPECT_EQ(l.x(0), x3);
    EXPECT_EQ(l.g(0), g3);

    EXPECT_EQ(l.s(2), x3 - x2);
    EXPECT_EQ(l.y(2), x3 - x2 + x̂2 - x̂3);

    EXPECT_EQ(l.x(1), x4);
    EXPECT_EQ(l.g(1), g4);
    EXPECT_EQ(l.x̂(), x̂4);

    EXPECT_EQ(l.s(0), x4 - x3);
    EXPECT_EQ(l.y(0), x4 - x3 + x̂3 - x̂4);

    // ------ GAMMA CHANGE

    vec x5(2), g5(2);
    x5 << 11, 10;
    g5 << 6, -2;
    γ      = 2;
    x̂2     = pa::project(x2 - γ * g2, C);
    x̂3     = pa::project(x3 - γ * g3, C);
    x̂4     = pa::project(x4 - γ * g4, C);
    vec x̂5 = pa::project(x5 - γ * g5, C);

    l.update(x5, g5, x̂5, C, γ);
    EXPECT_EQ(l.x(0), x3);
    EXPECT_EQ(l.g(0), g3);

    EXPECT_EQ(l.s(2), x3 - x2);
    EXPECT_EQ(l.y(2), x3 - x2 + x̂2 - x̂3);

    EXPECT_EQ(l.x(1), x4);
    EXPECT_EQ(l.g(1), g4);

    EXPECT_EQ(l.s(0), x4 - x3);
    EXPECT_EQ(l.y(0), x4 - x3 + x̂3 - x̂4);

    EXPECT_EQ(l.x(2), x5);
    EXPECT_EQ(l.g(2), g5);
    EXPECT_EQ(l.x̂(), x̂5);

    EXPECT_EQ(l.s(1), x5 - x4);
    EXPECT_EQ(l.y(1), x5 - x4 + x̂4 - x̂5);
}