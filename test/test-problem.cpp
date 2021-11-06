#include <iomanip>
#include <alpaqa-ref/fd.hpp>
#include <alpaqa/util/problem.hpp>

#include "eigen-matchers.hpp"

using namespace alpaqa;

TEST(Problem, onlyD) {
    Problem original;
    original.n = 3;
    original.m = 2;
    original.C.lowerbound.resize(original.n);
    original.C.upperbound.resize(original.n);
    original.D.lowerbound.resize(original.m);
    original.D.upperbound.resize(original.m);
    original.C.lowerbound << -13, -7, -3;
    original.C.upperbound << 11, 17, 23;
    original.D.lowerbound << -51, 5;
    original.D.upperbound << 41, 37;
    original.f      = [](crvec x) { return 43 * x(0) + 41 * x(1) + 61 * x(2); };
    original.grad_f = [](crvec, rvec g) {
        g(0) = 43;
        g(1) = 41;
        g(2) = 61;
    };
    original.g = [](crvec x, rvec g) {
        g(0) = 71 * x(0) * x(1) + 83 * x(1) * x(2);
        g(1) = 97 * x(1) * x(1) + 29 * x(0) * x(2);
    };
    original.grad_g_prod = [](crvec x, crvec y, rvec g) {
        mat grad(3, 2);
        grad << 71 * x(1), 29 * x(2),             //
            71 * x(0) + 83 * x(2), 2 * 97 * x(1), //
            83 * x(1), 29 * x(0);                 //
        g = grad * y;
    };

    auto g_component = [](const Problem &problem) {
        return [&problem](size_t component) {
            return [&problem, component, g = vec(problem.m)](crvec x) mutable {
                problem.g(x, g);
                return g(component);
            };
        };
    };

    vec x(3);
    x << 13, 23, 31;
    vec fd_f = pa_ref::finite_diff(original.f, x);
    vec gr_f(3);
    original.grad_f(x, gr_f);
    vec g(2);
    original.g(x, g);
    ASSERT_THAT(print_wrap(gr_f), EigenAlmostEqual(print_wrap(fd_f), 1e-6));
    auto orig_g_component = g_component(original);
    vec fd_g0             = pa_ref::finite_diff(orig_g_component(0), x);
    vec fd_g1             = pa_ref::finite_diff(orig_g_component(1), x);
    vec gr_g0(3), gr_g1(3);
    original.grad_g_prod(x, mat::Identity(2, 2).col(0), gr_g0);
    original.grad_g_prod(x, mat::Identity(2, 2).col(1), gr_g1);
    ASSERT_THAT(print_wrap(gr_g0), EigenAlmostEqual(print_wrap(fd_g0), 1e-4));
    ASSERT_THAT(print_wrap(gr_g1), EigenAlmostEqual(print_wrap(fd_g1), 1e-4));
    // std::cout << std::setprecision(17);
    // std::cout << "gr f: " << gr_f.transpose() << std::endl;
    // std::cout << "fd f: " << fd_f.transpose() << std::endl;
    // std::cout << "gr g0: " << gr_g0.transpose() << std::endl;
    // std::cout << "fd g0: " << fd_g0.transpose() << std::endl;
    // std::cout << "gr g1: " << gr_g1.transpose() << std::endl;
    // std::cout << "fd g1: " << fd_g1.transpose() << std::endl;

    auto onlyD = ProblemOnlyD(original);

    ASSERT_EQ(onlyD.n, 3);
    ASSERT_EQ(onlyD.m, 2 + 3);

    vec Cu_expected = vec::Constant(3, inf);
    vec Cl_expected = vec::Constant(3, -inf);
    vec Du_expected(2 + 3);
    vec Dl_expected(2 + 3);
    Du_expected << 41, 37, 11, 17, 23;
    Dl_expected << -51, 5, -13, -7, -3;
    EXPECT_THAT(print_wrap(onlyD.D.lowerbound),
                EigenEqual(print_wrap(Dl_expected)));
    EXPECT_THAT(print_wrap(onlyD.D.upperbound),
                EigenEqual(print_wrap(Du_expected)));

    EXPECT_EQ(original.f(x), onlyD.f(x));
    vec gr_fD(3);
    onlyD.grad_f(x, gr_fD);
    EXPECT_THAT(print_wrap(gr_fD), EigenEqual(print_wrap(gr_f)));
    vec gD(2 + 3);
    vec gD_expected(2 + 3);
    gD_expected << g(0), g(1), x(0), x(1), x(2);
    onlyD.g(x, gD);
    EXPECT_THAT(print_wrap(gD), EigenEqual(print_wrap(gD_expected)));

    auto gD_component = g_component(onlyD);
    vec fd_gD0        = pa_ref::finite_diff(gD_component(0), x);
    vec fd_gD1        = pa_ref::finite_diff(gD_component(1), x);
    vec fd_gD2        = pa_ref::finite_diff(gD_component(2), x);
    vec fd_gD3        = pa_ref::finite_diff(gD_component(3), x);
    vec fd_gD4        = pa_ref::finite_diff(gD_component(4), x);
    vec gr_gD0(3), gr_gD1(3), gr_gD2(3), gr_gD3(3), gr_gD4(3);
    onlyD.grad_g_prod(x, mat::Identity(5, 5).col(0), gr_gD0);
    onlyD.grad_g_prod(x, mat::Identity(5, 5).col(1), gr_gD1);
    onlyD.grad_g_prod(x, mat::Identity(5, 5).col(2), gr_gD2);
    onlyD.grad_g_prod(x, mat::Identity(5, 5).col(3), gr_gD3);
    onlyD.grad_g_prod(x, mat::Identity(5, 5).col(4), gr_gD4);
    EXPECT_THAT(print_wrap(gr_gD0), EigenAlmostEqual(print_wrap(fd_gD0), 1e-4));
    EXPECT_THAT(print_wrap(gr_gD1), EigenAlmostEqual(print_wrap(fd_gD1), 1e-4));
    EXPECT_THAT(print_wrap(gr_gD2), EigenAlmostEqual(print_wrap(fd_gD2), 1e-4));
    EXPECT_THAT(print_wrap(gr_gD3), EigenAlmostEqual(print_wrap(fd_gD3), 1e-4));
    EXPECT_THAT(print_wrap(gr_gD4), EigenAlmostEqual(print_wrap(fd_gD4), 1e-4));
}
