#include <gtest/gtest.h>

#include <alpaqa/inner/decl/lbfgs-stepsize.hpp>
#include <alpaqa/inner/directions/lbfgs.hpp>
#include <alpaqa/standalone/panoc.hpp>
#include <alpaqa/util/alloc.hpp>

#include <chrono>
#include <cmath>

TEST(PANOCStandalone, cosh) {
    using namespace alpaqa;

    auto f = [](crvec x_) {
        auto x = x_(0);
        return std::cosh(x) - x * x + x;
    };
    auto grad_f = [](crvec x_, rvec g) {
        auto x = x_(0);
        g(0)   = std::sinh(x) - 2 * x + 1;
    };
    Box C;
    C.lowerbound = vec::Constant(1, -2.5);
    C.upperbound = vec::Constant(1, 10);

    alpaqa::PANOCParams params;
    params.max_iter       = 32;
    params.τ_min          = 1. / 16;
    params.print_interval = 0;
    alpaqa::LBFGSParams lbfgsparams;
    lbfgsparams.memory = 20;
    real_t ε           = 1e-14;

    vec_allocator alloc{50, 1};
    vec x = vec(1);
    x(0)  = 5;
    std::chrono::microseconds total_duration{0};
    for (size_t i = 0; i < 10'000; ++i) {
        auto stats =
            alpaqa::panoc(f, grad_f, C, x, ε, params, {lbfgsparams}, alloc);
        total_duration += stats.elapsed_time;
        x(0) = x(0) + 5 - x(0);
    }

    auto stats = alpaqa::panoc(f, grad_f, C, x, ε, params, {lbfgsparams}, alloc);

    EXPECT_EQ(alloc.highwater(), panoc_min_alloc_size);
    EXPECT_EQ(alloc.size(), 50);
    EXPECT_EQ(alloc.used_space(), 0);

    EXPECT_DOUBLE_EQ(x(0), -2.4877411259719642);

    std::cout << std::setprecision(16) << std::fixed;
    std::cout << "\n===========\n" << std::endl;
    std::cout << "f(x)     = " << f(x) << std::endl;
    std::cout << "x        = " << x.transpose() << std::endl;
    std::cout << "Iter:      " << stats.iterations << std::endl;
    std::cout << "Highwater: " << alloc.highwater() << std::endl;
    std::cout << "Time:      " << total_duration.count() << std::endl;
    std::cout << "Average τ: " << stats.sum_τ / stats.count_τ << std::endl;
    std::cout << "Status:    " << stats.status << std::endl << std::endl;
}