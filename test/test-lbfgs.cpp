#include <alpaqa/inner/directions/lbfgs.hpp>

#include <limits>

#include "eigen-matchers.hpp"
#include <Eigen/LU>

TEST(LBFGS, quadratic) {
    alpaqa::mat H(2, 2);
    H << 2, -1, -1, 3;

    std::cout << "Inverse Hessian: \n" << H.inverse() << std::endl;

    auto f      = [&H](alpaqa::crvec v) { return 0.5 * v.dot(H * v); };
    auto grad_f = [&H](alpaqa::crvec v) { return alpaqa::vec(H * v); };

    alpaqa::mat B = alpaqa::mat::Identity(2, 2);

    alpaqa::LBFGSParams param;
    param.memory = 5;
    alpaqa::LBFGS lbfgs(param, 2);
    alpaqa::vec x(2);
    x << 10, -5;
    auto r                = grad_f(x);
    unsigned update_count = 0;
    for (size_t i = 0; i < 10; ++i) {
        { // Print L-BFGS inverse Hessian estimate
            std::cout << std::endl << i << std::endl;
            std::cout << "x:    " << x.transpose() << std::endl;
            std::cout << "f(x): " << f(x) << std::endl;

            alpaqa::mat H⁻¹ = alpaqa::mat::Identity(2, 2);
            if (i > 0) {
                lbfgs.apply(H⁻¹.col(0), 1);
                lbfgs.apply(H⁻¹.col(1), 1);
            }
            std::cout << std::endl << "LB⁻¹ = \n" << H⁻¹ << std::endl;
            std::cout << "B⁻¹  = \n" << B.inverse() << std::endl;
            std::cout << "   " << update_count << std::endl;
        }

        alpaqa::vec d = r;
        if (i > 0)
            lbfgs.apply(d, 1);
        alpaqa::vec x_new = x - d;
        alpaqa::vec r_new = grad_f(x_new);
        lbfgs.update(x, x_new, r, r_new, alpaqa::LBFGS::Sign::Positive);
        ++update_count;

        alpaqa::vec y = r_new - r;
        alpaqa::vec s = -d;
        B         = B + y * y.transpose() / y.dot(s) -
            (B * s) * (s.transpose() * B.transpose()) / (s.transpose() * B * s);

        r = std::move(r_new);
        x = std::move(x_new);
    }
    std::cout << std::endl << "final" << std::endl;
    std::cout << "x:    " << x.transpose() << std::endl;
    std::cout << "f(x): " << f(x) << std::endl;

    EXPECT_NEAR(x(0), 0, 1e-10);
    EXPECT_NEAR(x(1), 0, 1e-10);
}
