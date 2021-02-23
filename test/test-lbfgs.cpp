#include <panoc-alm/box.hpp>
#include <panoc-alm/lbfgs.hpp>
#include <panoc-alm/vec.hpp>

#include <limits>

#include "eigen-matchers.hpp"
#include <Eigen/LU>

TEST(LBFGS, quadratic) {
    pa::mat H(2, 2);
    H << 2, -1, -1, 3;

    std::cout << "Inverse Hessian: \n" << H.inverse() << std::endl;

    auto f      = [&H](const pa::vec &v) { return 0.5 * v.dot(H * v); };
    auto grad_f = [&H](const pa::vec &v) { return pa::vec(H * v); };

    pa::mat B = pa::mat::Identity(2, 2);

    pa::LBFGS lbfgs(2, 5);
    pa::vec x(2);
    x << 10, -5;
    auto r = grad_f(x);
    for (size_t i = 0; i < 10; ++i) {
        { // Print L-BFGS inverse Hessian estimate
            std::cout << std::endl << i << std::endl;
            std::cout << "x:    " << x.transpose() << std::endl;
            std::cout << "f(x): " << f(x) << std::endl;

            pa::mat H⁻¹ = pa::mat::Identity(2, 2);
            lbfgs.apply(H⁻¹.col(0));
            lbfgs.apply(H⁻¹.col(1));
            std::cout << std::endl << "LB⁻¹ = \n" << H⁻¹ << std::endl;
            std::cout << "B⁻¹  = \n" << B.inverse() << std::endl;
            std::cout << "B    = \n" << B << std::endl;
        }

        pa::vec d = r;
        lbfgs.apply(d);
        x -= d;
        auto r_new = grad_f(x);
        auto y     = r_new - r;
        lbfgs.update(-d, y);

        auto s = -d;
        B      = B + y * y.transpose() / y.dot(s) -
            (B * s) * (s.transpose() * B.transpose()) / (s.transpose() * B * s);

        r = std::move(r_new);
    }
    std::cout << std::endl << "final" << std::endl;
    std::cout << "x:    " << x.transpose() << std::endl;
    std::cout << "f(x): " << f(x) << std::endl;

    EXPECT_NEAR(x(0), 0, 1e-10);
    EXPECT_NEAR(x(1), 0, 1e-10);
}
