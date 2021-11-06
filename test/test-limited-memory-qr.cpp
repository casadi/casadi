#include "eigen-matchers.hpp"
#include "alpaqa/util/vec.hpp"
#include <gtest/gtest.h>

#include <iomanip>
#include <iostream>
#include <alpaqa/inner/detail/limited-memory-qr.hpp>

TEST(LimitedMemoryQR, adding) {
    using alpaqa::mat;
    using alpaqa::real_t;
    using alpaqa::vec;

    mat A(4, 3);
    A << 3, 3, 2, //
        4, 4, 1,  //
        0, 6, 2,  //
        0, 8, 1;

    constexpr real_t ε = 1e-15;
    std::cout << std::setprecision(16) << std::scientific;
    alpaqa::LimitedMemoryQR qr(4, 3);

    qr.add_column(A.col(0));
    EXPECT_THAT(print_wrap(A.block(0, 0, 4, 1)),
                EigenAlmostEqual(print_wrap(qr.get_Q() * qr.get_R()), ε));

    std::cout << qr.get_Q() * qr.get_R() << "\r\n\n";

    qr.add_column(A.col(1));
    EXPECT_THAT(print_wrap(A.block(0, 0, 4, 2)),
                EigenAlmostEqual(print_wrap(qr.get_Q() * qr.get_R()), ε));

    std::cout << std::setprecision(16);
    std::cout << qr.get_Q() * qr.get_R() << "\r\n\n";

    qr.add_column(A.col(2));
    EXPECT_THAT(print_wrap(A.block(0, 0, 4, 3)),
                EigenAlmostEqual(print_wrap(qr.get_Q() * qr.get_R()), ε));

    std::cout << std::setprecision(16);
    std::cout << qr.get_Q() * qr.get_R() << "\r\n\n";

    const real_t root2 = std::sqrt(2);
    mat Q_exp(4, 3);
    Q_exp << 3. / 5, 0, 0.4 * root2, //
        4. / 5, 0, -0.3 * root2,     //
        0, 3. / 5, 0.4 * root2,      //
        0, 4. / 5, -0.3 * root2;
    mat R_exp(3, 3);
    R_exp << 5, 5, 2, //
        0, 10, 2,     //
        0, 0, root2;

    EXPECT_THAT(print_wrap(Q_exp), EigenAlmostEqual(print_wrap(qr.get_Q()), ε));
    EXPECT_THAT(print_wrap(R_exp), EigenAlmostEqual(print_wrap(qr.get_R()), ε));
}

TEST(LimitedMemoryQR, removing) {
    using alpaqa::mat;
    using alpaqa::real_t;
    using alpaqa::vec;

    mat A(4, 3);
    A << 3, 3, 2, //
        4, 4, 1,  //
        0, 6, 2,  //
        0, 8, 1;

    constexpr real_t ε = 1e-15;
    std::cout << std::setprecision(16) << std::scientific;
    alpaqa::LimitedMemoryQR qr(4, 3);

    qr.add_column(A.col(0));
    qr.add_column(A.col(1));
    qr.add_column(A.col(2));
    EXPECT_THAT(print_wrap(A.block(0, 0, 4, 3)),
                EigenAlmostEqual(print_wrap(qr.get_Q() * qr.get_R()), ε));

    qr.remove_column();
    EXPECT_THAT(print_wrap(A.block(0, 1, 4, 2)),
                EigenAlmostEqual(print_wrap(qr.get_Q() * qr.get_R()), ε));
    std::cout << qr.get_Q() * qr.get_R() << "\r\n\n";

    qr.remove_column();
    EXPECT_THAT(print_wrap(A.block(0, 2, 4, 1)),
                EigenAlmostEqual(print_wrap(qr.get_Q() * qr.get_R()), ε));
    std::cout << qr.get_Q() * qr.get_R() << "\r\n\n";
}

TEST(LimitedMemoryQR, addingremoving) {
    using alpaqa::mat;
    using alpaqa::real_t;
    using alpaqa::vec;

    mat A(4, 11);
    A << 3, 3, 2, 2, 1, 1, 0, 3, 7, 1, 2,  //
        4, 4, 1, 1, 3, 5, 1, 7, 4, 2, 2,   //
        0, 6, 2, 4, 1, 6, 1, -3, 6, 5, -5, //
        0, 8, 1, 7, 6, -4, 3, 1, 1, 5, 5;

    constexpr real_t ε = 1e-14;
    std::cout << std::setprecision(16) << std::scientific;
    alpaqa::LimitedMemoryQR qr(4, 3);

    qr.add_column(A.col(0));
    qr.add_column(A.col(1));
    qr.add_column(A.col(2));

    qr.remove_column();
    qr.add_column(A.col(3));
    EXPECT_THAT(print_wrap(A.block(0, 1, 4, 3)),
                EigenAlmostEqual(print_wrap(qr.get_Q() * qr.get_R()), ε));
    std::cout << qr.get_Q() * qr.get_R() << "\r\n\n";

    qr.remove_column();
    qr.add_column(A.col(4));
    EXPECT_THAT(print_wrap(A.block(0, 2, 4, 3)),
                EigenAlmostEqual(print_wrap(qr.get_Q() * qr.get_R()), ε));
    std::cout << qr.get_Q() * qr.get_R() << "\r\n\n";

    qr.remove_column();
    qr.add_column(A.col(5));
    EXPECT_THAT(print_wrap(A.block(0, 3, 4, 3)),
                EigenAlmostEqual(print_wrap(qr.get_Q() * qr.get_R()), ε));
    std::cout << qr.get_Q() * qr.get_R() << "\r\n\n";

    qr.remove_column();
    qr.add_column(A.col(6));
    EXPECT_THAT(print_wrap(A.block(0, 4, 4, 3)),
                EigenAlmostEqual(print_wrap(qr.get_Q() * qr.get_R()), ε));
    std::cout << qr.get_Q() * qr.get_R() << "\r\n\n";

    qr.remove_column();
    qr.add_column(A.col(7));
    EXPECT_THAT(print_wrap(A.block(0, 5, 4, 3)),
                EigenAlmostEqual(print_wrap(qr.get_Q() * qr.get_R()), ε));
    std::cout << qr.get_Q() * qr.get_R() << "\r\n\n";

    qr.remove_column();
    qr.add_column(A.col(8));
    EXPECT_THAT(print_wrap(A.block(0, 6, 4, 3)),
                EigenAlmostEqual(print_wrap(qr.get_Q() * qr.get_R()), ε));
    std::cout << qr.get_Q() * qr.get_R() << "\r\n\n";

    qr.remove_column();
    qr.add_column(A.col(9));
    EXPECT_THAT(print_wrap(A.block(0, 7, 4, 3)),
                EigenAlmostEqual(print_wrap(qr.get_Q() * qr.get_R()), ε));
    std::cout << qr.get_Q() * qr.get_R() << "\r\n\n";

    qr.remove_column();
    qr.add_column(A.col(10));
    EXPECT_THAT(print_wrap(A.block(0, 8, 4, 3)),
                EigenAlmostEqual(print_wrap(qr.get_Q() * qr.get_R()), ε));
    std::cout << qr.get_Q() * qr.get_R() << "\r\n\n";
}

TEST(LimitedMemoryQR, solve) {
    using alpaqa::mat;
    using alpaqa::real_t;
    using alpaqa::vec;

    mat A(4, 3);
    A << 3, 3, 2, //
        4, 4, 1,  //
        0, 6, 2,  //
        0, 8, 1;

    constexpr real_t ε = 1e-15;
    std::cout << std::setprecision(16) << std::scientific;
    alpaqa::LimitedMemoryQR qr(4, 3);

    qr.add_column(vec::Constant(4, 13));
    qr.add_column(A.col(0));
    qr.add_column(A.col(1));
    qr.remove_column();
    qr.add_column(A.col(2));
    EXPECT_THAT(print_wrap(qr.get_Q() * qr.get_R()),
                EigenAlmostEqual(print_wrap(A.block(0, 0, 4, 3)), ε));

    mat b(4, 2);
    b << 1, 2, 3, 4, //
        4, 3, 2, 1;
    mat x1(3, 2);
    qr.solve(b, x1);
    auto x2 = qr.solve(b);
    std::cout << "x1:\n" << x1.transpose() << std::endl;
    std::cout << "x2:\n" << x1.transpose() << std::endl;

    mat x_exp = A.householderQr().solve(b);
    std::cout << "x_exp:\n" << x_exp.transpose() << std::endl;

    EXPECT_THAT(print_wrap(x1), EigenAlmostEqual(print_wrap(x_exp), ε));
    EXPECT_THAT(print_wrap(x2), EigenAlmostEqual(print_wrap(x_exp), ε));
}

TEST(LimitedMemoryQR, scale_R) {
    using alpaqa::mat;
    using alpaqa::real_t;
    using alpaqa::vec;

    mat A(4, 3);
    A << 3, 3, 2, //
        4, 4, 1,  //
        0, 6, 2,  //
        0, 8, 1;

    constexpr real_t ε = 1e-15;
    std::cout << std::setprecision(16) << std::scientific;
    alpaqa::LimitedMemoryQR qr(4, 3);

    qr.add_column(vec::Constant(4, 13));
    qr.add_column(A.col(0));
    qr.add_column(A.col(1));
    qr.remove_column();
    qr.add_column(A.col(2));
    qr.scale_R(0.5);
    EXPECT_THAT(print_wrap(qr.get_Q() * qr.get_R()),
                EigenAlmostEqual(print_wrap(0.5 * A.block(0, 0, 4, 3)), ε));
}