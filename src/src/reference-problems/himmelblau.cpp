#include <memory>
#include <alpaqa/reference-problems/himmelblau.hpp>

namespace alpaqa {
namespace problems {

Problem himmelblau_problem() {
    auto sq = [](auto x) { return x * x; };
    Box C{vec(2), vec(2)};
    C.lowerbound << -1, -1;
    C.upperbound << 4, 1.8;
    return Problem{
        2,
        0,
        C,
        Box{},
        [sq](crvec x) {
            return sq(sq(x(0)) + x(1) - 11) + sq(x(0) + sq(x(1)) - 7);
        },
        [sq](crvec x, rvec g) {
            g(0) =
                2 * (2 * x(0) * (sq(x(0)) + x(1) - 11) + x(0) + sq(x(1)) - 7);
            g(1) =
                2 * (sq(x(0)) + 2 * x(1) * (x(0) + sq(x(1)) - 7) + x(1) - 11);
        },
        [](crvec , rvec ) {},
        [](crvec , crvec , rvec grad) { grad.setZero(); },
        [](crvec , unsigned, rvec grad_gi) { grad_gi.setZero(); },
        [sq](crvec x, crvec , crvec v, rvec Hv) {
            real_t H00 = 4 * (sq(x(0)) + x(1) - 11) + 8 * sq(x(0)) + 2;
            real_t H01 = 4 * x(0) + 4 * x(1);
            real_t H10 = 4 * x(0) + 4 * x(1);
            real_t H11 = 4 * (x(0) + sq(x(1)) - 7) + 8 * sq(x(1)) + 2;
            Hv(0) = H00 * v(0) + H01 * v(1);
            Hv(1) = H10 * v(0) + H11 * v(1);
        },
        [sq](crvec x, crvec , rmat H) {
            H(0, 0) = 4 * (sq(x(0)) + x(1) - 11) + 8 * sq(x(0)) + 2;
            H(0, 1) = 4 * x(0) + 4 * x(1);
            H(1, 0) = 4 * x(0) + 4 * x(1);
            H(1, 1) = 4 * (x(0) + sq(x(1)) - 7) + 8 * sq(x(1)) + 2;
        },
    };
}

} // namespace problems
} // namespace alpaqa