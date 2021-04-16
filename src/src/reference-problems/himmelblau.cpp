#include <memory>
#include <panoc-alm/reference-problems/himmelblau.hpp>

namespace pa {
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
        [sq](const vec &x) {
            return sq(sq(x(0)) + x(1) - 11) + sq(x(0) + sq(x(1)) - 7);
        },
        [sq](const vec &x, vec &g) {
            g(0) =
                2 * (2 * x(0) * (sq(x(0)) + x(1) - 11) + x(0) + sq(x(1)) - 7);
            g(1) =
                2 * (sq(x(0)) + 2 * x(1) * (x(0) + sq(x(1)) - 7) + x(1) - 11);
        },
        [](const vec &, vec &) {},
        [](const vec &, const vec &, vec &grad) { grad.setZero(); },
        {},
        {},
    };
}

} // namespace problems
} // namespace pa