#include <cmath>
#include <limits>

#include <alpaqa-ref/fd.hpp>

namespace pa_ref {

vec finite_diff(std::function<real_t(crvec )> f, crvec x) {
    const auto n = x.size();
    vec grad(n);
    vec h        = vec::Zero(n);
    const auto ε = std::sqrt(std::numeric_limits<real_t>::epsilon());
    const auto δ = std::numeric_limits<real_t>::min() / ε;
    for (unsigned i = 0; i < n; ++i) {
        real_t hh        = std::abs(x(i)) > δ ? x(i) * ε : δ;
        h(i)             = hh;
        grad.coeffRef(i) = (f(x + h) - f(x)) / hh;
        h(i)             = 0;
    }
    return grad;
}

} // namespace pa_ref
