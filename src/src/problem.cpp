#include <panoc-alm/problem.hpp>

namespace pa {

void ProblemWithCounters::attach_counters(ProblemWithCounters &wc) {
    wc.f = [&wc, f{std::move(wc.f)}](const vec &x) {
        ++wc.evaluations.f;
        return f(x);
    };
    wc.grad_f = [&wc, grad_f{std::move(wc.grad_f)}](const vec &x, vec &grad) {
        ++wc.evaluations.grad_f;
        grad_f(x, grad);
    };
    wc.g = [&wc, g{std::move(wc.g)}](const vec &x, vec &gx) {
        ++wc.evaluations.g;
        g(x, gx);
    };
    wc.grad_g = [&wc, grad_g{std::move(wc.grad_g)}](const vec &x, const vec &y,
                                                    vec &grad) {
        ++wc.evaluations.grad_g;
        grad_g(x, y, grad);
    };
}

} // namespace pa