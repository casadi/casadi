#include <panoc-alm/util/problem.hpp>

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

void ProblemOnlyD::transform() {
    work.resize(original.m);
    this->n      = original.n;
    this->m      = original.m + original.n;
    this->f      = [this](const vec &x) { return original.f(x); };
    this->grad_f = [this](const vec &x, vec &grad) {
        original.grad_f(x, grad);
    };
    this->g = [this](const vec &x, vec &gg) {
        original.g(x, work);
        gg.topRows(original.m)    = work;
        gg.bottomRows(original.n) = x;
    };
    this->grad_g = [this](const vec &x, const vec &y, vec &gg) {
        work = y.topRows(original.m);
        original.grad_g(x, work, gg);
        gg += y.bottomRows(original.n);
    };
    this->C.lowerbound = vec::Constant(original.n, -inf);
    this->C.upperbound = vec::Constant(original.n, +inf);
    this->D.lowerbound.resize(this->m);
    this->D.upperbound.resize(this->m);
    this->D.lowerbound.topRows(original.m)    = original.D.lowerbound;
    this->D.lowerbound.bottomRows(original.n) = original.C.lowerbound;
    this->D.upperbound.topRows(original.m)    = original.D.upperbound;
    this->D.upperbound.bottomRows(original.n) = original.C.upperbound;
}

} // namespace pa