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
    wc.grad_g_prod = [&wc, grad_g_prod{std::move(wc.grad_g_prod)}](
                         const vec &x, const vec &y, vec &grad) {
        ++wc.evaluations.grad_g_prod;
        grad_g_prod(x, y, grad);
    };
    wc.grad_gi = [&wc, grad_gi{std::move(wc.grad_gi)}](const vec &x, unsigned i,
                                                       vec &grad) {
        ++wc.evaluations.grad_g_prod;
        grad_gi(x, i, grad);
    };
    wc.hess_L_prod = [&wc, hess_L_prod{std::move(wc.hess_L_prod)}](
                         const vec &x, const vec &y, const vec &v, vec &Hv) {
        ++wc.evaluations.hess_L_prod;
        hess_L_prod(x, y, v, Hv);
    };
    wc.hess_L = [&wc, hess_L{std::move(wc.hess_L)}](const vec &x, const vec &y,
                                                    mat &H) {
        ++wc.evaluations.hess_L;
        hess_L(x, y, H);
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
    this->grad_g_prod = [this](const vec &x, const vec &y, vec &gg) {
        work = y.topRows(original.m);
        original.grad_g_prod(x, work, gg);
        gg += y.bottomRows(original.n);
    };
    this->grad_gi = [](const vec &, unsigned, vec &) {
        throw std::logic_error("ProblemOnlyD::grad_gi: Not implemented");
    };
    this->hess_L_prod = [](const vec &, const vec &, const vec &, vec &) {
        throw std::logic_error("ProblemOnlyD::hess_L_prod: Not implemented");
    };
    this->hess_L = [](const vec &, const vec &, mat &) {
        throw std::logic_error("ProblemOnlyD::hess_L: Not implemented");
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