#include <panoc-alm/util/problem.hpp>

namespace pa {

void ProblemWithCounters::attach_counters(ProblemWithCounters &wc) {
    wc.f = [&wc, f{std::move(wc.f)}](crvec x) {
        ++wc.evaluations.f;
        return f(x);
    };
    wc.grad_f = [&wc, grad_f{std::move(wc.grad_f)}](crvec x, rvec grad) {
        ++wc.evaluations.grad_f;
        grad_f(x, grad);
    };
    wc.g = [&wc, g{std::move(wc.g)}](crvec x, rvec gx) {
        ++wc.evaluations.g;
        g(x, gx);
    };
    wc.grad_g_prod = [&wc, grad_g_prod{std::move(wc.grad_g_prod)}](
                         crvec x, crvec y, rvec grad) {
        ++wc.evaluations.grad_g_prod;
        grad_g_prod(x, y, grad);
    };
    wc.grad_gi = [&wc, grad_gi{std::move(wc.grad_gi)}](crvec x, unsigned i,
                                                       rvec grad) {
        ++wc.evaluations.grad_g_prod;
        grad_gi(x, i, grad);
    };
    wc.hess_L_prod = [&wc, hess_L_prod{std::move(wc.hess_L_prod)}](
                         crvec x, crvec y, crvec v, rvec Hv) {
        ++wc.evaluations.hess_L_prod;
        hess_L_prod(x, y, v, Hv);
    };
    wc.hess_L = [&wc, hess_L{std::move(wc.hess_L)}](crvec x, crvec y, rmat H) {
        ++wc.evaluations.hess_L;
        hess_L(x, y, H);
    };
}

void ProblemOnlyD::transform() {
    work.resize(original.m);
    this->n      = original.n;
    this->m      = original.m + original.n;
    this->f      = [this](crvec x) { return original.f(x); };
    this->grad_f = [this](crvec x, rvec grad) { original.grad_f(x, grad); };
    this->g      = [this](crvec x, rvec gg) {
        original.g(x, work);
        gg.topRows(original.m)    = work;
        gg.bottomRows(original.n) = x;
    };
    this->grad_g_prod = [this](crvec x, crvec y, rvec gg) {
        work = y.topRows(original.m);
        original.grad_g_prod(x, work, gg);
        gg += y.bottomRows(original.n);
    };
    this->grad_gi = [](crvec, unsigned, rvec) {
        throw std::logic_error("ProblemOnlyD::grad_gi: Not implemented");
    };
    this->hess_L_prod = [](crvec, crvec, crvec, rvec) {
        throw std::logic_error("ProblemOnlyD::hess_L_prod: Not implemented");
    };
    this->hess_L = [](crvec, crvec, rmat) {
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