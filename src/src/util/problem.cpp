#include <alpaqa/util/problem.hpp>

namespace alpaqa {

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

} // namespace alpaqa