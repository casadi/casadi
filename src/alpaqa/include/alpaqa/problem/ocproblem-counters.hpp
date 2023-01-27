#pragma once

#include <alpaqa/export.h>

#include <chrono>
#include <iosfwd>

namespace alpaqa {

struct OCPEvalCounter {
    unsigned f{};
    unsigned jac_f{};
    unsigned grad_f_prod{};
    unsigned h{};
    unsigned h_N{};
    unsigned l{};
    unsigned l_N{};
    unsigned qr{};
    unsigned q_N{};
    unsigned add_Q{};
    unsigned add_Q_N{};
    unsigned add_R_masked{};
    unsigned add_S_masked{};
    unsigned add_R_prod_masked{};
    unsigned add_S_prod_masked{};
    unsigned constr{};
    unsigned constr_N{};
    unsigned grad_constr_prod{};
    unsigned grad_constr_prod_N{};
    unsigned add_gn_hess_constr{};
    unsigned add_gn_hess_constr_N{};

    struct OCPEvalTimer {
        std::chrono::nanoseconds f{};
        std::chrono::nanoseconds jac_f{};
        std::chrono::nanoseconds grad_f_prod{};
        std::chrono::nanoseconds h{};
        std::chrono::nanoseconds h_N{};
        std::chrono::nanoseconds l{};
        std::chrono::nanoseconds l_N{};
        std::chrono::nanoseconds qr{};
        std::chrono::nanoseconds q_N{};
        std::chrono::nanoseconds add_Q{};
        std::chrono::nanoseconds add_Q_N{};
        std::chrono::nanoseconds add_R_masked{};
        std::chrono::nanoseconds add_S_masked{};
        std::chrono::nanoseconds add_R_prod_masked{};
        std::chrono::nanoseconds add_S_prod_masked{};
        std::chrono::nanoseconds constr{};
        std::chrono::nanoseconds constr_N{};
        std::chrono::nanoseconds grad_constr_prod{};
        std::chrono::nanoseconds grad_constr_prod_N{};
        std::chrono::nanoseconds add_gn_hess_constr{};
        std::chrono::nanoseconds add_gn_hess_constr_N{};
    } time;

    void reset() { *this = {}; }
};

ALPAQA_EXPORT std::ostream &operator<<(std::ostream &, const OCPEvalCounter &);

inline OCPEvalCounter::OCPEvalTimer &operator+=(OCPEvalCounter::OCPEvalTimer &a,
                                                const OCPEvalCounter::OCPEvalTimer &b) {
    a.f += b.f;
    a.jac_f += b.jac_f;
    a.grad_f_prod += b.grad_f_prod;
    a.h += b.h;
    a.h_N += b.h_N;
    a.l += b.l;
    a.l_N += b.l_N;
    a.qr += b.qr;
    a.q_N += b.q_N;
    a.add_Q += b.add_Q;
    a.add_Q_N += b.add_Q_N;
    a.add_R_masked += b.add_R_masked;
    a.add_S_masked += b.add_S_masked;
    a.add_R_prod_masked += b.add_R_prod_masked;
    a.add_S_prod_masked += b.add_S_prod_masked;
    a.constr += b.constr;
    a.constr_N += b.constr_N;
    a.grad_constr_prod += b.grad_constr_prod;
    a.grad_constr_prod_N += b.grad_constr_prod_N;
    a.add_gn_hess_constr += b.add_gn_hess_constr;
    a.add_gn_hess_constr_N += b.add_gn_hess_constr_N;
    return a;
}

inline OCPEvalCounter &operator+=(OCPEvalCounter &a, const OCPEvalCounter &b) {
    a.f += b.f;
    a.jac_f += b.jac_f;
    a.grad_f_prod += b.grad_f_prod;
    a.h += b.h;
    a.h_N += b.h_N;
    a.l += b.l;
    a.l_N += b.l_N;
    a.qr += b.qr;
    a.q_N += b.q_N;
    a.add_Q += b.add_Q;
    a.add_Q_N += b.add_Q_N;
    a.add_R_masked += b.add_R_masked;
    a.add_S_masked += b.add_S_masked;
    a.add_R_prod_masked += b.add_R_prod_masked;
    a.add_S_prod_masked += b.add_S_prod_masked;
    a.constr += b.constr;
    a.constr_N += b.constr_N;
    a.grad_constr_prod += b.grad_constr_prod;
    a.grad_constr_prod_N += b.grad_constr_prod_N;
    a.add_gn_hess_constr += b.add_gn_hess_constr;
    a.add_gn_hess_constr_N += b.add_gn_hess_constr_N;
    a.time += b.time;
    return a;
}

inline OCPEvalCounter operator+(OCPEvalCounter a, const OCPEvalCounter &b) { return a += b; }

} // namespace alpaqa
