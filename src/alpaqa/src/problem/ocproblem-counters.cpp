#include <alpaqa/problem/ocproblem-counters.hpp>

#include <iomanip>
#include <iostream>

namespace alpaqa {

namespace {
struct CountResult {
    unsigned count;
    std::chrono::nanoseconds time;
};
std::ostream &operator<<(std::ostream &os, const CountResult &t) {
    auto sec = [](auto t) { return std::chrono::duration<double>(t).count(); };
    os << std::setw(8);
    if (t.count > 0) {
        os << t.count << "  (";
        auto old  = os.flags();
        auto prec = os.precision(3);
        os << std::scientific << std::setw(9) << 1e6 * sec(t.time) << " µs, "
           << std::setw(9) << 1e6 * sec(t.time) / static_cast<double>(t.count)
           << " µs/call)\r\n";
        os.precision(prec);
        os.flags(old);
    } else {
        os << '-' << "\r\n";
    }
    return os;
}
} // namespace

std::ostream &operator<<(std::ostream &os, const OCPEvalCounter &c) {
    os << "                   f:" //
       << CountResult{c.f, c.time.f};
    os << "               jac_f:" //
       << CountResult{c.jac_f, c.time.jac_f};
    os << "         grad_f_prod:" //
       << CountResult{c.grad_f_prod, c.time.grad_f_prod};
    os << "                   h:" //
       << CountResult{c.h, c.time.h};
    os << "                 h_N:" //
       << CountResult{c.h_N, c.time.h_N};
    os << "                   l:" //
       << CountResult{c.l, c.time.l};
    os << "                 l_N:" //
       << CountResult{c.l_N, c.time.l_N};
    os << "                  qr:" //
       << CountResult{c.qr, c.time.qr};
    os << "                 q_N:" //
       << CountResult{c.q_N, c.time.q_N};
    os << "               add_Q:" //
       << CountResult{c.add_Q, c.time.add_Q};
    os << "             add_Q_N:" //
       << CountResult{c.add_Q_N, c.time.add_Q_N};
    os << "        add_R_masked:" //
       << CountResult{c.add_R_masked, c.time.add_R_masked};
    os << "        add_S_masked:" //
       << CountResult{c.add_S_masked, c.time.add_S_masked};
    os << "   add_R_prod_masked:" //
       << CountResult{c.add_R_prod_masked, c.time.add_R_prod_masked};
    os << "   add_S_prod_masked:" //
       << CountResult{c.add_S_prod_masked, c.time.add_S_prod_masked};
    os << "              constr:" //
       << CountResult{c.constr, c.time.constr};
    os << "            constr_N:" //
       << CountResult{c.constr_N, c.time.constr_N};
    os << "    grad_constr_prod:" //
       << CountResult{c.grad_constr_prod, c.time.grad_constr_prod};
    os << "  grad_constr_prod_N:" //
       << CountResult{c.grad_constr_prod_N, c.time.grad_constr_prod_N};
    os << "  add_gn_hess_constr:" //
       << CountResult{c.add_gn_hess_constr, c.time.add_gn_hess_constr};
    os << "add_gn_hess_constr_N:" //
       << CountResult{c.add_gn_hess_constr_N, c.time.add_gn_hess_constr_N};
    return os;
}

} // namespace alpaqa