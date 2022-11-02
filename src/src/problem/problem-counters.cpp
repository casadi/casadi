#include <alpaqa/problem/problem-counters.hpp>

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
        os << t.count << "  (" << sec(t.time) << " s)\r\n";
    } else {
        os << '-' << "\r\n";
    }
    return os;
}
} // namespace

std::ostream &operator<<(std::ostream &os, const EvalCounter &c) {
    os << "        proj_diff_g:" << std::setw(6) //
       << CountResult{c.proj_diff_g, c.time.proj_diff_g};
    os << "   proj_multipliers:" << std::setw(6) //
       << CountResult{c.proj_multipliers, c.time.proj_multipliers};
    os << "     prox_grad_step:" << std::setw(6) //
       << CountResult{c.prox_grad_step, c.time.prox_grad_step};
    os << "                  f:" << std::setw(6) //
       << CountResult{c.f, c.time.f};
    os << "             grad_f:" << std::setw(6) //
       << CountResult{c.grad_f, c.time.grad_f};
    os << "           f_grad_f:" << std::setw(6) //
       << CountResult{c.f_grad_f, c.time.f_grad_f};
    os << "                f_g:" << std::setw(6) //
       << CountResult{c.f_g, c.time.f_g};
    os << "         f_grad_f_g:" << std::setw(6) //
       << CountResult{c.f_grad_f_g, c.time.f_grad_f_g};
    os << " grad_f_grad_g_prod:" << std::setw(6) //
       << CountResult{c.grad_f_grad_g_prod, c.time.grad_f_grad_g_prod};
    os << "                  g:" << std::setw(6) //
       << CountResult{c.g, c.time.g};
    os << "        grad_g_prod:" << std::setw(6) //
       << CountResult{c.grad_g_prod, c.time.grad_g_prod};
    os << "            grad_gi:" << std::setw(6) //
       << CountResult{c.grad_gi, c.time.grad_gi};
    os << "             grad_L:" << std::setw(6) //
       << CountResult{c.grad_L, c.time.grad_L};
    os << "        hess_L_prod:" << std::setw(6) //
       << CountResult{c.hess_L_prod, c.time.hess_L_prod};
    os << "             hess_L:" << std::setw(6) //
       << CountResult{c.hess_L, c.time.hess_L};
    os << "                  ψ:" << std::setw(6) //
       << CountResult{c.ψ, c.time.ψ};
    os << "             grad_ψ:" << std::setw(6) //
       << CountResult{c.grad_ψ, c.time.grad_ψ};
    os << "      grad_ψ_from_ŷ:" << std::setw(6) //
       << CountResult{c.grad_ψ_from_ŷ, c.time.grad_ψ_from_ŷ};
    os << "           ψ_grad_ψ:" << std::setw(6) //
       << CountResult{c.ψ_grad_ψ, c.time.ψ_grad_ψ};
    return os;
}

} // namespace alpaqa