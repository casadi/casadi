#include <alpaqa/problem/problem-counters.hpp>

#include <iomanip>
#include <iostream>

namespace alpaqa {

std::ostream &operator<<(std::ostream &os, const EvalCounter &c) {
    auto cnt = [](auto t) { return std::chrono::duration<double>(t).count(); };
    os << "                 f:" << std::setw(6) << c.f << "  (" << cnt(c.time.f)
       << " s)\r\n";
    os << "            grad_f:" << std::setw(6) << c.grad_f << "  ("
       << cnt(c.time.grad_f) << " s)\r\n";
    os << "          f_grad_f:" << std::setw(6) << c.f_grad_f << "  ("
       << cnt(c.time.f_grad_f) << " s)\r\n";
    os << "               f_g:" << std::setw(6) << c.f_g << "  ("
       << cnt(c.time.f_g) << " s)\r\n";
    os << "        f_grad_f_g:" << std::setw(6) << c.f_grad_f_g << "  ("
       << cnt(c.time.f_grad_f_g) << " s)\r\n";
    os << "grad_f_grad_g_prod:" << std::setw(6) << c.grad_f_grad_g_prod << "  ("
       << cnt(c.time.grad_f_grad_g_prod) << " s)\r\n";
    os << "                 g:" << std::setw(6) << c.g << "  (" << cnt(c.time.g)
       << " s)\r\n";
    os << "       grad_g_prod:" << std::setw(6) << c.grad_g_prod << "  ("
       << cnt(c.time.grad_g_prod) << " s)\r\n";
    os << "           grad_gi:" << std::setw(6) << c.grad_gi << "  ("
       << cnt(c.time.grad_gi) << " s)\r\n";
    os << "            grad_L:" << std::setw(6) << c.grad_L << "  ("
       << cnt(c.time.grad_L) << " s)\r\n";
    os << "       hess_L_prod:" << std::setw(6) << c.hess_L_prod << "  ("
       << cnt(c.time.hess_L_prod) << " s)\r\n";
    os << "            hess_L:" << std::setw(6) << c.hess_L << "  ("
       << cnt(c.time.hess_L) << " s)\r\n";
    os << "                 ψ:" << std::setw(6) << c.ψ << "  (" << cnt(c.time.ψ)
       << " s)\r\n";
    os << "            grad_ψ:" << std::setw(6) << c.grad_ψ << "  ("
       << cnt(c.time.grad_ψ) << " s)\r\n";
    os << "     grad_ψ_from_ŷ:" << std::setw(6) << c.grad_ψ_from_ŷ << "  ("
       << cnt(c.time.grad_ψ_from_ŷ) << " s)\r\n";
    os << "          ψ_grad_ψ:" << std::setw(6) << c.ψ_grad_ψ << "  ("
       << cnt(c.time.ψ_grad_ψ) << " s)";
    return os;
}

} // namespace alpaqa