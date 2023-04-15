#pragma once

#include <atomic>
#include <csignal>
#include <memory>
#include <stdexcept>

/**
 * Attach SIGINT and SIGTERM handlers to stop the given solver.
 * @tparam  solver_to_stop
 *          Reference to atomic pointer variable that will store a pointer to
 *          the solver, should be static, so it can be accessed from the signal
 *          handler without causing lifetime issues.
 * @param   solver
 *          The solver that should be stopped by the handler.
 * @return  A RAII object that detaches the handler when destroyed.
 */
template <auto &solver_to_stop>
auto attach_cancellation(auto &solver) {
    if constexpr (requires { solver.stop(); }) {
        auto *old = solver_to_stop.exchange(&solver, std::memory_order_release);
        if (old) {
            old->stop();
            throw std::runtime_error(
                "alpaqa-driver:attach_cancellation is not reentrant");
        }
        struct sigaction action;
        action.sa_handler = [](int) {
            if (auto *s = solver_to_stop.load(std::memory_order::acquire))
                s->stop();
        };
        sigemptyset(&action.sa_mask);
        action.sa_flags = 0;
        sigaction(SIGINT, &action, nullptr);
        sigaction(SIGTERM, &action, nullptr);
    }
    using solver_to_stop_t = std::remove_reference_t<decltype(solver_to_stop)>;
    auto detach_solver     = [](solver_to_stop_t *p) {
        struct sigaction action;
        action.sa_handler = SIG_DFL;
        sigemptyset(&action.sa_mask);
        action.sa_flags = 0;
        sigaction(SIGINT, &action, nullptr);
        sigaction(SIGTERM, &action, nullptr);
        p->store(nullptr, std::memory_order_relaxed);
    };
    return std::unique_ptr<solver_to_stop_t, decltype(detach_solver)>{
        &solver_to_stop};
}
