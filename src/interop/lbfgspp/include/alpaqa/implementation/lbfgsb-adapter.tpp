#pragma once

#include <alpaqa/lbfgsb-adapter.hpp>
#include <cmath>
#include <exception>

namespace alpaqa::lbfgspp {

template <Config Conf>
std::string LBFGSBSolver<Conf>::get_name() const {
    return "LBFGSBSolver";
}

template <Config Conf>
auto LBFGSBSolver<Conf>::operator()(
    /// [in]    Problem description
    const Problem &problem,
    /// [in]    Solve options
    const SolveOptions &opts,
    /// [inout] Decision variable @f$ x @f$
    rvec x,
    /// [inout] Lagrange multipliers @f$ y @f$
    rvec y,
    /// [in]    Constraint weights @f$ \Sigma @f$
    crvec Σ,
    /// [out]   Slack variable error @f$ g(x) - \Pi_D(g(x) + \Sigma^{-1} y) @f$
    rvec err_z) -> Stats {

    using std::chrono::nanoseconds;
    auto start_time = std::chrono::steady_clock::now();
    auto os         = opts.os ? opts.os : this->os;
    auto max_time   = nanoseconds::max();
    if (opts.max_time)
        max_time = std::min(max_time, *opts.max_time);
    Stats s;

    const auto n  = problem.get_n();
    const auto m  = problem.get_m();
    const auto &C = problem.get_box_C();

    vec work_n(n), work_m(m);
    vec x_solve = x;

    struct BreakException {
        SolverStatus status = SolverStatus::Busy;
    };

    ::LBFGSpp::LBFGSBParam<real_t> effective_params = params;
    effective_params.epsilon                        = opts.tolerance;
    effective_params.epsilon_rel                    = 0;
    ::LBFGSpp::LBFGSBSolver<real_t> solver{effective_params};

    // Evaluate cost and its gradient, checking for termination
    auto eval_f_grad_f = [&](crvec xk, rvec grad) {
        // Check user interrupt
        if (stop_signal.stop_requested()) {
            x_solve = xk;
            throw BreakException{SolverStatus::Interrupted};
        }
        // Check maximum run time
        auto time_elapsed = std::chrono::steady_clock::now() - start_time;
        bool out_of_time  = time_elapsed > max_time;
        if (out_of_time) {
            if (opts.always_overwrite_results)
                x_solve = xk;
            else
                s.final_ψ = problem.eval_ψ(xk, y, Σ, work_m);
            throw BreakException{SolverStatus::MaxTime};
        }
        // Perform the actual evaluation
        const auto ψx = problem.eval_ψ_grad_ψ(xk, y, Σ, grad, work_n, work_m);
        // Check that the function value is finite
        if (!std::isfinite(ψx)) {
            if (opts.always_overwrite_results)
                x_solve = xk;
            s.final_ψ = ψx;
            throw BreakException{SolverStatus::NotFinite};
        }
        return ψx;
    };

    // Solve problem
    try {
        s.iterations = solver.minimize(eval_f_grad_f, x_solve, s.final_ψ,
                                       C.lowerbound, C.upperbound);
        s.status     = SolverStatus::Converged;
        if (static_cast<int>(s.iterations) == effective_params.max_iterations)
            s.status = SolverStatus::MaxIter;
    } catch (const BreakException &e) {
        s.status = e.status;
    } catch (const std::exception &e) {
        s.status = SolverStatus::Exception;
        *os << "[LBFGSB] Exception: " << e.what() << std::endl;
    }
    auto time_elapsed = std::chrono::steady_clock::now() - start_time;
    s.elapsed_time    = duration_cast<nanoseconds>(time_elapsed);

    // Check final error
    s.ε = solver.final_grad_norm();
    // Update result vectors
    if (s.status == SolverStatus::Converged ||
        s.status == SolverStatus::Interrupted ||
        opts.always_overwrite_results) {
        auto &ŷ   = work_m;
        s.final_ψ = problem.eval_ψ(x_solve, y, Σ, ŷ);
        if (err_z.size() > 0)
            err_z = Σ.asDiagonal().inverse() * (ŷ - y);
        x = x_solve;
        y = ŷ;
    }
    return s;
}

} // namespace alpaqa::lbfgspp
