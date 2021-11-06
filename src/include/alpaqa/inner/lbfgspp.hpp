#pragma once

#include "alpaqa/util/atomic_stop_signal.hpp"
#include <alpaqa/inner/detail/panoc-helpers.hpp>
#include <alpaqa/util/problem.hpp>
#include <alpaqa/util/solverstatus.hpp>

#include <atomic>
#include <chrono>
#include <stdexcept>
#include <string>

#include <LBFGS.h>
#include <LBFGSB.h>

namespace alpaqa {

using std::chrono::duration_cast;
using std::chrono::microseconds;

/// Unconstrained LBFGS solver for ALM.
/// @ingroup    grp_InnerSolvers
template <template <class> class LineSearchT = LBFGSpp::LineSearchBacktracking>
class LBFGSSolver {
  public:
    using Params     = LBFGSpp::LBFGSParam<real_t>;
    using LineSearch = LineSearchT<real_t>;

    struct Stats {
        unsigned iterations = 0;
        real_t ε            = inf;
        std::chrono::microseconds elapsed_time;
        SolverStatus status = SolverStatus::Unknown;

        // TODO: find more generic way of handling this
        unsigned linesearch_failures = 0;
        unsigned lbfgs_failures      = 0; // TODO: more generic name
        unsigned lbfgs_rejected      = 0;
    };

    LBFGSSolver(Params params) : params(params) {}

    Stats operator()(const Problem &problem,        // in
                     crvec Σ,                       // in
                     real_t ε,                      // in
                     bool always_overwrite_results, // in
                     rvec x,                        // inout
                     rvec y,                        // inout
                     rvec err_z                     // out
    ) {
        Stats s;
        auto start_time = std::chrono::steady_clock::now();
        work_n.resize(problem.n);
        work_m.resize(problem.m);
        auto calc_ψ_grad_ψ = [&](crvec x, rvec grad) {
            return detail::calc_ψ_grad_ψ(problem, x, y, Σ, // in
                                         grad,             // out
                                         work_n, work_m);  // work
        };
        real_t ψ٭;
        params.epsilon = ε;
        vec x٭         = x;
        try {
            LBFGSpp::LBFGSSolver<double> solver(params);
            s.iterations = solver.minimize(calc_ψ_grad_ψ, x٭, ψ٭);
            s.status     = SolverStatus::Converged; // TODO: <?>
        } catch (std::runtime_error &) {
            s.status = SolverStatus::MaxIter;
        } catch (std::logic_error &) {
            // TODO: when does the following error occur?
            // “the moving direction increases the objective function value”
            s.status = SolverStatus::NotFinite; // <?>
        }
        auto time_elapsed = std::chrono::steady_clock::now() - start_time;
        s.elapsed_time    = duration_cast<microseconds>(time_elapsed);
        s.ε               = ε; // TODO: <?>
        if (stop_signal.stop_requested())
            s.status = SolverStatus::Interrupted;
        bool conv        = s.status == SolverStatus::Converged;
        bool interrupted = s.status == SolverStatus::Interrupted;
        if (conv || interrupted || always_overwrite_results) {
            x٭.swap(x);
            work_m.swap(y);
            // TODO: could be optimized
            (void)detail::calc_ψ_ŷ(problem, x, work_m, Σ, y);
            detail::calc_err_z(problem, x, work_m, Σ, err_z);
        }
        return s;
    }

    std::string get_name() const { return "LBFGSSolver<?>"; }

    const Params &get_params() const { return params; }

    void stop() { stop_signal.stop(); }

  private:
    Params params;
    AtomicStopSignal stop_signal;
    vec work_n, work_m;
};

struct LBFGSBStats {
    SolverStatus status = SolverStatus::Unknown;
    real_t ε            = inf;
    std::chrono::microseconds elapsed_time;
    unsigned iterations = 0;
};

/// Box-constrained LBFGS solver for ALM.
/// @ingroup    grp_InnerSolvers
template <template <class> class LineSearchT = LBFGSpp::LineSearchBacktracking>
class LBFGSBSolver {
  public:
    using Params     = LBFGSpp::LBFGSBParam<real_t>;
    using LineSearch = LineSearchT<real_t>;
    using Stats      = LBFGSBStats;

    LBFGSBSolver(Params params) : params(params) {}

    Stats operator()(const Problem &problem,        // in
                     crvec Σ,                       // in
                     real_t ε,                      // in
                     bool always_overwrite_results, // in
                     rvec x,                        // inout
                     rvec y,                        // inout
                     rvec err_z                     // out
    ) {
        Stats s;
        auto start_time = std::chrono::steady_clock::now();
        work_n.resize(problem.n);
        work_m.resize(problem.m);
        auto calc_ψ_grad_ψ = [&](crvec x, rvec grad) {
            return detail::calc_ψ_grad_ψ(problem, x, y, Σ, // in
                                         grad,             // out
                                         work_n, work_m);  // work
        };
        real_t ψ٭;
        params.epsilon = ε;
        vec x٭         = x;
        try {
            LBFGSpp::LBFGSBSolver<double> solver(params);
            s.iterations =
                solver.minimize(calc_ψ_grad_ψ, x٭, ψ٭, problem.C.lowerbound,
                                problem.C.upperbound);
            s.status = SolverStatus::Converged; // TODO: <?>
        } catch (std::runtime_error &) {
            s.status = SolverStatus::MaxIter;
        } catch (std::logic_error &) {
            // TODO: when does the following error occur?
            // “the moving direction increases the objective function value”
            s.status = SolverStatus::NotFinite; // <?>
        }
        auto time_elapsed = std::chrono::steady_clock::now() - start_time;
        s.elapsed_time    = duration_cast<microseconds>(time_elapsed);
        s.ε               = ε; // TODO: <?>
        if (stop_signal.stop_requested())
            s.status = SolverStatus::Interrupted;
        bool conv        = s.status == SolverStatus::Converged;
        bool interrupted = s.status == SolverStatus::Interrupted;
        if (conv || interrupted || always_overwrite_results) {
            x٭.swap(x);
            work_m.swap(y);
            // TODO: could be optimized
            (void)detail::calc_ψ_ŷ(problem, x, work_m, Σ, y);
            detail::calc_err_z(problem, x, work_m, Σ, err_z);
        }
        return s;
    }

    std::string get_name() const { return "LBFGSBSolver<?>"; }

    const Params &get_params() const { return params; }

    void stop() { stop_signal.stop(); }

  private:
    Params params;
    AtomicStopSignal stop_signal;
    vec work_n, work_m;
};

template <class InnerSolverStats>
struct InnerStatsAccumulator;

template <>
struct InnerStatsAccumulator<LBFGSBStats> {
    std::chrono::microseconds elapsed_time;
    unsigned iterations = 0;
};

inline InnerStatsAccumulator<LBFGSBStats> &
operator+=(InnerStatsAccumulator<LBFGSBStats> &acc, const LBFGSBStats &s) {
    acc.iterations += s.iterations;
    acc.elapsed_time += s.elapsed_time;
    return acc;
}

} // namespace alpaqa
