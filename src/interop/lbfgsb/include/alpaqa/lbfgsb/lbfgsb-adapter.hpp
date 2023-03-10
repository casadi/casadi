#pragma once

#include <alpaqa/config/config.hpp>
#include <alpaqa/inner/inner-solve-options.hpp>
#include <alpaqa/inner/internal/panoc-stop-crit.hpp>
#include <alpaqa/inner/internal/solverstatus.hpp>
#include <alpaqa/lbfgsb-adapter-export.h>
#include <alpaqa/problem/type-erased-problem.hpp>
#include <alpaqa/util/atomic-stop-signal.hpp>

#include <chrono>
#include <iostream>

namespace alpaqa::lbfgsb {

struct LBFGSB_ADAPTER_EXPORT LBFGSBParams {
    USING_ALPAQA_CONFIG(alpaqa::DefaultConfig);

    unsigned memory                   = 10;
    unsigned max_iter                 = 1000;
    std::chrono::nanoseconds max_time = std::chrono::minutes(5);
    PANOCStopCrit stop_crit           = PANOCStopCrit::ProjGradUnitNorm;
    int print                         = -1;
    unsigned print_interval           = 0;
    int print_precision = std::numeric_limits<real_t>::max_digits10 / 2;
};

struct LBFGSB_ADAPTER_EXPORT LBFGSBStats {
    USING_ALPAQA_CONFIG(alpaqa::DefaultConfig);

    SolverStatus status = SolverStatus::Busy;
    real_t ε            = inf<config_t>;
    std::chrono::nanoseconds elapsed_time{};
    unsigned iterations     = 0;
    real_t final_ψ          = 0;
    unsigned lbfgs_rejected = 0;
};

/// L-BFGS-B solver for ALM.
/// @ingroup    grp_InnerSolvers
class LBFGSB_ADAPTER_EXPORT LBFGSBSolver {
  public:
    USING_ALPAQA_CONFIG(alpaqa::DefaultConfig);

    using Problem      = TypeErasedProblem<config_t>;
    using Params       = LBFGSBParams;
    using Stats        = LBFGSBStats;
    using SolveOptions = InnerSolveOptions<config_t>;

    LBFGSBSolver(const Params &params) : params(params) {}

    Stats operator()(const Problem &problem,   // in
                     const SolveOptions &opts, // in
                     rvec x,                   // inout
                     rvec y,                   // inout
                     crvec Σ,                  // in
                     rvec err_z);              // out

    template <class P>
    Stats operator()(const P &problem, const SolveOptions &opts, rvec u, rvec y,
                     crvec Σ, rvec e) {
        return operator()(Problem::template make<P>(problem), opts, u, y, Σ, e);
    }

    std::string get_name() const;

    void stop() { stop_signal.stop(); }

    const Params &get_params() const { return params; }

  private:
    Params params;
    AtomicStopSignal stop_signal;

  public:
    std::ostream *os = &std::cout;
};

} // namespace alpaqa::lbfgsb

namespace alpaqa {

template <class InnerSolverStats>
struct InnerStatsAccumulator;

template <>
struct InnerStatsAccumulator<lbfgsb::LBFGSBStats> {
    USING_ALPAQA_CONFIG(alpaqa::DefaultConfig);

    /// Total elapsed time in the inner solver.
    std::chrono::nanoseconds elapsed_time{};
    /// Total number of inner PANOC iterations.
    unsigned iterations = 0;
    /// Final value of the smooth cost @f$ \psi(\hat x) @f$.
    real_t final_ψ = 0;
    /// Total number of times that the L-BFGS update was rejected (i.e. it
    /// could have resulted in a non-positive definite Hessian estimate).
    unsigned lbfgs_rejected = 0;
};

inline InnerStatsAccumulator<lbfgsb::LBFGSBStats> &
operator+=(InnerStatsAccumulator<lbfgsb::LBFGSBStats> &acc,
           const lbfgsb::LBFGSBStats &s) {
    acc.iterations += s.iterations;
    acc.elapsed_time += s.elapsed_time;
    acc.final_ψ = s.final_ψ;
    acc.lbfgs_rejected += s.lbfgs_rejected;
    return acc;
}

} // namespace alpaqa