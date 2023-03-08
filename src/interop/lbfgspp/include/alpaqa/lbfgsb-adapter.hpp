#pragma once

#include <alpaqa/config/config.hpp>
#include <alpaqa/inner/inner-solve-options.hpp>
#include <alpaqa/inner/internal/solverstatus.hpp>
#include <alpaqa/lbfgspp-adapter-export.hpp>
#include <alpaqa/problem/type-erased-problem.hpp>
#include <alpaqa/util/atomic-stop-signal.hpp>

#include <chrono>
#include <iostream>

#include <LBFGSB.h>

namespace alpaqa::lbfgspp {

template <Config Conf = DefaultConfig>
struct LBFGSBStats {
    USING_ALPAQA_CONFIG(Conf);

    SolverStatus status = SolverStatus::Busy;
    real_t ε            = inf<config_t>;
    std::chrono::nanoseconds elapsed_time{};
    unsigned iterations = 0;
    real_t final_ψ      = 0;
};

/// L-BFGS-B solver for ALM.
/// @ingroup    grp_InnerSolvers
template <Config Conf>
class LBFGSBSolver {
  public:
    USING_ALPAQA_CONFIG(Conf);

    using Problem      = TypeErasedProblem<config_t>;
    using Params       = ::LBFGSpp::LBFGSBParam<real_t>;
    using Stats        = LBFGSBStats<config_t>;
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

ALPAQA_LBFGSPP_EXPORT_EXTERN_TEMPLATE(struct, LBFGSBStats, DefaultConfig);
ALPAQA_LBFGSPP_EXPORT_EXTERN_TEMPLATE(struct, LBFGSBStats, EigenConfigf);
ALPAQA_LBFGSPP_EXPORT_EXTERN_TEMPLATE(struct, LBFGSBStats, EigenConfigd);
ALPAQA_LBFGSPP_EXPORT_EXTERN_TEMPLATE(struct, LBFGSBStats, EigenConfigl);
#ifdef ALPAQA_WITH_QUAD_PRECISION
ALPAQA_LBFGSPP_EXPORT_EXTERN_TEMPLATE(struct, LBFGSBStats, EigenConfigq);
#endif

ALPAQA_LBFGSPP_EXPORT_EXTERN_TEMPLATE(class, LBFGSBSolver, DefaultConfig);
ALPAQA_LBFGSPP_EXPORT_EXTERN_TEMPLATE(class, LBFGSBSolver, EigenConfigf);
ALPAQA_LBFGSPP_EXPORT_EXTERN_TEMPLATE(class, LBFGSBSolver, EigenConfigd);
ALPAQA_LBFGSPP_EXPORT_EXTERN_TEMPLATE(class, LBFGSBSolver, EigenConfigl);
#ifdef ALPAQA_WITH_QUAD_PRECISION
ALPAQA_LBFGSPP_EXPORT_EXTERN_TEMPLATE(class, LBFGSBSolver, EigenConfigq);
#endif

} // namespace alpaqa::lbfgspp

namespace alpaqa {

template <class InnerSolverStats>
struct InnerStatsAccumulator;

template <Config Conf>
struct InnerStatsAccumulator<lbfgspp::LBFGSBStats<Conf>> {
    USING_ALPAQA_CONFIG(Conf);

    /// Total elapsed time in the inner solver.
    std::chrono::nanoseconds elapsed_time{};
    /// Total number of inner PANOC iterations.
    unsigned iterations = 0;
    /// Final value of the smooth cost @f$ \psi(\hat x) @f$.
    real_t final_ψ = 0;
};

template <Config Conf>
InnerStatsAccumulator<lbfgspp::LBFGSBStats<Conf>> &
operator+=(InnerStatsAccumulator<lbfgspp::LBFGSBStats<Conf>> &acc,
           const lbfgspp::LBFGSBStats<Conf> &s) {
    acc.iterations += s.iterations;
    acc.elapsed_time += s.elapsed_time;
    acc.final_ψ = s.final_ψ;
    return acc;
}

} // namespace alpaqa