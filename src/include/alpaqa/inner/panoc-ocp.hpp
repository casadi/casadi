#pragma once

#include <alpaqa/inner/panoc.hpp>
#include "alpaqa/problem/dynamics.hpp"

#include <chrono>
#include <limits>
#include <string>

namespace alpaqa {

template <Config Conf = DefaultConfig>
struct PANOCOCPProgressInfo {
    USING_ALPAQA_CONFIG(Conf);

    unsigned k;
    crvec xu;
    crvec p;
    real_t norm_sq_p;
    crvec x̂u;
    real_t φγ;
    real_t ψ;
    crvec grad_ψ;
    real_t ψ_hat;
    crvec q;
    real_t L;
    real_t γ;
    real_t τ;
    real_t ε;
    const TypeErasedControlProblem<config_t> &problem;
    const PANOCParams<config_t> &params;
};

template <Config Conf>
class PANOCOCPSolver {
  public:
    USING_ALPAQA_CONFIG(Conf);

    using Problem      = alpaqa::TypeErasedControlProblem<config_t>;
    using Params       = PANOCParams<config_t>;
    using Stats        = PANOCStats<config_t>;
    using ProgressInfo = PANOCOCPProgressInfo<config_t>;

    PANOCOCPSolver(const Params &params) : params(params) {}

    Stats operator()(const Problem &problem, // in
                     real_t ε,               // in
                     rvec u);                // inout

    /// Specify a callable that is invoked with some intermediate results on
    /// each iteration of the algorithm.
    /// @see @ref ProgressInfo
    PANOCOCPSolver &
    set_progress_callback(std::function<void(const ProgressInfo &)> cb) {
        this->progress_cb = cb;
        return *this;
    }

    std::string get_name() const;

    void stop() { stop_signal.stop(); }

    const Params &get_params() const { return params; }

  private:
    Params params;
    AtomicStopSignal stop_signal;
    std::function<void(const ProgressInfo &)> progress_cb;
    using Helpers = detail::PANOCHelpers<config_t>;
};

ALPAQA_EXPORT_EXTERN_TEMPLATE(struct, PANOCOCPProgressInfo, DefaultConfig);
ALPAQA_EXPORT_EXTERN_TEMPLATE(struct, PANOCOCPProgressInfo, EigenConfigf);
ALPAQA_EXPORT_EXTERN_TEMPLATE(struct, PANOCOCPProgressInfo, EigenConfigd);
ALPAQA_EXPORT_EXTERN_TEMPLATE(struct, PANOCOCPProgressInfo, EigenConfigl);
#ifdef ALPAQA_WITH_QUAD_PRECISION
ALPAQA_EXPORT_EXTERN_TEMPLATE(struct, PANOCOCPProgressInfo, EigenConfigq);
#endif

ALPAQA_EXPORT_EXTERN_TEMPLATE(class, PANOCOCPSolver, DefaultConfig);
ALPAQA_EXPORT_EXTERN_TEMPLATE(class, PANOCOCPSolver, EigenConfigf);
ALPAQA_EXPORT_EXTERN_TEMPLATE(class, PANOCOCPSolver, EigenConfigd);
ALPAQA_EXPORT_EXTERN_TEMPLATE(class, PANOCOCPSolver, EigenConfigl);
#ifdef ALPAQA_WITH_QUAD_PRECISION
ALPAQA_EXPORT_EXTERN_TEMPLATE(class, PANOCOCPSolver, EigenConfigq);
#endif

} // namespace alpaqa
