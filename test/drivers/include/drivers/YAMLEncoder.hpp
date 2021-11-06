#include <alpaqa/inner/decl/panoc-stop-crit.hpp>
#include <alpaqa/inner/decl/panoc.hpp>
#include <alpaqa/inner/decl/second-order-panoc.hpp>
#include <alpaqa/inner/decl/structured-panoc-lbfgs.hpp>
#include <alpaqa/inner/guarded-aa-pga.hpp>
#include <alpaqa/inner/lbfgspp.hpp>
#include <alpaqa/inner/pga.hpp>

#include <alpaqa/decl/alm.hpp>
#include <alpaqa/interop/cutest/CUTEstLoader.hpp>
#include <alpaqa/util/solverstatus.hpp>

#include <yaml-cpp/emitter.h>
#include <yaml-cpp/emittermanip.h>
#include <yaml-cpp/yaml.h>

#include <sstream>

inline YAML::Emitter &operator<<(YAML::Emitter &out, alpaqa::crvec v) {
    out << YAML::Flow;
    out << YAML::BeginSeq;
    for (alpaqa::vec::Index i = 0; i < v.size(); ++i)
        out << v[i];
    out << YAML::EndSeq;
    return out;
}

inline YAML::Emitter &operator<<(YAML::Emitter &out,
                                 CUTEstProblem::Report::Status s) {
    return out << enum_name(s);
}

inline YAML::Emitter &operator<<(YAML::Emitter &out,
                                 const CUTEstProblem::Report &r) {
    out << YAML::BeginMap;
    out << YAML::Key << r.name;
    out << YAML::Value << YAML::BeginMap;
    out << YAML::Key << "nvar" << YAML::Value << r.nvar;
    out << YAML::Key << "ncon" << YAML::Value << r.ncon;
    out << YAML::Key << "status" << YAML::Value << r.status;
    out << YAML::Key << "calls" << YAML::Value << YAML::BeginMap;
    out << YAML::Key << "objective" << YAML::Value << r.calls.objective
        << YAML::Key << "objective_grad" << YAML::Value
        << r.calls.objective_grad //
        << YAML::Key << "objective_hess" << YAML::Value
        << r.calls.objective_hess //
        << YAML::Key << "hessian_times_vector" << YAML::Value
        << r.calls.hessian_times_vector;
    if (r.ncon > 0)
        out << YAML::Key << "constraints" << YAML::Value
            << r.calls.constraints //
            << YAML::Key << "constraints_grad" << YAML::Value
            << r.calls.constraints_grad //
            << YAML::Key << "constraints_hess" << YAML::Value
            << r.calls.constraints_hess;
    out << YAML::EndMap;
    out << YAML::Key << "time_setup" << YAML::Value << r.time_setup //
        << YAML::Key << "time" << YAML::Value << r.time;
    out << YAML::EndMap << YAML::EndMap;
    return out;
}

inline YAML::Emitter &operator<<(YAML::Emitter &out, alpaqa::EvalCounter ctr) {
    out << YAML::BeginMap;
    out << YAML::Key << "f" << YAML::Value << ctr.f;
    out << YAML::Key << "grad_f" << YAML::Value << ctr.grad_f;
    out << YAML::Key << "g" << YAML::Value << ctr.g;
    out << YAML::Key << "grad_g_prod" << YAML::Value << ctr.grad_g_prod;
    out << YAML::Key << "grad_gi" << YAML::Value << ctr.grad_gi;
    out << YAML::Key << "hess_L_prod" << YAML::Value << ctr.hess_L_prod;
    out << YAML::Key << "hess_L" << YAML::Value << ctr.hess_L;
    out << YAML::EndMap;
    return out;
}

inline YAML::Emitter &operator<<(YAML::Emitter &out, alpaqa::SolverStatus s) {
    return out << enum_name(s);
}

inline YAML::Emitter &operator<<(YAML::Emitter &out, alpaqa::PANOCStopCrit p) {
    return out << enum_name(p);
}

inline YAML::Emitter &operator<<(YAML::Emitter &out, const alpaqa::LBFGSParams &p) {
    out << YAML::BeginMap;
    out << YAML::Key << "memory" << YAML::Value << p.memory;
    out << YAML::BeginMap;
    out << YAML::Key << "cbfgs" << YAML::Value << YAML::BeginMap;
    out << YAML::Key << "α" << YAML::Value << p.cbfgs.α;
    out << YAML::Key << "ϵ" << YAML::Value << p.cbfgs.ϵ;
    out << YAML::EndMap;
    out << YAML::Key << "rescale_when_γ_changes" << YAML::Value
        << p.rescale_when_γ_changes;
    out << YAML::EndMap;
    return out;
}

inline YAML::Emitter &operator<<(YAML::Emitter &out, const alpaqa::PANOCParams &p) {
    out << YAML::BeginMap;
    out << YAML::Key << "Lipschitz" << YAML::Value << YAML::BeginMap;
    out << YAML::Key << "ε" << YAML::Value << p.Lipschitz.ε;
    out << YAML::Key << "δ" << YAML::Value << p.Lipschitz.δ;
    out << YAML::Key << "Lγ_factor" << YAML::Value << p.Lipschitz.Lγ_factor;
    out << YAML::EndMap;
    out << YAML::Key << "max_iter" << YAML::Value << p.max_iter;
    out << YAML::Key << "max_time" << YAML::Value << p.max_time.count();
    out << YAML::Key << "τ_min" << YAML::Value << p.τ_min;
    out << YAML::Key << "L_min" << YAML::Value << p.L_min;
    out << YAML::Key << "L_max" << YAML::Value << p.L_max;
    out << YAML::Key << "stop_crit" << YAML::Value << p.stop_crit;
    out << YAML::Key << "update_lipschitz_in_linesearch" << YAML::Value
        << p.update_lipschitz_in_linesearch;
    out << YAML::Key << "alternative_linesearch_cond" << YAML::Value
        << p.alternative_linesearch_cond;
    out << YAML::EndMap;
    return out;
}

inline YAML::Emitter &operator<<(YAML::Emitter &out, const alpaqa::ALMParams &p) {
    out << YAML::BeginMap;
    out << YAML::Key << "ε" << YAML::Value << p.ε;
    out << YAML::Key << "δ" << YAML::Value << p.δ;
    out << YAML::Key << "Δ" << YAML::Value << p.Δ;
    out << YAML::Key << "Δ_lower" << YAML::Value << p.Δ_lower;
    out << YAML::Key << "Σ₀" << YAML::Value << p.Σ₀;
    out << YAML::Key << "σ₀" << YAML::Value << p.σ₀;
    out << YAML::Key << "Σ₀_lower" << YAML::Value << p.Σ₀_lower;
    out << YAML::Key << "ε₀" << YAML::Value << p.ε₀;
    out << YAML::Key << "ε₀_increase" << YAML::Value << p.ε₀_increase;
    out << YAML::Key << "ρ" << YAML::Value << p.ρ;
    out << YAML::Key << "ρ_increase" << YAML::Value << p.ρ_increase;
    out << YAML::Key << "θ" << YAML::Value << p.θ;
    out << YAML::Key << "M" << YAML::Value << p.M;
    out << YAML::Key << "Σ_max" << YAML::Value << p.Σ_max;
    out << YAML::Key << "Σ_min" << YAML::Value << p.Σ_min;
    out << YAML::Key << "max_iter" << YAML::Value << p.max_iter;
    out << YAML::Key << "max_time" << YAML::Value << p.max_time.count();
    out << YAML::Key << "max_num_initial_retries" << YAML::Value
        << p.max_num_initial_retries;
    out << YAML::Key << "max_num_retries" << YAML::Value << p.max_num_retries;
    out << YAML::Key << "max_total_num_retries" << YAML::Value
        << p.max_total_num_retries;
    out << YAML::Key << "print_interval" << YAML::Value << p.print_interval;
    out << YAML::Key << "preconditioning" << YAML::Value << p.preconditioning;
    out << YAML::Key << "single_penalty_factor" << YAML::Value
        << p.single_penalty_factor;
    out << YAML::EndMap;
    return out;
}

inline YAML::Emitter &
operator<<(YAML::Emitter &out,
           const alpaqa::InnerStatsAccumulator<alpaqa::PANOCStats> &s) {
    out << YAML::BeginMap;
    out << YAML::Key << "elapsed_time" << YAML::Value << s.elapsed_time.count();
    out << YAML::Key << "iterations" << YAML::Value << s.iterations;
    out << YAML::Key << "linesearch_failures" << YAML::Value
        << s.linesearch_failures;
    out << YAML::Key << "lbfgs_failures" << YAML::Value << s.lbfgs_failures;
    out << YAML::Key << "lbfgs_rejected" << YAML::Value << s.lbfgs_rejected;
    out << YAML::Key << "τ_1_accepted" << YAML::Value << s.τ_1_accepted;
    out << YAML::Key << "sum_τ" << YAML::Value << s.sum_τ;
    out << YAML::Key << "count_τ" << YAML::Value << s.count_τ;
    out << YAML::EndMap;
    return out;
}

inline YAML::Emitter &operator<<(
    YAML::Emitter &out,
    const alpaqa::InnerStatsAccumulator<alpaqa::StructuredPANOCLBFGSSolver::Stats> &s) {
    out << YAML::BeginMap;
    out << YAML::Key << "elapsed_time" << YAML::Value << s.elapsed_time.count();
    out << YAML::Key << "iterations" << YAML::Value << s.iterations;
    out << YAML::Key << "linesearch_failures" << YAML::Value
        << s.linesearch_failures;
    out << YAML::Key << "lbfgs_failures" << YAML::Value << s.lbfgs_failures;
    out << YAML::Key << "lbfgs_rejected" << YAML::Value << s.lbfgs_rejected;
    out << YAML::Key << "τ_1_accepted" << YAML::Value << s.τ_1_accepted;
    out << YAML::Key << "sum_τ" << YAML::Value << s.sum_τ;
    out << YAML::Key << "count_τ" << YAML::Value << s.count_τ;
    out << YAML::EndMap;
    return out;
}

inline YAML::Emitter &
operator<<(YAML::Emitter &out,
           const alpaqa::InnerStatsAccumulator<alpaqa::PGASolver::Stats> &s) {
    out << YAML::BeginMap;
    out << YAML::Key << "elapsed_time" << YAML::Value << s.elapsed_time.count();
    out << YAML::Key << "iterations" << YAML::Value << s.iterations;
    out << YAML::EndMap;
    return out;
}

inline YAML::Emitter &
operator<<(YAML::Emitter &out,
           const alpaqa::InnerStatsAccumulator<alpaqa::GAAPGASolver::Stats> &s) {
    out << YAML::BeginMap;
    out << YAML::Key << "elapsed_time" << YAML::Value << s.elapsed_time.count();
    out << YAML::Key << "iterations" << YAML::Value << s.iterations;
    out << YAML::Key << "accelerated_steps_accepted" << YAML::Value
        << s.accelerated_steps_accepted;
    out << YAML::EndMap;
    return out;
}

inline YAML::Emitter &operator<<(
    YAML::Emitter &out,
    const alpaqa::InnerStatsAccumulator<alpaqa::SecondOrderPANOCSolver::Stats> &s) {
    out << YAML::BeginMap;
    out << YAML::Key << "elapsed_time" << YAML::Value << s.elapsed_time.count();
    out << YAML::Key << "iterations" << YAML::Value << s.iterations;
    out << YAML::Key << "newton_failures" << YAML::Value << s.newton_failures;
    out << YAML::Key << "linesearch_failures" << YAML::Value
        << s.linesearch_failures;
    out << YAML::Key << "τ_1_accepted" << YAML::Value << s.τ_1_accepted;
    out << YAML::Key << "sum_τ" << YAML::Value << s.sum_τ;
    out << YAML::Key << "count_τ" << YAML::Value << s.count_τ;
    out << YAML::EndMap;
    return out;
}

inline YAML::Emitter &
operator<<(YAML::Emitter &out,
           const alpaqa::InnerStatsAccumulator<alpaqa::LBFGSBStats> &s) {
    out << YAML::BeginMap;
    out << YAML::Key << "elapsed_time" << YAML::Value << s.elapsed_time.count();
    out << YAML::Key << "iterations" << YAML::Value << s.iterations;
    out << YAML::EndMap;
    return out;
}