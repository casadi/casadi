#pragma once

#include <alpaqa/inner/decl/structured-panoc-lbfgs.hpp>
#include <alpaqa/inner/directions/lbfgs.hpp>
#include <alpaqa/inner/guarded-aa-pga.hpp>
#include <alpaqa/inner/panoc.hpp>
#include <alpaqa/inner/pga.hpp>
#include <alpaqa/util/solverstatus.hpp>

#include <memory>
#include <type_traits>

#include <pybind11/cast.h>
#include <pybind11/pybind11.h>
namespace py = pybind11;

namespace alpaqa {

template <class InnerSolver>
auto InnerSolverCallWrapper() {
    return [](InnerSolver &solver, const alpaqa::Problem &p, alpaqa::crvec Σ,
              alpaqa::real_t ε, alpaqa::vec x,
              alpaqa::vec y) -> std::tuple<alpaqa::vec, alpaqa::vec, alpaqa::vec, py::dict> {
        alpaqa::vec z(p.m);
        auto stats = solver(p, Σ, ε, true, x, y, z);
        return std::make_tuple(std::move(x), std::move(y), std::move(z),
                               stats.ptr->to_dict());
    };
}

class PolymorphicInnerSolverStatsAccumulatorBase
    : public std::enable_shared_from_this<
          PolymorphicInnerSolverStatsAccumulatorBase> {
  public:
    virtual ~PolymorphicInnerSolverStatsAccumulatorBase() = default;
    virtual py::dict to_dict() const                      = 0;
    virtual void accumulate(const class PolymorphicInnerSolverStatsBase &) = 0;
};

class PolymorphicInnerSolverStatsBase
    : public std::enable_shared_from_this<PolymorphicInnerSolverStatsBase> {
  public:
    virtual ~PolymorphicInnerSolverStatsBase() = default;
    virtual py::dict to_dict() const           = 0;
    virtual std::shared_ptr<PolymorphicInnerSolverStatsAccumulatorBase>
    accumulator() const = 0;
};

class PolymorphicInnerSolverBase
    : public std::enable_shared_from_this<PolymorphicInnerSolverBase> {
  public:
    struct Stats {
        std::shared_ptr<PolymorphicInnerSolverStatsBase> ptr;
        SolverStatus status;
        real_t ε;
        unsigned iterations;

        static Stats from_dict(py::dict d) {
            using PolyStats    = alpaqa::PolymorphicInnerSolverStatsBase;
            using PolyAccStats = alpaqa::PolymorphicInnerSolverStatsAccumulatorBase;
            using InnerStats   = alpaqa::PolymorphicInnerSolverBase::Stats;
            struct AccStats : PolyAccStats {
                AccStats(py::dict dict) : dict(std::move(dict)) {}
                py::dict dict;
                py::dict to_dict() const override { return dict; }
                void accumulate(const PolyStats &s) override {
                    if (this->dict.contains("accumulate"))
                        this->dict["accumulate"](this->dict, s.to_dict());
                    else
                        throw py::key_error("Stats accumulator does not define "
                                            "an accumulate function");
                }
            };
            struct Stats : PolyStats {
                Stats(py::dict dict) : dict(std::move(dict)) {}
                py::dict dict;
                py::dict to_dict() const override { return dict; }
                std::shared_ptr<PolyAccStats> accumulator() const override {
                    if (this->dict.contains("accumulator"))
                        return {
                            std::make_shared<AccStats>(
                                dict["accumulator"].cast<py::dict>()),
                        };
                    else
                        throw py::key_error(
                            "Stats do not define an accumulator");
                }
            };
            bool ok = d.contains("status") && d.contains("ε") &&
                      d.contains("iterations");
            if (not ok)
                throw py::key_error(
                    "Stats should contain status, ε and iterations");
            return {
                std::static_pointer_cast<PolyStats>(std::make_shared<Stats>(d)),
                d["status"].cast<decltype(InnerStats::status)>(),
                d["ε"].cast<decltype(InnerStats::ε)>(),
                d["iterations"].cast<decltype(InnerStats::iterations)>(),
            };
        }
    };

    virtual ~PolymorphicInnerSolverBase() = default;
    virtual Stats operator()(
        /// [in]    Problem description
        const Problem &problem,
        /// [in]    Constraint weights @f$ \Sigma @f$
        crvec Σ,
        /// [in]    Tolerance @f$ \varepsilon @f$
        real_t ε,
        /// [in]    Overwrite @p x, @p y and @p err_z even if not converged
        bool always_overwrite_results,
        /// [inout] Decision variable @f$ x @f$
        rvec x,
        /// [inout] Lagrange multipliers @f$ y @f$
        rvec y,
        /// [out]   Slack variable error @f$ g(x) - z @f$
        rvec err_z)                       = 0;
    virtual void stop()                   = 0;
    virtual std::string get_name() const  = 0;
    virtual py::object get_params() const = 0;
};

struct PolymorphicInnerSolverWrapper {
    using Stats = PolymorphicInnerSolverBase::Stats;
    std::shared_ptr<PolymorphicInnerSolverBase> solver;
    PolymorphicInnerSolverWrapper(
        std::shared_ptr<PolymorphicInnerSolverBase> &&solver)
        : solver(std::move(solver)) {}

    Stats operator()(const Problem &problem, crvec Σ, real_t ε,
                     bool always_overwrite_results, rvec x, rvec y,
                     rvec err_z) {
        return solver->operator()(problem, Σ, ε, always_overwrite_results, x, y,
                                  err_z);
    }
    void stop() { solver->stop(); }
    std::string get_name() const { return solver->get_name(); }
    py::object get_params() const { return solver->get_params(); }
};

template <class InnerSolverStats>
struct InnerStatsAccumulator;

template <>
struct InnerStatsAccumulator<PolymorphicInnerSolverWrapper::Stats> {
    std::shared_ptr<PolymorphicInnerSolverStatsAccumulatorBase> ptr;
    py::dict to_dict() const { return ptr->to_dict(); }
};

inline InnerStatsAccumulator<PolymorphicInnerSolverWrapper::Stats> &
operator+=(InnerStatsAccumulator<PolymorphicInnerSolverWrapper::Stats> &acc,
           const PolymorphicInnerSolverWrapper::Stats &s) {
    assert(s.ptr);
    if (not acc.ptr)
        acc.ptr = s.ptr->accumulator();
    acc.ptr->accumulate(*s.ptr);
    return acc;
}

class PolymorphicInnerSolverTrampoline : public PolymorphicInnerSolverBase {
  public:
    Stats operator()(const Problem &problem, crvec Σ, real_t ε,
                     bool always_overwrite_results, rvec x, rvec y,
                     rvec err_z) override {
        py::dict stats;
        std::tie(x, y, err_z, stats) =
            call(problem, Σ, ε, always_overwrite_results, x, y);
        return Stats::from_dict(stats);
    }
    virtual std::tuple<alpaqa::vec, alpaqa::vec, alpaqa::vec, py::dict>
    call(const alpaqa::Problem &problem, alpaqa::crvec Σ, alpaqa::real_t ε,
         bool always_overwrite_results, alpaqa::vec x, alpaqa::vec y) {
        using ret = std::tuple<alpaqa::vec, alpaqa::vec, alpaqa::vec, py::dict>;
        PYBIND11_OVERRIDE_PURE_NAME(ret, PolymorphicInnerSolverBase, "__call__",
                                    call, problem, Σ, ε,
                                    always_overwrite_results, x, y);
    }
    std::string get_name() const override {
        PYBIND11_OVERRIDE_PURE(std::string, PolymorphicInnerSolverBase,
                               get_name, );
    }
    py::object get_params() const override {
        PYBIND11_OVERRIDE_PURE(py::object, PolymorphicInnerSolverBase,
                               get_params, );
    }
    void stop() override {
        PYBIND11_OVERRIDE_PURE(void, PolymorphicInnerSolverBase, stop, );
    }
};

inline py::dict stats_to_dict(const PANOCStats &s) {
    using py::operator""_a;
    return py::dict{
        "status"_a              = s.status,
        "ε"_a                   = s.ε,
        "elapsed_time"_a        = s.elapsed_time,
        "iterations"_a          = s.iterations,
        "linesearch_failures"_a = s.linesearch_failures,
        "lbfgs_failures"_a      = s.lbfgs_failures,
        "lbfgs_rejected"_a      = s.lbfgs_rejected,
        "τ_1_accepted"_a        = s.τ_1_accepted,
        "count_τ"_a             = s.count_τ,
        "sum_τ"_a               = s.sum_τ,
    };
}

inline py::dict stats_to_dict(const InnerStatsAccumulator<PANOCStats> &s) {
    using py::operator""_a;
    return py::dict{
        "elapsed_time"_a        = s.elapsed_time,
        "iterations"_a          = s.iterations,
        "linesearch_failures"_a = s.linesearch_failures,
        "lbfgs_failures"_a      = s.lbfgs_failures,
        "lbfgs_rejected"_a      = s.lbfgs_rejected,
        "τ_1_accepted"_a        = s.τ_1_accepted,
        "count_τ"_a             = s.count_τ,
        "sum_τ"_a               = s.sum_τ,
    };
}

inline py::dict stats_to_dict(const StructuredPANOCLBFGSSolver::Stats &s) {
    using py::operator""_a;
    return py::dict{
        "status"_a              = s.status,
        "ε"_a                   = s.ε,
        "elapsed_time"_a        = s.elapsed_time,
        "iterations"_a          = s.iterations,
        "linesearch_failures"_a = s.linesearch_failures,
        "lbfgs_failures"_a      = s.lbfgs_failures,
        "lbfgs_rejected"_a      = s.lbfgs_rejected,
        "τ_1_accepted"_a        = s.τ_1_accepted,
        "count_τ"_a             = s.count_τ,
        "sum_τ"_a               = s.sum_τ,
    };
}

inline py::dict stats_to_dict(const PGASolver::Stats &s) {
    using py::operator""_a;
    return py::dict{
        "status"_a       = s.status,
        "ε"_a            = s.ε,
        "elapsed_time"_a = s.elapsed_time,
        "iterations"_a   = s.iterations,
    };
}

inline py::dict stats_to_dict(const GAAPGASolver::Stats &s) {
    using py::operator""_a;
    return py::dict{
        "status"_a                     = s.status,
        "ε"_a                          = s.ε,
        "elapsed_time"_a               = s.elapsed_time,
        "iterations"_a                 = s.iterations,
        "accelerated_steps_accepted"_a = s.accelerated_steps_accepted,
    };
}

inline py::dict stats_to_dict(
    const InnerStatsAccumulator<StructuredPANOCLBFGSSolver::Stats> &s) {
    using py::operator""_a;
    return py::dict{
        "elapsed_time"_a        = s.elapsed_time,
        "iterations"_a          = s.iterations,
        "linesearch_failures"_a = s.linesearch_failures,
        "lbfgs_failures"_a      = s.lbfgs_failures,
        "lbfgs_rejected"_a      = s.lbfgs_rejected,
        "τ_1_accepted"_a        = s.τ_1_accepted,
        "count_τ"_a             = s.count_τ,
        "sum_τ"_a               = s.sum_τ,
    };
}

inline py::dict
stats_to_dict(const InnerStatsAccumulator<PGASolver::Stats> &s) {
    using py::operator""_a;
    return py::dict{
        "elapsed_time"_a = s.elapsed_time,
        "iterations"_a   = s.iterations,
    };
}

inline py::dict
stats_to_dict(const InnerStatsAccumulator<GAAPGASolver::Stats> &s) {
    using py::operator""_a;
    return py::dict{
        "elapsed_time"_a               = s.elapsed_time,
        "iterations"_a                 = s.iterations,
        "accelerated_steps_accepted"_a = s.accelerated_steps_accepted,
    };
}

template <class InnerSolver>
class PolymorphicInnerSolver : public PolymorphicInnerSolverBase {
  public:
    PolymorphicInnerSolver(InnerSolver &&innersolver)
        : innersolver(std::forward<InnerSolver>(innersolver)) {}
    PolymorphicInnerSolver(const InnerSolver &innersolver)
        : innersolver(innersolver) {}
    template <class... Args>
    PolymorphicInnerSolver(Args... args)
        : innersolver(InnerSolver{std::forward<Args>(args)...}) {}

    struct WrappedStatsAccumulator
        : PolymorphicInnerSolverStatsAccumulatorBase {
        InnerStatsAccumulator<typename InnerSolver::Stats> acc;
        void
        accumulate(const PolymorphicInnerSolverStatsBase &bstats) override {
            auto &stats = dynamic_cast<const WrappedStats &>(bstats).stats;
            acc += stats;
        }
        py::dict to_dict() const override { return stats_to_dict(acc); }
    };
    struct WrappedStats : PolymorphicInnerSolverStatsBase {
        using Stats = typename InnerSolver::Stats;
        WrappedStats(const Stats &stats) : stats(stats) {}
        Stats stats;
        std::shared_ptr<PolymorphicInnerSolverStatsAccumulatorBase>
        accumulator() const override {
            return std::static_pointer_cast<
                PolymorphicInnerSolverStatsAccumulatorBase>(
                std::make_shared<WrappedStatsAccumulator>());
        }
        py::dict to_dict() const override { return stats_to_dict(stats); }
    };

    Stats operator()(
        /// [in]    Problem description
        const Problem &problem,
        /// [in]    Constraint weights @f$ \Sigma @f$
        crvec Σ,
        /// [in]    Tolerance @f$ \varepsilon @f$
        real_t ε,
        /// [in]    Overwrite @p x, @p y and @p err_z even if not converged
        bool always_overwrite_results,
        /// [inout] Decision variable @f$ x @f$
        rvec x,
        /// [inout] Lagrange multipliers @f$ y @f$
        rvec y,
        /// [out]   Slack variable error @f$ g(x) - z @f$
        rvec err_z) override {
        auto stats =
            innersolver(problem, Σ, ε, always_overwrite_results, x, y, err_z);
        return {
            std::static_pointer_cast<PolymorphicInnerSolverStatsBase>(
                std::make_shared<WrappedStats>(stats)),
            stats.status,
            stats.ε,
            stats.iterations,
        };
    }
    void stop() override { innersolver.stop(); }
    std::string get_name() const override { return innersolver.get_name(); }
    py::object get_params() const override {
        return py::cast(innersolver.get_params());
    }

    void set_progress_callback(
        std::function<void(const typename InnerSolver::ProgressInfo &)> cb) {
        this->innersolver.set_progress_callback(std::move(cb));
    }

    InnerSolver innersolver;
};

} // namespace alpaqa

#include "polymorphic-panoc-direction.hpp"
#include <alpaqa/alm.hpp>

namespace alpaqa {

using PolymorphicPGASolver    = PolymorphicInnerSolver<PGASolver>;
using PolymorphicGAAPGASolver = PolymorphicInnerSolver<GAAPGASolver>;
using PolymorphicPANOCSolver =
    PolymorphicInnerSolver<PANOCSolver<PolymorphicPANOCDirectionBase>>;
using PolymorphicStructuredPANOCLBFGSSolver =
    PolymorphicInnerSolver<StructuredPANOCLBFGSSolver>;

using PolymorphicALMSolver = ALMSolver<PolymorphicInnerSolverWrapper>;

inline py::dict stats_to_dict(const PolymorphicALMSolver::Stats &s) {
    using py::operator""_a;
    return py::dict{
        "outer_iterations"_a           = s.outer_iterations,
        "elapsed_time"_a               = s.elapsed_time,
        "initial_penalty_reduced"_a    = s.initial_penalty_reduced,
        "penalty_reduced"_a            = s.penalty_reduced,
        "inner_convergence_failures"_a = s.inner_convergence_failures,
        "ε"_a                          = s.ε,
        "δ"_a                          = s.δ,
        "norm_penalty"_a               = s.norm_penalty,
        "status"_a                     = s.status,
        "inner"_a                      = s.inner.to_dict(),
    };
}

} // namespace alpaqa