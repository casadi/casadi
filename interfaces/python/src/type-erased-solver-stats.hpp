#pragma once

#include <alpaqa/config/config.hpp>
#include <alpaqa/util/type-erasure.hpp>

#include <any>
#include <new>
#include <stdexcept>
#include <type_traits>

#include "stats-to-dict.hpp"

#include <pybind11/pybind11.h>
namespace py = pybind11;

namespace alpaqa {

template <class InnerSolverStats>
struct InnerStatsAccumulator;

template <Config Conf>
struct TypeErasedInnerSolverStats;

template <Config Conf>
struct InnerStatsAccumulator<TypeErasedInnerSolverStats<Conf>> {
    std::any accumulator;
    py::dict as_dict;
};

template <Config Conf, class Stats>
InnerStatsAccumulator<TypeErasedInnerSolverStats<Conf>> &
operator+=(InnerStatsAccumulator<TypeErasedInnerSolverStats<Conf>> &acc, const Stats &stats) {
    using ActualAccumulator = InnerStatsAccumulator<Stats>;
    if (!acc.accumulator.has_value())
        acc.accumulator = ActualAccumulator{};
    auto *act_acc = std::any_cast<ActualAccumulator>(&acc.accumulator);
    if (!act_acc)
        throw std::logic_error("Cannot combine different types of solver stats");
    *act_acc += stats;
    acc.as_dict = conv::stats_to_dict(*act_acc);
    return acc;
}

template <Config Conf>
struct TypeErasedInnerSolverStats {
    USING_ALPAQA_CONFIG(Conf);
    using Accumulator = InnerStatsAccumulator<TypeErasedInnerSolverStats<Conf>>;
    void (*combine_p)(Accumulator &acc, const std::any &stats) = nullptr;
    py::dict (*to_dict_p)(const std::any &self)                = nullptr;
    SolverStatus status;
    real_t ε;
    unsigned iterations;

    std::any stats;

    template <class StatsR>
    TypeErasedInnerSolverStats(StatsR &&stats)
        : status(stats.status), ε(stats.ε), iterations(stats.iterations),
          stats(std::forward<StatsR>(stats)) {
        using Stats = std::remove_cvref_t<StatsR>;
        combine_p   = [](Accumulator &acc, const std::any &stats) {
            auto *act_stats = std::any_cast<Stats>(&stats);
            assert(act_stats);
            acc += *act_stats;
        };
        to_dict_p = [](const std::any &self) {
            auto *act_self = std::any_cast<Stats>(&self);
            assert(act_self);
            return conv::stats_to_dict(*act_self);
        };
    }

    void combine(Accumulator &acc) const { return combine_p(acc, stats); }
    py::dict to_dict() const { return to_dict_p(stats); }
};

template <Config Conf>
InnerStatsAccumulator<TypeErasedInnerSolverStats<Conf>> &
operator+=(InnerStatsAccumulator<TypeErasedInnerSolverStats<Conf>> &acc,
           const TypeErasedInnerSolverStats<Conf> &stats) {
    stats.combine(acc);
    return acc;
}

} // namespace alpaqa