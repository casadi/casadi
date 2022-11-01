/**
 * @file
 * This file defines mappings from Python dicts (kwargs) to simple parameter
 * structs.
 */

#pragma once

#include <alpaqa/inner/panoc.hpp>
// #include <alpaqa/inner/structured-panoc.hpp>

#include <pybind11/chrono.h>
#include <pybind11/pybind11.h>
namespace py = pybind11;

namespace alpaqa {

template <class Inner>
struct ALMSolver;

} // namespace alpaqa

namespace alpaqa::conv {

template <Config Conf>
py::dict stats_to_dict(const PANOCStats<Conf> &s) {
    using namespace py::literals;
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
        "final_γ"_a             = s.final_γ,
    };
}

template <Config Conf>
py::dict stats_to_dict(const InnerStatsAccumulator<PANOCStats<Conf>> &s) {
    using namespace py::literals;
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

// template <Config Conf>
// py::dict stats_to_dict(const StructuredPANOCLBFGSStats<Conf> &s) {
//     using namespace py::literals;
//     return py::dict{
//         "status"_a              = s.status,
//         "ε"_a                   = s.ε,
//         "elapsed_time"_a        = s.elapsed_time,
//         "iterations"_a          = s.iterations,
//         "linesearch_failures"_a = s.linesearch_failures,
//         "lbfgs_failures"_a      = s.lbfgs_failures,
//         "lbfgs_rejected"_a      = s.lbfgs_rejected,
//         "τ_1_accepted"_a        = s.τ_1_accepted,
//         "count_τ"_a             = s.count_τ,
//         "sum_τ"_a               = s.sum_τ,
//         "fpr_shortcuts"_a       = s.fpr_shortcuts,
//         "final_γ"_a             = s.final_γ,
//     };
// }

// template <Config Conf>
// py::dict stats_to_dict(const InnerStatsAccumulator<StructuredPANOCLBFGSStats<Conf>> &s) {
//     using namespace py::literals;
//     return py::dict{
//         "elapsed_time"_a        = s.elapsed_time,
//         "iterations"_a          = s.iterations,
//         "linesearch_failures"_a = s.linesearch_failures,
//         "lbfgs_failures"_a      = s.lbfgs_failures,
//         "lbfgs_rejected"_a      = s.lbfgs_rejected,
//         "τ_1_accepted"_a        = s.τ_1_accepted,
//         "count_τ"_a             = s.count_τ,
//         "sum_τ"_a               = s.sum_τ,
//         "fpr_shortcuts"_a       = s.fpr_shortcuts,
//     };
// }

template <class Inner>
py::dict stats_to_dict(const typename ALMSolver<Inner>::Stats &s) {
    using namespace py::literals;
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
        "inner"_a                      = s.inner.as_dict,
    };
}

} // namespace alpaqa::conv