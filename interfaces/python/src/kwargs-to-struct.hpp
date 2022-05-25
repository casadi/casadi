/**
 * @file
 * This file defines mappings from Python dicts (kwargs) to simple parameter
 * structs.
 */

#pragma once

#include <functional>
#include <map>
#include <variant>

#include <pybind11/detail/typeid.h>
#include <pybind11/pybind11.h>
namespace py = pybind11;

struct cast_error_with_types : py::cast_error {
    cast_error_with_types(const py::cast_error &e, std::string from, std::string to)
        : py::cast_error(e), from(std::move(from)), to(std::move(to)) {}
    std::string from;
    std::string to;
};

template <class T, class A>
auto attr_setter(A T::*attr) {
    return [attr](T &t, const py::handle &h) {
        try {
            t.*attr = h.cast<A>();
        } catch (const py::cast_error &e) {
            throw cast_error_with_types(e, py::str(py::type::handle_of(h)), py::type_id<A>());
        }
    };
}
template <class T, class A>
auto attr_getter(A T::*attr) {
    return [attr](const T &t) { return py::cast(t.*attr); };
}

template <class T>
struct attr_setter_fun_t {
    template <class A>
    attr_setter_fun_t(A T::*attr) : set(attr_setter(attr)), get(attr_getter(attr)) {}

    std::function<void(T &, const py::handle &)> set;
    std::function<py::object(const T &)> get;
};

template <class T>
using kwargs_to_struct_table_t = std::map<std::string, attr_setter_fun_t<T>>;

template <class T>
struct kwargs_to_struct_table;

template <class T>
void kwargs_to_struct_helper(T &t, const py::kwargs &kwargs) {
    const auto &m = kwargs_to_struct_table<T>::table;
    for (auto &&[key, val] : kwargs) {
        auto skey = key.template cast<std::string>();
        auto it   = m.find(skey);
        if (it == m.end())
            throw py::key_error("Unknown parameter " + skey);
        try {
            it->second.set(t, val);
        } catch (const cast_error_with_types &e) {
            throw std::runtime_error("Error converting parameter '" + skey + "' from " + e.from +
                                     " to '" + e.to + "': " + e.what());
        } catch (const std::runtime_error &e) {
            throw std::runtime_error("Error setting parameter '" + skey + "': " + e.what());
        }
    }
}

template <class T>
py::dict struct_to_dict_helper(const T &t) {
    const auto &m = kwargs_to_struct_table<T>::table;
    py::dict d;
    for (auto &&[key, val] : m) {
        py::object o = val.get(t);
        if (py::hasattr(o, "to_dict"))
            o = o.attr("to_dict")();
        d[key.c_str()] = std::move(o);
    }
    return d;
}

template <class T>
T kwargs_to_struct(const py::kwargs &kwargs) {
    T t{};
    kwargs_to_struct_helper(t, kwargs);
    return t;
}

template <class T>
py::dict struct_to_dict(const T &t) {
    return struct_to_dict_helper<T>(t);
}

template <class Params>
using params_or_dict = std::variant<Params, py::dict>;

template <class T>
T var_kwargs_to_struct(const params_or_dict<T> &p) {
    return std::holds_alternative<T>(p) ? std::get<T>(p)
                                        : kwargs_to_struct<T>(std::get<py::dict>(p));
}

#include <alpaqa/inner/panoc.hpp>

template <alpaqa::Config Conf>
struct kwargs_to_struct_table<alpaqa::PANOCParams<Conf>> {
    inline const static kwargs_to_struct_table_t<alpaqa::PANOCParams<Conf>> table{
        {"Lipschitz", &alpaqa::PANOCParams<Conf>::Lipschitz},
        {"max_iter", &alpaqa::PANOCParams<Conf>::max_iter},
        {"max_time", &alpaqa::PANOCParams<Conf>::max_time},
        {"τ_min", &alpaqa::PANOCParams<Conf>::τ_min},
        {"L_min", &alpaqa::PANOCParams<Conf>::L_min},
        {"L_max", &alpaqa::PANOCParams<Conf>::L_max},
        {"stop_crit", &alpaqa::PANOCParams<Conf>::stop_crit},
        {"max_no_progress", &alpaqa::PANOCParams<Conf>::max_no_progress},
        {"print_interval", &alpaqa::PANOCParams<Conf>::print_interval},
        {"quadratic_upperbound_tolerance_factor",
         &alpaqa::PANOCParams<Conf>::quadratic_upperbound_tolerance_factor},
        {"update_lipschitz_in_linesearch",
         &alpaqa::PANOCParams<Conf>::update_lipschitz_in_linesearch},
        {"alternative_linesearch_cond", &alpaqa::PANOCParams<Conf>::alternative_linesearch_cond},
        {"lbfgs_stepsize", &alpaqa::PANOCParams<Conf>::lbfgs_stepsize},
    };
};

template <alpaqa::Config Conf>
struct kwargs_to_struct_table<alpaqa::LipschitzEstimateParams<Conf>> {
    inline const static kwargs_to_struct_table_t<alpaqa::LipschitzEstimateParams<Conf>> table{
        {"L_0", &alpaqa::LipschitzEstimateParams<Conf>::L₀},
        {"δ", &alpaqa::LipschitzEstimateParams<Conf>::δ},
        {"ε", &alpaqa::LipschitzEstimateParams<Conf>::ε},
        {"Lγ_factor", &alpaqa::LipschitzEstimateParams<Conf>::Lγ_factor},
    };
};

#if 0
#include <alpaqa/inner/pga.hpp>

template <>
inline const kwargs_to_struct_table_t<alpaqa::PGAParams>
    kwargs_to_struct_table<alpaqa::PGAParams>{
        {"Lipschitz", &alpaqa::PGAParams::Lipschitz},
        {"max_iter", &alpaqa::PGAParams::max_iter},
        {"max_time", &alpaqa::PGAParams::max_time},
        {"L_min", &alpaqa::PGAParams::L_min},
        {"L_max", &alpaqa::PGAParams::L_max},
        {"stop_crit", &alpaqa::PGAParams::stop_crit},
        {"print_interval", &alpaqa::PGAParams::print_interval},
        {"quadratic_upperbound_tolerance_factor",
         &alpaqa::PGAParams::quadratic_upperbound_tolerance_factor},
    };

#include <alpaqa/inner/guarded-aa-pga.hpp>

template <>
inline const kwargs_to_struct_table_t<alpaqa::GAAPGAParams>
    kwargs_to_struct_table<alpaqa::GAAPGAParams>{
        {"Lipschitz", &alpaqa::GAAPGAParams::Lipschitz},
        {"limitedqr_mem", &alpaqa::GAAPGAParams::limitedqr_mem},
        {"max_iter", &alpaqa::GAAPGAParams::max_iter},
        {"max_time", &alpaqa::GAAPGAParams::max_time},
        {"L_min", &alpaqa::GAAPGAParams::L_min},
        {"L_max", &alpaqa::GAAPGAParams::L_max},
        {"stop_crit", &alpaqa::GAAPGAParams::stop_crit},
        {"print_interval", &alpaqa::GAAPGAParams::print_interval},
        {"quadratic_upperbound_tolerance_factor",
         &alpaqa::GAAPGAParams::quadratic_upperbound_tolerance_factor},
        {"max_no_progress", &alpaqa::GAAPGAParams::max_no_progress},
        {"full_flush_on_γ_change",
         &alpaqa::GAAPGAParams::full_flush_on_γ_change},
    };

#include <alpaqa/inner/decl/structured-panoc-lbfgs.hpp>

template <>
inline const kwargs_to_struct_table_t<alpaqa::StructuredPANOCLBFGSParams>
    kwargs_to_struct_table<alpaqa::StructuredPANOCLBFGSParams>{
        {"Lipschitz", &alpaqa::StructuredPANOCLBFGSParams::Lipschitz},
        {"max_iter", &alpaqa::StructuredPANOCLBFGSParams::max_iter},
        {"max_time", &alpaqa::StructuredPANOCLBFGSParams::max_time},
        {"τ_min", &alpaqa::StructuredPANOCLBFGSParams::τ_min},
        {"L_min", &alpaqa::StructuredPANOCLBFGSParams::L_min},
        {"L_max", &alpaqa::StructuredPANOCLBFGSParams::L_max},
        {"nonmonotone_linesearch",
         &alpaqa::StructuredPANOCLBFGSParams::nonmonotone_linesearch},
        {"fpr_shortcut_accept_factor",
         &alpaqa::StructuredPANOCLBFGSParams::fpr_shortcut_accept_factor},
        {"fpr_shortcut_history",
         &alpaqa::StructuredPANOCLBFGSParams::fpr_shortcut_history},
        {"stop_crit", &alpaqa::StructuredPANOCLBFGSParams::stop_crit},
        {"max_no_progress",
         &alpaqa::StructuredPANOCLBFGSParams::max_no_progress},
        {"print_interval", &alpaqa::StructuredPANOCLBFGSParams::print_interval},
        {"quadratic_upperbound_tolerance_factor",
         &alpaqa::StructuredPANOCLBFGSParams::
             quadratic_upperbound_tolerance_factor},
        {"update_lipschitz_in_linesearch",
         &alpaqa::StructuredPANOCLBFGSParams::update_lipschitz_in_linesearch},
        {"alternative_linesearch_cond",
         &alpaqa::StructuredPANOCLBFGSParams::alternative_linesearch_cond},
        {"hessian_vec", &alpaqa::StructuredPANOCLBFGSParams::hessian_vec},
        {"hessian_vec_finite_differences",
         &alpaqa::StructuredPANOCLBFGSParams::hessian_vec_finite_differences},
        {"full_augmented_hessian",
         &alpaqa::StructuredPANOCLBFGSParams::full_augmented_hessian},
        {"hessian_step_size_heuristic",
         &alpaqa::StructuredPANOCLBFGSParams::hessian_step_size_heuristic},
        {"lbfgs_stepsize", &alpaqa::StructuredPANOCLBFGSParams::lbfgs_stepsize},
    };
#endif

#include <alpaqa/inner/directions/panoc/lbfgs.hpp>

template <alpaqa::Config Conf>
struct kwargs_to_struct_table<alpaqa::LBFGSParams<Conf>> {
    inline const static kwargs_to_struct_table_t<alpaqa::LBFGSParams<Conf>> table{
        {"memory", &alpaqa::LBFGSParams<Conf>::memory},
        {"cbfgs", &alpaqa::LBFGSParams<Conf>::cbfgs},
    };
};

template <alpaqa::Config Conf>
struct kwargs_to_struct_table<alpaqa::CBFGSParams<Conf>> {
    inline const static kwargs_to_struct_table_t<alpaqa::CBFGSParams<Conf>> table{
        {"α", &alpaqa::CBFGSParams<Conf>::α},
        {"ϵ", &alpaqa::CBFGSParams<Conf>::ϵ},
    };
};

#include <alpaqa/outer/alm.hpp>

template <alpaqa::Config Conf>
struct kwargs_to_struct_table<alpaqa::ALMParams<Conf>> {
    inline const static kwargs_to_struct_table_t<alpaqa::ALMParams<Conf>> table{
        {"ε", &alpaqa::ALMParams<Conf>::ε},
        {"δ", &alpaqa::ALMParams<Conf>::δ},
        {"Δ", &alpaqa::ALMParams<Conf>::Δ},
        {"Δ_lower", &alpaqa::ALMParams<Conf>::Δ_lower},
        {"Σ_0", &alpaqa::ALMParams<Conf>::Σ₀},
        {"σ_0", &alpaqa::ALMParams<Conf>::σ₀},
        {"Σ_0_lower", &alpaqa::ALMParams<Conf>::Σ₀_lower},
        {"ε_0", &alpaqa::ALMParams<Conf>::ε₀},
        {"ε_0_increase", &alpaqa::ALMParams<Conf>::ε₀_increase},
        {"ρ", &alpaqa::ALMParams<Conf>::ρ},
        {"ρ_increase", &alpaqa::ALMParams<Conf>::ρ_increase},
        {"θ", &alpaqa::ALMParams<Conf>::θ},
        {"M", &alpaqa::ALMParams<Conf>::M},
        {"Σ_max", &alpaqa::ALMParams<Conf>::Σ_max},
        {"Σ_min", &alpaqa::ALMParams<Conf>::Σ_min},
        {"max_iter", &alpaqa::ALMParams<Conf>::max_iter},
        {"max_time", &alpaqa::ALMParams<Conf>::max_time},
        {"max_num_initial_retries", &alpaqa::ALMParams<Conf>::max_num_initial_retries},
        {"max_num_retries", &alpaqa::ALMParams<Conf>::max_num_retries},
        {"max_total_num_retries", &alpaqa::ALMParams<Conf>::max_total_num_retries},
        {"print_interval", &alpaqa::ALMParams<Conf>::print_interval},
        {"single_penalty_factor", &alpaqa::ALMParams<Conf>::single_penalty_factor},
    };
};
