/**
 * @file
 * This file defines mappings from Python dicts (kwargs) to simple parameter
 * structs.
 */

#pragma once

#include <functional>
#include <map>

#include <pybind11/pybind11.h>
namespace py = pybind11;

template <class T, class A>
auto attr_setter(A T::*attr) {
    return [attr](T &t, const py::handle &h) { t.*attr = h.cast<A>(); };
}
template <class T, class A>
auto attr_getter(A T::*attr) {
    return [attr](const T &t) { return py::cast(t.*attr); };
}

template <class T>
class attr_setter_fun_t {
  public:
    template <class A>
    attr_setter_fun_t(A T::*attr)
        : set(attr_setter(attr)), get(attr_getter(attr)) {}

    std::function<void(T &, const py::handle &)> set;
    std::function<py::object(const T &)> get;
};

template <class T>
using kwargs_to_struct_table_t = std::map<std::string, attr_setter_fun_t<T>>;

template <class T>
kwargs_to_struct_table_t<T> kwargs_to_struct_table;

template <class T>
void kwargs_to_struct_helper(T &t, const py::kwargs &kwargs) {
    const auto &m = kwargs_to_struct_table<T>;
    for (auto &&[key, val] : kwargs) {
        auto skey = key.template cast<std::string>();
        auto it   = m.find(skey);
        if (it == m.end())
            throw py::key_error("Unknown parameter " + skey);
        it->second.set(t, val);
    }
}

template <class T>
py::dict struct_to_dict_helper(const T &t) {
    const auto &m = kwargs_to_struct_table<T>;
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

#include <panoc-alm/inner/decl/panoc.hpp>

template <>
inline const kwargs_to_struct_table_t<pa::PANOCParams>
    kwargs_to_struct_table<pa::PANOCParams>{
        {"Lipschitz", &pa::PANOCParams::Lipschitz},
        {"max_iter", &pa::PANOCParams::max_iter},
        {"max_time", &pa::PANOCParams::max_time},
        {"τ_min", &pa::PANOCParams::τ_min},
        {"γ_min", &pa::PANOCParams::γ_min},
        {"stop_crit", &pa::PANOCParams::stop_crit},
        {"max_no_progress", &pa::PANOCParams::max_no_progress},
        {"print_interval", &pa::PANOCParams::print_interval},
        {"quadratic_upperbound_tolerance_factor",
         &pa::PANOCParams::quadratic_upperbound_tolerance_factor},
        {"update_lipschitz_in_linesearch",
         &pa::PANOCParams::update_lipschitz_in_linesearch},
        {"alternative_linesearch_cond",
         &pa::PANOCParams::alternative_linesearch_cond},
        {"lbfgs_stepsize", &pa::PANOCParams::lbfgs_stepsize},
    };

template <>
inline const kwargs_to_struct_table_t<decltype(pa::PANOCParams::Lipschitz)>
    kwargs_to_struct_table<decltype(pa::PANOCParams::Lipschitz)>{
        {"L_0", &decltype(pa::PANOCParams::Lipschitz)::L₀},
        {"δ", &decltype(pa::PANOCParams::Lipschitz)::δ},
        {"ε", &decltype(pa::PANOCParams::Lipschitz)::ε},
        {"Lγ_factor", &decltype(pa::PANOCParams::Lipschitz)::Lγ_factor},
    };

#include <panoc-alm/inner/pga.hpp>

template <>
inline const kwargs_to_struct_table_t<decltype(pa::PGAParams::Lipschitz)>
    kwargs_to_struct_table<decltype(pa::PGAParams::Lipschitz)>{
        {"L_0", &decltype(pa::PGAParams::Lipschitz)::L₀},
        {"δ", &decltype(pa::PGAParams::Lipschitz)::δ},
        {"ε", &decltype(pa::PGAParams::Lipschitz)::ε},
        {"Lγ_factor", &decltype(pa::PGAParams::Lipschitz)::Lγ_factor},
    };

template <>
inline const kwargs_to_struct_table_t<pa::PGAParams>
    kwargs_to_struct_table<pa::PGAParams>{
        {"Lipschitz", &pa::PGAParams::Lipschitz},
        {"max_iter", &pa::PGAParams::max_iter},
        {"max_time", &pa::PGAParams::max_time},
        {"γ_min", &pa::PGAParams::γ_min},
        {"stop_crit", &pa::PGAParams::stop_crit},
        {"print_interval", &pa::PGAParams::print_interval},
        {"quadratic_upperbound_tolerance_factor",
         &pa::PGAParams::quadratic_upperbound_tolerance_factor},
    };

#include <panoc-alm/inner/guarded-aa-pga.hpp>

template <>
inline const kwargs_to_struct_table_t<pa::GuardedAAPGAParams>
    kwargs_to_struct_table<pa::GuardedAAPGAParams>{
        {"Lipschitz", &pa::GuardedAAPGAParams::Lipschitz},
        {"limitedqr_mem", &pa::GuardedAAPGAParams::limitedqr_mem},
        {"max_iter", &pa::GuardedAAPGAParams::max_iter},
        {"max_time", &pa::GuardedAAPGAParams::max_time},
        {"γ_min", &pa::GuardedAAPGAParams::γ_min},
        {"stop_crit", &pa::GuardedAAPGAParams::stop_crit},
        {"print_interval", &pa::GuardedAAPGAParams::print_interval},
        {"quadratic_upperbound_tolerance_factor",
         &pa::GuardedAAPGAParams::quadratic_upperbound_tolerance_factor},
        {"max_no_progress", &pa::GuardedAAPGAParams::max_no_progress},
        {"full_flush_on_γ_change",
         &pa::GuardedAAPGAParams::full_flush_on_γ_change},
    };

#include <panoc-alm/inner/decl/second-order-panoc-lbfgs.hpp>

// TODO: move Lipschitz to its own reusable struct

template <>
inline const kwargs_to_struct_table_t<decltype(
    pa::SecondOrderPANOCLBFGSParams::Lipschitz)>
    kwargs_to_struct_table<decltype(
        pa::SecondOrderPANOCLBFGSParams::Lipschitz)>{
        {"L_0", &decltype(pa::SecondOrderPANOCLBFGSParams::Lipschitz)::L₀},
        {"δ", &decltype(pa::SecondOrderPANOCLBFGSParams::Lipschitz)::δ},
        {"ε", &decltype(pa::SecondOrderPANOCLBFGSParams::Lipschitz)::ε},
        {"Lγ_factor",
         &decltype(pa::SecondOrderPANOCLBFGSParams::Lipschitz)::Lγ_factor},
    };

template <>
inline const kwargs_to_struct_table_t<pa::SecondOrderPANOCLBFGSParams>
    kwargs_to_struct_table<pa::SecondOrderPANOCLBFGSParams>{
        {"Lipschitz", &pa::SecondOrderPANOCLBFGSParams::Lipschitz},
        {"max_iter", &pa::SecondOrderPANOCLBFGSParams::max_iter},
        {"max_time", &pa::SecondOrderPANOCLBFGSParams::max_time},
        {"τ_min", &pa::SecondOrderPANOCLBFGSParams::τ_min},
        {"γ_min", &pa::SecondOrderPANOCLBFGSParams::γ_min},
        {"nonmonotone_linesearch",
         &pa::SecondOrderPANOCLBFGSParams::nonmonotone_linesearch},
        {"stop_crit", &pa::SecondOrderPANOCLBFGSParams::stop_crit},
        {"max_no_progress", &pa::SecondOrderPANOCLBFGSParams::max_no_progress},
        {"print_interval", &pa::SecondOrderPANOCLBFGSParams::print_interval},
        {"quadratic_upperbound_tolerance_factor",
         &pa::SecondOrderPANOCLBFGSParams::
             quadratic_upperbound_tolerance_factor},
        {"update_lipschitz_in_linesearch",
         &pa::SecondOrderPANOCLBFGSParams::update_lipschitz_in_linesearch},
        {"alternative_linesearch_cond",
         &pa::SecondOrderPANOCLBFGSParams::alternative_linesearch_cond},
        {"hessian_vec_finited_differences",
         &pa::SecondOrderPANOCLBFGSParams::hessian_vec_finited_differences},
        {"full_augmented_hessian",
         &pa::SecondOrderPANOCLBFGSParams::full_augmented_hessian},
        {"lbfgs_stepsize", &pa::SecondOrderPANOCLBFGSParams::lbfgs_stepsize},
    };

#include <panoc-alm/inner/directions/decl/lbfgs.hpp>

template <>
inline const kwargs_to_struct_table_t<pa::LBFGSParams>
    kwargs_to_struct_table<pa::LBFGSParams>{
        {"memory", &pa::LBFGSParams::memory},
        {"cbfgs", &pa::LBFGSParams::cbfgs},
        {"rescale_when_γ_changes", &pa::LBFGSParams::rescale_when_γ_changes},
    };

template <>
inline const kwargs_to_struct_table_t<decltype(pa::LBFGSParams::cbfgs)>
    kwargs_to_struct_table<decltype(pa::LBFGSParams::cbfgs)>{
        {"α", &decltype(pa::LBFGSParams::cbfgs)::α},
        {"ϵ", &decltype(pa::LBFGSParams::cbfgs)::ϵ},
    };

#include <panoc-alm/decl/alm.hpp>

template <>
inline const kwargs_to_struct_table_t<pa::ALMParams>
    kwargs_to_struct_table<pa::ALMParams>{
        {"ε", &pa::ALMParams::ε},
        {"δ", &pa::ALMParams::δ},
        {"Δ", &pa::ALMParams::Δ},
        {"Δ_lower", &pa::ALMParams::Δ_lower},
        {"Σ_0", &pa::ALMParams::Σ₀},
        {"σ_0", &pa::ALMParams::σ₀},
        {"Σ_0_lower", &pa::ALMParams::Σ₀_lower},
        {"ε_0", &pa::ALMParams::ε₀},
        {"ε_0_increase", &pa::ALMParams::ε₀_increase},
        {"ρ", &pa::ALMParams::ρ},
        {"ρ_increase", &pa::ALMParams::ρ_increase},
        {"θ", &pa::ALMParams::θ},
        {"M", &pa::ALMParams::M},
        {"Σ_max", &pa::ALMParams::Σ_max},
        {"Σ_min", &pa::ALMParams::Σ_min},
        {"max_iter", &pa::ALMParams::max_iter},
        {"max_time", &pa::ALMParams::max_time},
        {"max_num_initial_retries", &pa::ALMParams::max_num_initial_retries},
        {"max_num_retries", &pa::ALMParams::max_num_retries},
        {"max_total_num_retries", &pa::ALMParams::max_total_num_retries},
        {"print_interval", &pa::ALMParams::print_interval},
        {"preconditioning", &pa::ALMParams::preconditioning},
        {"single_penalty_factor", &pa::ALMParams::single_penalty_factor},
    };