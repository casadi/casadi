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

template <class T>
class attr_setter_fun_t {
  public:
    template <class A>
    attr_setter_fun_t(A T::*attr) : f(attr_setter(attr)) {}
    void operator()(T &t, const py::handle &h) const { f(t, h); }

  private:
    std::function<void(T &, const py::handle &)> f;
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
        it->second(t, val);
    }
}

template <class T>
T kwargs_to_struct(const py::kwargs &kwargs) {
    T t{};
    kwargs_to_struct_helper(t, kwargs);
    return t;
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

#include <panoc-alm/inner/decl/second-order-panoc-lbfgs.hpp>

// TODO: move Lipschitz to its own reusable struct

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