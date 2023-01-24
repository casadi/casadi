/**
 * @file
 * This file defines mappings from Python dicts (kwargs) to simple parameter
 * structs.
 */

#pragma once

#include <alpaqa/util/demangled-typename.hpp>
#include <functional>
#include <map>
#include <typeinfo>
#include <variant>

#include <pybind11/chrono.h>
#include <pybind11/detail/typeid.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
namespace py = pybind11;

template <class T>
T dict_to_struct(const py::dict &);

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
            if (py::isinstance<py::dict>(h))
                t.*attr = dict_to_struct<A>(py::cast<py::dict>(h));
            else
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
using dict_to_struct_table_t = std::map<std::string, attr_setter_fun_t<T>>;

template <class T>
struct dict_to_struct_table {};

template <class T>
    requires requires { dict_to_struct_table<T>::table; }
void dict_to_struct_helper(T &t, const py::dict &kwargs) {
    const auto &m = dict_to_struct_table<T>::table;
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
    requires(!requires { dict_to_struct_table<T>::table; })
void dict_to_struct_helper(T &, const py::dict &) {
    throw std::runtime_error("No known conversion from Python dict to C++ type '" +
                             demangled_typename(typeid(T)) + '\'');
}

template <class T>
py::dict struct_to_dict_helper(const T &t) {
    const auto &m = dict_to_struct_table<T>::table;
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
T dict_to_struct(const py::dict &kwargs) {
    T t{};
    dict_to_struct_helper(t, kwargs);
    return t;
}

template <class T>
T kwargs_to_struct(const py::kwargs &kwargs) {
    return dict_to_struct<T>(kwargs);
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

#if 0
#include <alpaqa/inner/pga.hpp>

template <>
inline const dict_to_struct_table_t<alpaqa::PGAParams>
    dict_to_struct_table<alpaqa::PGAParams>{
        {"Lipschitz", &alpaqa::PGAParams::Lipschitz},
        {"max_iter", &alpaqa::PGAParams::max_iter},
        {"max_time", &alpaqa::PGAParams::max_time},
        {"L_min", &alpaqa::PGAParams::L_min},
        {"L_max", &alpaqa::PGAParams::L_max},
        {"stop_crit", &alpaqa::PGAParams::stop_crit},
        {"print_interval", &alpaqa::PGAParams::print_interval},
        {"quadratic_upperbound_tolerance_factor",
         &alpaqa::PGAParams::quadratic_upperbound_tolerance_factor},
        {"linesearch_tolerance_factor",
         &alpaqa::PGAParams::linesearch_tolerance_factor},
    };

#include <alpaqa/inner/guarded-aa-pga.hpp>

template <>
inline const dict_to_struct_table_t<alpaqa::GAAPGAParams>
    dict_to_struct_table<alpaqa::GAAPGAParams>{
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
        {"linesearch_tolerance_factor",
         &alpaqa::GAAPGAParams::linesearch_tolerance_factor},
        {"max_no_progress", &alpaqa::GAAPGAParams::max_no_progress},
        {"full_flush_on_γ_change",
         &alpaqa::GAAPGAParams::full_flush_on_γ_change},
    };

#endif

template <class T>
void make_dataclass(py::class_<T> &cls) {
    cls //
        .def(py::init(&dict_to_struct<T>))
        .def(py::init(&kwargs_to_struct<T>))
        .def("to_dict", &struct_to_dict<T>);
    for (auto &&[key, getset] : dict_to_struct_table<T>::table)
        cls.def_property(key.c_str(), getset.get, getset.set);
}

template <class T, class... Extra>
auto register_dataclass(py::handle scope, const char *name, const Extra &...extra) {
    py::class_<T> cls(scope, name, extra...);
    make_dataclass(cls);
    return cls;
}
