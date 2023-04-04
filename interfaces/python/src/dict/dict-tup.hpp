#pragma once

#include "kwargs-to-struct.hpp"

#include <tuple>
#include <utility>

template <class P>
py::object to_dict_tup(const P &params) {
    return struct_to_dict(params);
}
template <class... Ps>
py::object to_dict_tup(const std::tuple<Ps...> tup) {
    return py::cast([&]<size_t... Is>(std::index_sequence<Is...>) {
        return std::make_tuple(to_dict_tup(std::get<Is>(tup))...);
    }(std::make_index_sequence<sizeof...(Ps)>()));
}
inline py::object to_dict_tup(py::object o) { return o; }
