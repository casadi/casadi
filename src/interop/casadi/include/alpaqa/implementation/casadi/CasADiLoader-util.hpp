#pragma once

#include <casadi/core/casadi_types.hpp>
#include <casadi/core/external.hpp>
#include <array>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>

namespace alpaqa::casadi_loader {

template <class F>
auto wrap_load(const std::string &so_name, const char *name, F f) {
    try {
        return f();
    } catch (const std::invalid_argument &e) {
        throw std::invalid_argument("Unable to load function '" + so_name +
                                    ":" + name + "': " + e.what());
    }
}

template <class T, class... Args>
auto wrapped_load(const std::string &so_name, const char *name,
                  Args &&...args) {
    return wrap_load(so_name, name, [&] {
        return T(casadi::external(name, so_name), std::forward<Args>(args)...);
    });
}

template <class T, class... Args>
std::optional<T> try_load(const std::string &so_name, const char *name,
                          Args &&...args) {
    try {
        return std::make_optional(
            wrapped_load<T>(so_name, name, std::forward<Args>(args)...));
    } catch (casadi::CasadiException &e) {
        return std::nullopt;
    }
}

using dim = std::pair<casadi_int, casadi_int>;
inline constexpr auto dims(auto... a) {
    if constexpr ((... && std::is_constructible_v<dim, decltype(a)>))
        return std::array{a...};
    else
        return std::array{dim{a, 1}...};
}

} // namespace alpaqa::casadi_loader