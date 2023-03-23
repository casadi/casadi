#include <alpaqa/params/params.hpp>

#include <alpaqa/config/config.hpp>
#include <alpaqa/util/demangled-typename.hpp>

#include <charconv>
#include <chrono>
#include <concepts>
#include <map>
#include <numeric>
#include <stdexcept>
#include <string_view>
#include <system_error>
#include <tuple>
#include <type_traits>

namespace alpaqa::params {

using config_t = DefaultConfig;

/// Throw a meaningful error when `s.key` is not empty, to indicate that
/// the given type @p T is not of struct type and cannot be indexed into.
template <class T>
void assert_key_empty(ParamString s) {
    if (!s.key.empty())
        throw invalid_param("Type '" + demangled_typename(typeid(T)) +
                            "' cannot be indexed in '" +
                            std::string(s.full_key) + "'");
}

/// Split the string @p full on the first occurrence of @p tok.
/// Returns (s, "") if tok was not found.
inline auto split_key(std::string_view full, char tok = '.') {
    auto tok_pos = full.find(tok);
    if (tok_pos == full.npos)
        return std::make_tuple(full, std::string_view{});
    std::string_view key{full.begin(), full.begin() + tok_pos};
    std::string_view rem{full.begin() + tok_pos + 1, full.end()};
    return std::make_tuple(key, rem);
}

/// Update/overwrite the first argument based on the option in @p s.
template <class T>
void set_param(T &,
               [[maybe_unused]] ParamString s); // deliberately undefined

/// Throw a meaningful error to indicate that parameters of type @p T are not
/// supported or implemented.
template <class T>
void unsupported_type(T &, [[maybe_unused]] ParamString s) {
    throw invalid_param("Unknown parameter type '" +
                        demangled_typename(typeid(T)) + "' in '" +
                        std::string(s.full_key) + "'");
}

template <>
void ALPAQA_EXPORT set_param(bool &b, ParamString s);

template <>
void ALPAQA_EXPORT set_param(std::string_view &v, ParamString s);

template <>
void ALPAQA_EXPORT set_param(std::string &v, ParamString s);

template <class T>
    requires((std::floating_point<T> || std::integral<T>) && !std::is_enum_v<T>)
void ALPAQA_EXPORT set_param(T &f, ParamString s);

template <>
void ALPAQA_EXPORT set_param(vec<config_t> &v, ParamString s);

template <class Rep, class Period>
void ALPAQA_EXPORT set_param(std::chrono::duration<Rep, Period> &t,
                             ParamString s);

/// Return a function that applies @ref set_param to the given attribute of a
/// value of type @p T.
template <class T, class A>
auto param_setter(A T::*attr) {
    return [attr](T &t, ParamString s) { return set_param(t.*attr, s); };
}

/// Function wrapper to set attributes of a struct, type-erasing the type of the
/// attribute.
template <class T>
struct param_setter_fun_t {
    template <class A>
    param_setter_fun_t(A T::*attr) : set(param_setter(attr)) {}
    std::function<void(T &, ParamString)> set;
    void operator()(T &t, ParamString s) const { return set(t, s); }
};

/// Dictionary that maps struct attribute names to type-erased functions that
/// set those attributes.
template <class T>
using dict_to_struct_table_t =
    std::map<std::string_view, param_setter_fun_t<T>>;

/// Specialize this type to define the attribute name to attribute setters
/// dictionaries for a struct type @p T.
template <class T>
struct dict_to_struct_table {};

/// Return a string enumerating the possible attribute names for the struct type
/// @p T.
template <class T>
auto possible_keys() {
    const auto &tbl = dict_to_struct_table<T>::table;
    if (tbl.empty())
        return std::string{};
    auto penult       = std::prev(tbl.end());
    auto quote_concat = [](std::string &&a, auto b) {
        return a + "'" + std::string(b.first) + "', ";
    };
    return std::accumulate(tbl.begin(), penult, std::string{}, quote_concat) +
           "'" + std::string(penult->first) + "'";
}

/// Use @p s to index into the struct type @p T and overwrite the attribute
/// given by @p s.key.
template <class T>
    requires requires { dict_to_struct_table<T>::table; }
void ALPAQA_EXPORT set_param(T &t, ParamString s) {
    const auto &m         = dict_to_struct_table<T>::table;
    auto [key, remainder] = split_key(s.key);
    auto it               = m.find(key);
    if (it == m.end())
        throw invalid_param("Invalid key '" + std::string(key) +
                            "' for type '" + demangled_typename(typeid(T)) +
                            "' in '" + std::string(s.full_key) +
                            "',\n  possible keys are: " + possible_keys<T>());
    s.key = remainder;
    it->second.set(t, s);
}

/// Helper macro to easily specialize @ref dict_to_struct_table.
#define PARAMS_TABLE(type_, ...)                                               \
    template <>                                                                \
    struct dict_to_struct_table<type_> {                                       \
        using type = type_;                                                    \
        inline static const dict_to_struct_table_t<type> table{__VA_ARGS__};   \
    }

/// Helper macro to easily initialize a @ref dict_to_struct_table_t.
#define PARAMS_MEMBER(name)                                                    \
    {                                                                          \
#name, &type::name                                                     \
    }

template <class T>
void ALPAQA_EXPORT set_params(T &t, std::string_view prefix,
                              std::span<const std::string_view> options,
                              std::optional<std::span<bool>> used) {
    size_t index = 0;
    for (const auto &kv : options) {
        auto [key, value]     = split_key(kv, '=');
        auto [pfx, remainder] = split_key(key);
        auto curr_index       = index++;
        if (pfx != prefix)
            continue;
        if (used)
            (*used)[curr_index] = true;
        set_param(t, {.full_key = kv, .key = remainder, .value = value});
    }
}

} // namespace alpaqa::params