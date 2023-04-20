#pragma once

#include <alpaqa/config/config.hpp>
#include <alpaqa/export.hpp>

#include <optional>
#include <span>
#include <stdexcept>
#include <string_view>
#include <tuple>

namespace alpaqa::params {

/// Represents a parameter value encoded as a string in the format
/// `abc.def.key=value`.
struct ALPAQA_EXPORT ParamString {
    /// Full key string, used for diagnostics.
    std::string_view full_key;
    /// The subkey to resolve next.
    std::string_view key;
    /// The value of the parameter to store.
    std::string_view value;
};

/// Custom parameter parsing exception.
struct ALPAQA_EXPORT invalid_param : std::invalid_argument {
    using std::invalid_argument::invalid_argument;
};

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
void set_param(T &, ParamString); /* deliberately undefined */

/// Overwrites @p t based on the @p options that start with @p prefix.
/// If @p used is not `nullopt`, sets corresponding flag of the options that
/// were used.
template <class T>
void ALPAQA_EXPORT set_params(
    T &t, std::string_view prefix, std::span<const std::string_view> options,
    std::optional<std::span<bool>> used = std::nullopt) {

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

template <Config Conf>
struct ALPAQA_EXPORT vec_from_file {
    USING_ALPAQA_CONFIG(Conf);
    length_t expected_size;
    std::optional<vec> value = std::nullopt;
};

} // namespace alpaqa::params
