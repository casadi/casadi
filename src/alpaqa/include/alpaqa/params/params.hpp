#pragma once

#include <alpaqa/export.hpp>

#include <optional>
#include <span>
#include <stdexcept>
#include <string_view>

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

/// Overwrites @p t based on the @p options that start with @p prefix.
/// If @p used is not `nullopt`, sets corresponding flag of the options that
/// were used.
template <class T>
void ALPAQA_EXPORT set_params(
    T &t, std::string_view prefix, std::span<const std::string_view> options,
    std::optional<std::span<bool>> used = std::nullopt);

/// Custom parameter parsing exception.
struct ALPAQA_EXPORT invalid_param : std::invalid_argument {
    using std::invalid_argument::invalid_argument;
};

} // namespace alpaqa::params
