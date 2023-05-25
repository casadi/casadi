#pragma once

#include <alpaqa/params/params.hpp>
#include <alpaqa/util/demangled-typename.hpp>
#if __has_include(<charconv>)
#include <charconv>
#endif

#if __cpp_lib_to_chars && (!defined(__clang__)) &&                             \
    (!defined(__GNUC__) || __GNUC__ > 1)
#define ALPAQA_USE_FROM_CHARS_FLOAT 1
#else
#pragma message "Using std::stod as a fallback to replace std::from_chars"
#define ALPAQA_USE_FROM_CHARS_FLOAT 0
#endif
#define ALPAQA_USE_FROM_CHARS_INT 1

namespace {
using namespace alpaqa::params;

template <class T>
    requires(
#if ALPAQA_USE_FROM_CHARS_FLOAT
        std::floating_point<T> ||
#endif
#if ALPAQA_USE_FROM_CHARS_INT
        std::integral<T> ||
#endif
        false) // NOLINT(readability-simplify-boolean-expr)
const char *set_param_float_int(T &f, ParamString s) {
    const auto *val_end = s.value.data() + s.value.size();
    auto [ptr, ec]      = std::from_chars(s.value.data(), val_end, f);
    if (ec != std::errc())
        throw std::invalid_argument(
            "Invalid value '" + std::string(s.value) + "' for type '" +
            demangled_typename(typeid(T)) + "' in '" + std::string(s.full_key) +
            "': " + std::make_error_code(ec).message());
    return ptr;
}
template <class T>
    requires(
#if !ALPAQA_USE_FROM_CHARS_FLOAT
        std::floating_point<T> ||
#endif
#if !ALPAQA_USE_FROM_CHARS_INT
        std::integral<T> ||
#endif
        false) // NOLINT(readability-simplify-boolean-expr)
const char *set_param_float_int(T &f, ParamString s) {
    size_t end_index;
    try {
        if constexpr (std::is_same_v<T, float>)
            f = std::stof(std::string(s.value), &end_index);
        else if constexpr (std::is_same_v<T, double>)
            f = std::stod(std::string(s.value), &end_index);
        else if constexpr (std::is_same_v<T, long double>)
            f = std::stold(std::string(s.value), &end_index);
        else if constexpr (std::is_same_v<T, signed char>)
            f = static_cast<signed char>(
                std::stoi(std::string(s.value), &end_index, 0));
        else if constexpr (std::is_same_v<T, short>)
            f = static_cast<short>(
                std::stoi(std::string(s.value), &end_index, 0));
        else if constexpr (std::is_same_v<T, int>)
            f = std::stoi(std::string(s.value), &end_index, 0);
        else if constexpr (std::is_same_v<T, long>)
            f = std::stol(std::string(s.value), &end_index, 0);
        else if constexpr (std::is_same_v<T, long long>)
            f = std::stoll(std::string(s.value), &end_index, 0);
        else if constexpr (std::is_same_v<T, unsigned char>)
            f = static_cast<unsigned char>(
                std::stoul(std::string(s.value), &end_index, 0));
        else if constexpr (std::is_same_v<T, unsigned short>)
            f = static_cast<unsigned short>(
                std::stoul(std::string(s.value), &end_index, 0));
        else if constexpr (std::is_same_v<T, unsigned int>)
            f = static_cast<unsigned int>(
                std::stoul(std::string(s.value), &end_index, 0));
        else if constexpr (std::is_same_v<T, unsigned long>)
            f = std::stoul(std::string(s.value), &end_index, 0);
        else if constexpr (std::is_same_v<T, unsigned long long>)
            f = std::stoull(std::string(s.value), &end_index, 0);
        else
            static_assert(std::is_same_v<T, void>); // false
    } catch (std::exception &e) {
        throw std::invalid_argument("Invalid value '" + std::string(s.value) +
                                    "' for type '" +
                                    demangled_typename(typeid(T)) + "' in '" +
                                    std::string(s.full_key) + "': " + e.what());
    }
    return s.value.data() + end_index;
}
} // namespace