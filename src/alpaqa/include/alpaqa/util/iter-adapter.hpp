#pragma once

#include <concepts>
#include <iterator>
#include <ranges>
#include <type_traits>
#include <utility>

#if defined(__clang_major__) && __clang_major__ <= 15 && !defined(__clangd__)
#error "Better ranges support required"
#endif

namespace alpaqa::util {

template <class It>
struct iter_range_adapter {
    // P2325R3: https://www.open-std.org/jtc1/sc22/wg21/docs/papers/2021/p2325r3.html
    iter_range_adapter() = default;
    iter_range_adapter(It it) : it{std::forward<It>(it)} {}
    It it;

    struct sentinel_t {};

    struct iter_t : std::remove_cvref_t<It> {
        // P2325R3: https://www.open-std.org/jtc1/sc22/wg21/docs/papers/2021/p2325r3.html
        iter_t() = default;
        iter_t(It it) : std::remove_cvref_t<It>{std::forward<It>(it)} {}

        bool operator!=(sentinel_t) const { return static_cast<bool>(*this); }
        bool operator==(sentinel_t) const { return !static_cast<bool>(*this); }
        // TODO: For Clang bug
        friend bool operator!=(sentinel_t s, const iter_t &i) { return i != s; }
        friend bool operator==(sentinel_t s, const iter_t &i) { return i == s; }

        iter_t &operator++() {
            this->std::remove_cvref_t<It>::operator++();
            return *this;
        }
        iter_t operator++(int i) const {
            this->std::remove_cvref_t<It>::operator++(i);
            return *this;
        }
        const iter_t &operator*() const { return *this; }
        using value_type      = iter_t;
        using pointer         = iter_t *;
        using reference       = iter_t &;
        using difference_type = std::ptrdiff_t;
    };

    auto begin() const & -> std::input_or_output_iterator auto{
        return iter_t{it};
    }
    auto begin() && -> std::input_or_output_iterator auto{
        return iter_t{std::forward<It>(it)};
    }
    auto end() const -> std::sentinel_for<iter_t> auto{ return sentinel_t{}; }
};

} // namespace alpaqa::util

// Make iter_range_adapter available as an std::ranges view.
// Iterators remain valid even if the range is destroyed.
// Assume that the user takes care of lifetime, similar to std::span.
#ifndef DOXYGEN
template <class It>
inline constexpr bool ::std::ranges::enable_borrowed_range<
    alpaqa::util::iter_range_adapter<It>> = true;
#endif
