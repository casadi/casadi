#pragma once

#include <ranges>
#include <utility>

namespace alpaqa::util {

template <class Rng>
struct enumerate_t : std::ranges::view_interface<enumerate_t<Rng>> {
    Rng rng;

    // P2325R3: https://www.open-std.org/jtc1/sc22/wg21/docs/papers/2021/p2325r3.html
    enumerate_t() = default;
    enumerate_t(Rng rng) : rng{std::forward<Rng>(rng)} {}

    using begin_t = decltype(std::ranges::begin(std::as_const(rng)));
    using end_t   = decltype(std::ranges::end(std::as_const(rng)));

    struct sentinel_t {
        end_t it;
    };

    struct iter_t {
        // P2325R3: https://www.open-std.org/jtc1/sc22/wg21/docs/papers/2021/p2325r3.html
        iter_t() = default;
        iter_t(begin_t &&it) : it{std::forward<begin_t>(it)} {}

        using index_t = std::ranges::range_difference_t<Rng>;
        index_t index{};
        begin_t it;

        using difference_type = std::ptrdiff_t;
        using value_type      = std::tuple<index_t, decltype(*it)>;

        bool operator!=(sentinel_t s) const { return s.it != it; }
        bool operator==(sentinel_t s) const { return s.it == it; }
        // TODO: For Clang bug
        friend bool operator!=(sentinel_t s, const iter_t &i) { return i != s; }
        friend bool operator==(sentinel_t s, const iter_t &i) { return i == s; }

        iter_t &operator++() {
            ++it;
            ++index;
            return *this;
        }
        iter_t operator++(int) const {
            auto tmp = *this;
            ++*this;
            return tmp;
        }

        value_type operator*() const { return {index, *it}; }
    };

    auto begin() const -> std::input_or_output_iterator auto{
        return iter_t{std::ranges::begin(rng)};
    }
    auto end() const
    // -> std::sentinel_for<iter_t> auto
    {
        return sentinel_t{std::ranges::end(rng)};
    }
};

template <class Rng>
auto enumerate(Rng &&rng) {
    return enumerate_t<std::views::all_t<Rng>>{std::forward<Rng>(rng)};
}

} // namespace alpaqa::util
