#pragma once

#include <algorithm>
#include <functional>
#include <iterator>
#include <ranges>

namespace alpaqa::util {

template <std::ranges::input_range R1, std::ranges::input_range R2,
          class Comp = std::ranges::less, class Proj1 = std::identity,
          class Proj2 = std::identity>
    requires(std::ranges::view<R1> && std::ranges::view<R2>)
struct set_intersection_iterable
    : std::ranges::view_interface<
          set_intersection_iterable<R1, R2, Comp, Proj1, Proj2>> {

    R1 range1;
    R2 range2;
    Comp comp;
    Proj1 proj1;
    Proj2 proj2;

    // P2325R3: https://www.open-std.org/jtc1/sc22/wg21/docs/papers/2021/p2325r3.html
    set_intersection_iterable() = default;
    set_intersection_iterable(R1 &&range1, R2 &&range2, Comp &&comp,
                              Proj1 &&proj1, Proj2 &&proj2)
        : range1{std::forward<R1>(range1)}, range2{std::forward<R2>(range2)},
          comp{std::forward<Comp>(comp)}, proj1{std::forward<Proj1>(proj1)},
          proj2{std::forward<Proj2>(proj2)} {}

    struct sentinel_t {};
    template <std::input_iterator I1, std::sentinel_for<I1> S1,
              std::input_iterator I2, std::sentinel_for<I2> S2>
    struct iter_t {
        // P2325R3: https://www.open-std.org/jtc1/sc22/wg21/docs/papers/2021/p2325r3.html
        iter_t() = default;
        iter_t(I1 first1, S1 last1, I2 first2, S2 last2, Comp comp, Proj1 proj1,
               Proj2 proj2)
            : first1{std::move(first1)}, last1{std::move(last1)},
              first2{std::move(first2)}, last2{std::move(last2)},
              comp{std::move(comp)}, proj1{std::move(proj1)}, proj2{std::move(
                                                                  proj2)} {}
        I1 first1;
        S1 last1;
        I2 first2;
        S2 last2;
        Comp comp;
        Proj1 proj1;
        Proj2 proj2;

        using difference_type = std::ptrdiff_t; // TODO
        using value_type = std::tuple<decltype(*first1), decltype(*first2)>;

        bool operator!=(sentinel_t) const {
            return first1 != last1 && first2 != last2;
        }
        bool operator==(sentinel_t s) const { return !(*this != s); }
        // TODO: For Clang bug
        friend bool operator!=(sentinel_t s, const iter_t &i) { return i != s; }
        friend bool operator==(sentinel_t s, const iter_t &i) { return i == s; }

        iter_t &operator++() {
            ++first1, ++first2;
            advance();
            return *this;
        }
        iter_t operator++(int) {
            auto tmp = *this;
            ++*this;
            return tmp;
        }
        value_type operator*() const { return {*first1, *first2}; }
        void advance() {
            while (*this != sentinel_t{}) {
                if (std::invoke(comp, std::invoke(proj1, *first1),
                                std::invoke(proj2, *first2)))
                    ++first1;
                else if (std::invoke(comp, std::invoke(proj2, *first2),
                                     std::invoke(proj1, *first1)))
                    ++first2;
                else
                    break;
            }
        }
    };

  private:
    template <class I1, class S1, class I2, class S2>
    iter_t<I1, S1, I2, S2> iter(I1 first1, S1 last1, I2 first2,
                                S2 last2) const {
        return {first1, last1, first2, last2, comp, proj1, proj2};
    }

  public:
    auto begin() const -> std::input_or_output_iterator auto{
        auto it = iter(std::ranges::begin(range1), std::ranges::end(range1),
                       std::ranges::begin(range2), std::ranges::end(range2));
        it.advance();
        return it;
    }
    auto end() const
    // -> std::sentinel_for< decltype(std::declval<set_intersection_iterable>().begin())> auto
    {
        return sentinel_t{};
    }
};

template <std::ranges::viewable_range R1, std::ranges::viewable_range R2,
          class Comp = std::ranges::less, class Proj1 = std::identity,
          class Proj2 = std::identity>
set_intersection_iterable<std::ranges::views::all_t<R1>,
                          std::ranges::views::all_t<R2>, Comp, Proj1, Proj2>
iter_set_intersection(R1 &&r1, R2 &&r2, Comp comp = {}, Proj1 proj1 = {},
                      Proj2 proj2 = {}) {
    static_assert(
        requires(set_intersection_iterable<std::ranges::views::all_t<R1>,
                                           std::ranges::views::all_t<R2>, Comp,
                                           Proj1, Proj2>
                     s) {
            { s.end() } -> std::sentinel_for<decltype(s.begin())>;
        });
    return {
        std::forward<R1>(r1), std::forward<R2>(r2), std::move(comp),
        std::move(proj1),     std::move(proj2),
    };
}

} // namespace alpaqa::util