#include <algorithm>
#include <gtest/gtest.h>

#include <iterator>
#include <alpaqa/util/ringbuffer.hpp>
#include <type_traits>

auto circular  = [](alpaqa::CircularIndices<int> i) { return i.circular; };
auto zerobased = [](alpaqa::CircularIndices<int> i) { return i.zerobased; };

TEST(CircularRange, mutforward) {
    alpaqa::CircularRange<int> r{4, 3, 2, 5};
    std::vector<int> result_z, result_c;
    std::transform(r.begin(), r.end(), std::back_insert_iterator(result_z),
                   zerobased);
    std::transform(r.begin(), r.end(), std::back_insert_iterator(result_c),
                   circular);
    std::vector<int> expected_z{0, 1, 2, 3}, expected_c{3, 4, 0, 1};
    EXPECT_EQ(result_z, expected_z);
    EXPECT_EQ(result_c, expected_c);
}

TEST(CircularRange, constforward) {
    const alpaqa::CircularRange<int> r{4, 3, 2, 5};
    std::vector<int> result_z, result_c;
    std::transform(r.begin(), r.end(), std::back_insert_iterator(result_z),
                   zerobased);
    std::transform(r.begin(), r.end(), std::back_insert_iterator(result_c),
                   circular);
    std::vector<int> expected_z{0, 1, 2, 3}, expected_c{3, 4, 0, 1};
    EXPECT_EQ(result_z, expected_z);
    EXPECT_EQ(result_c, expected_c);
}

TEST(CircularRange, forwardfull) {
    alpaqa::CircularRange<int> r{5, 3, 3, 5};
    std::vector<int> result_z, result_c;
    std::transform(r.begin(), r.end(), std::back_insert_iterator(result_z),
                   zerobased);
    std::transform(r.begin(), r.end(), std::back_insert_iterator(result_c),
                   circular);
    std::vector<int> expected_z{0, 1, 2, 3, 4}, expected_c{3, 4, 0, 1, 2};
    EXPECT_EQ(result_z, expected_z);
    EXPECT_EQ(result_c, expected_c);
}

TEST(CircularRange, mutreverse) {
    alpaqa::CircularRange<int> r{4, 3, 2, 5};
    std::vector<int> result_z, result_c;
    std::transform(r.rbegin(), r.rend(), std::back_insert_iterator(result_z),
                   zerobased);
    std::transform(r.rbegin(), r.rend(), std::back_insert_iterator(result_c),
                   circular);
    std::vector<int> expected_z{0, 1, 2, 3}, expected_c{3, 4, 0, 1};
    std::reverse(expected_z.begin(), expected_z.end());
    std::reverse(expected_c.begin(), expected_c.end());
    EXPECT_EQ(result_z, expected_z);
    EXPECT_EQ(result_c, expected_c);
}

TEST(CircularRange, constreverse) {
    const alpaqa::CircularRange<int> r{4, 3, 2, 5};
    std::vector<int> result_z, result_c;
    std::transform(r.rbegin(), r.rend(), std::back_insert_iterator(result_z),
                   zerobased);
    std::transform(r.rbegin(), r.rend(), std::back_insert_iterator(result_c),
                   circular);
    std::vector<int> expected_z{0, 1, 2, 3}, expected_c{3, 4, 0, 1};
    std::reverse(expected_z.begin(), expected_z.end());
    std::reverse(expected_c.begin(), expected_c.end());
    EXPECT_EQ(result_z, expected_z);
    EXPECT_EQ(result_c, expected_c);
}

TEST(CircularRange, reversefull) {
    alpaqa::CircularRange<int> r{5, 3, 3, 5};
    std::vector<int> result_z, result_c;
    std::transform(r.rbegin(), r.rend(), std::back_insert_iterator(result_z),
                   zerobased);
    std::transform(r.rbegin(), r.rend(), std::back_insert_iterator(result_c),
                   circular);
    std::vector<int> expected_z{0, 1, 2, 3, 4}, expected_c{3, 4, 0, 1, 2};
    std::reverse(expected_z.begin(), expected_z.end());
    std::reverse(expected_c.begin(), expected_c.end());
    EXPECT_EQ(result_z, expected_z);
    EXPECT_EQ(result_c, expected_c);
}

TEST(CircularRange, iteratortraits) {
    alpaqa::CircularRange<int> r{5, 3, 3, 5};
    auto a  = std::next(r.begin());
    using I = decltype(a);

    // // https://en.cppreference.com/w/cpp/named_req/BidirectionalIterator
    // EXPECT_EQ(--(++a), a);
    // EXPECT_EQ(&(--a), &a);
    // EXPECT_TRUE((std::is_same_v<decltype(--a), I &>));
    // EXPECT_TRUE((std::is_convertible_v<decltype(a--), const I &>));
    // EXPECT_TRUE(
    //     (std::is_same_v<decltype(*a--), std::iterator_traits<I>::reference>));

    // // https://en.cppreference.com/w/cpp/named_req/ForwardIterator
    // EXPECT_EQ([](I a) { return a++; }(a),
    //           [](I a) {
    //               I ip = a;
    //               ++a;
    //               return ip;
    //           }(a));
    // EXPECT_TRUE(
    //     (std::is_lvalue_reference_v<std::iterator_traits<I>::reference>));
    // EXPECT_TRUE((std::is_same_v<std::remove_cv_t<std::remove_reference_t<
    //                                 std::iterator_traits<I>::reference>>,
    //                             std::iterator_traits<I>::value_type>));
    // EXPECT_TRUE((std::is_convertible_v<decltype(a++), const I &>));
    // EXPECT_TRUE(
    //     (std::is_same_v<decltype(*a++), std::iterator_traits<I>::reference>));
    // EXPECT_TRUE((std::is_default_constructible_v<I>));
    // // Note: multipass guarantee is not satisfied because a == b does not imply
    // //       that the references *a and *b are bound to the same object.

    // https://en.cppreference.com/w/cpp/named_req/InputIterator
    auto b = a;
    EXPECT_EQ(a, b);
    EXPECT_EQ(a != b, !(a == b));
    b = std::next(b);
    EXPECT_NE(a, b);
    EXPECT_EQ(a != b, !(a == b));
    // EXPECT_EQ(&(a->zerobased), &((*a).zerobased));
    EXPECT_EQ([](I a) { return *a++; }(a),
              [](I a) {
                  std::iterator_traits<I>::value_type x = *a;
                  ++a;
                  return x;
              }(a));
    EXPECT_TRUE(std::is_signed_v<std::iterator_traits<I>::difference_type>);

    // https://en.cppreference.com/w/cpp/named_req/Iterator
    EXPECT_TRUE(std::is_copy_constructible_v<I>);
    EXPECT_TRUE(std::is_copy_assignable_v<I>);
    EXPECT_TRUE(std::is_destructible_v<I>);
    EXPECT_TRUE(std::is_swappable_v<I>);
    EXPECT_TRUE((std::is_same_v<decltype(++a), I &>));
}

TEST(ReverseCircularRange, mutforward) {
    alpaqa::ReverseCircularRange<int> r{4, 3, 2, 5};
    std::vector<int> result_z, result_c;
    std::transform(r.begin(), r.end(), std::back_insert_iterator(result_z),
                   zerobased);
    std::transform(r.begin(), r.end(), std::back_insert_iterator(result_c),
                   circular);
    std::vector<int> expected_z{0, 1, 2, 3}, expected_c{3, 4, 0, 1};
    std::reverse(expected_z.begin(), expected_z.end());
    std::reverse(expected_c.begin(), expected_c.end());
    EXPECT_EQ(result_z, expected_z);
    EXPECT_EQ(result_c, expected_c);
}

TEST(ReverseCircularRange, constforward) {
    const alpaqa::ReverseCircularRange<int> r{4, 3, 2, 5};
    std::vector<int> result_z, result_c;
    std::transform(r.begin(), r.end(), std::back_insert_iterator(result_z),
                   zerobased);
    std::transform(r.begin(), r.end(), std::back_insert_iterator(result_c),
                   circular);
    std::vector<int> expected_z{0, 1, 2, 3}, expected_c{3, 4, 0, 1};
    std::reverse(expected_z.begin(), expected_z.end());
    std::reverse(expected_c.begin(), expected_c.end());
    EXPECT_EQ(result_z, expected_z);
    EXPECT_EQ(result_c, expected_c);
}

TEST(ReverseCircularRange, forwardfull) {
    alpaqa::ReverseCircularRange<int> r{5, 3, 3, 5};
    std::vector<int> result_z, result_c;
    std::transform(r.begin(), r.end(), std::back_insert_iterator(result_z),
                   zerobased);
    std::transform(r.begin(), r.end(), std::back_insert_iterator(result_c),
                   circular);
    std::vector<int> expected_z{0, 1, 2, 3, 4}, expected_c{3, 4, 0, 1, 2};
    std::reverse(expected_z.begin(), expected_z.end());
    std::reverse(expected_c.begin(), expected_c.end());
    EXPECT_EQ(result_z, expected_z);
    EXPECT_EQ(result_c, expected_c);
}

TEST(ReverseCircularRange, mutreverse) {
    alpaqa::ReverseCircularRange<int> r{4, 3, 2, 5};
    std::vector<int> result_z, result_c;
    std::transform(r.rbegin(), r.rend(), std::back_insert_iterator(result_z),
                   zerobased);
    std::transform(r.rbegin(), r.rend(), std::back_insert_iterator(result_c),
                   circular);
    std::vector<int> expected_z{0, 1, 2, 3}, expected_c{3, 4, 0, 1};
    EXPECT_EQ(result_z, expected_z);
    EXPECT_EQ(result_c, expected_c);
}

TEST(ReverseCircularRange, constreverse) {
    const alpaqa::ReverseCircularRange<int> r{4, 3, 2, 5};
    std::vector<int> result_z, result_c;
    std::transform(r.rbegin(), r.rend(), std::back_insert_iterator(result_z),
                   zerobased);
    std::transform(r.rbegin(), r.rend(), std::back_insert_iterator(result_c),
                   circular);
    std::vector<int> expected_z{0, 1, 2, 3}, expected_c{3, 4, 0, 1};
    EXPECT_EQ(result_z, expected_z);
    EXPECT_EQ(result_c, expected_c);
}

TEST(ReverseCircularRange, reversefull) {
    alpaqa::ReverseCircularRange<int> r{5, 3, 3, 5};
    std::vector<int> result_z, result_c;
    std::transform(r.rbegin(), r.rend(), std::back_insert_iterator(result_z),
                   zerobased);
    std::transform(r.rbegin(), r.rend(), std::back_insert_iterator(result_c),
                   circular);
    std::vector<int> expected_z{0, 1, 2, 3, 4}, expected_c{3, 4, 0, 1, 2};
    EXPECT_EQ(result_z, expected_z);
    EXPECT_EQ(result_c, expected_c);
}