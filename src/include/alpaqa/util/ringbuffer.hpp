#pragma once

#include <cstddef>
#include <iterator>

namespace alpaqa {

template <class IndexT = size_t>
struct CircularIndices {
    using Index = IndexT;
    CircularIndices(Index zerobased, Index circular)
        : zerobased(zerobased), circular(circular) {}
    Index zerobased;
    Index circular;
};
/// @related    CircularIndices
/// @note   Only valid for two indices in the same range.
template <class IndexT>
bool operator==(CircularIndices<IndexT> a, CircularIndices<IndexT> b) {
    return a.zerobased == b.zerobased;
}
/// @related    CircularIndices
/// @note   Only valid for two indices in the same range.
template <class IndexT>
bool operator!=(CircularIndices<IndexT> a, CircularIndices<IndexT> b) {
    return !(a == b);
}

template <class IndexT = size_t>
struct CircularIndexIterator {
    using Index   = IndexT;
    using Indices = CircularIndices<Index>;

    CircularIndexIterator() : i{0, 0}, max{0} {}
    CircularIndexIterator(Indices i, Index max) : i(i), max(max) {}

    Indices i;
    Index max;

    using value_type        = Indices;
    using reference         = value_type;
    using difference_type   = std::ptrdiff_t; // This is required but not used
    using pointer           = void;
    using iterator_category = std::input_iterator_tag;

    reference operator*() const { return i; }
    CircularIndexIterator &operator++() {
        assert(i.zerobased < max);
        ++i.zerobased;
        i.circular = i.circular + 1 == max ? Index{0} : i.circular + 1;
        return *this;
    }
    CircularIndexIterator &operator--() {
        assert(i.zerobased > 0);
        --i.zerobased;
        i.circular = i.circular == Index{0} ? max - 1 : i.circular - 1;
        return *this;
    }
    CircularIndexIterator operator++(int) {
        auto r = *this;
        ++(*this);
        return r;
    }
    CircularIndexIterator operator--(int) {
        auto r = *this;
        --(*this);
        return r;
    }
};

/// @related    CircularIndexIterator
/// @note   Only valid for two indices in the same range.
template <class IndexT>
bool operator==(CircularIndexIterator<IndexT> a,
                CircularIndexIterator<IndexT> b) {
    assert(a.max == b.max);
    return a.i == b.i;
}
/// @related    CircularIndexIterator
/// @note   Only valid for two indices in the same range.
template <class IndexT>
bool operator!=(CircularIndexIterator<IndexT> a,
                CircularIndexIterator<IndexT> b) {
    return !(a == b);
}

template <class IndexT = size_t>
struct ReverseCircularIndexIterator {
    using ForwardIterator = CircularIndexIterator<IndexT>;
    using Index           = typename ForwardIterator::Index;
    using Indices         = typename ForwardIterator::Indices;

    ReverseCircularIndexIterator() : forwardit() {}
    ReverseCircularIndexIterator(Indices i, Index max) : forwardit(i, max) {}
    ReverseCircularIndexIterator(ForwardIterator forwardit)
        : forwardit(forwardit) {}

    ForwardIterator forwardit;

    using value_type        = Indices;
    using reference         = value_type;
    using difference_type   = std::ptrdiff_t; // This is required but not used
    using pointer           = void;
    using iterator_category = std::input_iterator_tag;

    reference operator*() const {
        auto tmp = forwardit;
        return *(--tmp);
    }
    ReverseCircularIndexIterator &operator++() {
        --forwardit;
        return *this;
    }
    ReverseCircularIndexIterator &operator--() {
        ++forwardit;
        return *this;
    }
    ReverseCircularIndexIterator operator++(int) {
        auto r = *this;
        ++(*this);
        return r;
    }
    ReverseCircularIndexIterator operator--(int) {
        auto r = *this;
        --(*this);
        return r;
    }
};

/// @related    ReverseCircularIndexIterator
/// @note   Only valid for two indices in the same range.
template <class IndexT>
bool operator==(ReverseCircularIndexIterator<IndexT> a,
                ReverseCircularIndexIterator<IndexT> b) {
    return a.forwardit == b.forwardit;
}
/// @related    ReverseCircularIndexIterator
/// @note   Only valid for two indices in the same range.
template <class IndexT>
bool operator!=(ReverseCircularIndexIterator<IndexT> a,
                ReverseCircularIndexIterator<IndexT> b) {
    return !(a == b);
}

template <class IndexT>
class CircularRange {
  public:
    using Index   = IndexT;
    using Indices = CircularIndices<Index>;

    CircularRange(Index size, Index idx1, Index idx2, Index max)
        : size(size), idx1(idx1), idx2(idx2), max(max) {}

    using const_iterator = CircularIndexIterator<Index>;
    using iterator       = const_iterator;

    using const_reverse_iterator = ReverseCircularIndexIterator<Index>;
    using reverse_iterator       = const_reverse_iterator;

    iterator begin() const { return {{Index{0}, idx1}, max}; }
    iterator end() const { return {{size, idx2}, max}; }
    const_iterator cbegin() const { return begin(); }
    const_iterator cend() const { return end(); }

    reverse_iterator rbegin() const { return reverse_iterator{end()}; }
    reverse_iterator rend() const { return reverse_iterator{begin()}; }
    const_reverse_iterator crbegin() const {
        return const_reverse_iterator{end()};
    }
    const_reverse_iterator crend() const {
        return const_reverse_iterator{begin()};
    }

  private:
    Index size;
    Index idx1, idx2;
    Index max;
};

template <class IndexT>
class ReverseCircularRange {
  public:
    using ForwardRange = CircularRange<IndexT>;
    using Index        = typename ForwardRange::Index;
    using Indices      = typename ForwardRange::Indices;

    ReverseCircularRange(const ForwardRange &forwardrange)
        : forwardrange(forwardrange) {}
    ReverseCircularRange(Index size, Index idx1, Index idx2, Index max)
        : forwardrange(size, idx1, idx2, max) {}

    using const_iterator = typename ForwardRange::const_reverse_iterator;
    using iterator       = typename ForwardRange::reverse_iterator;

    using const_reverse_iterator = typename ForwardRange::const_iterator;
    using reverse_iterator       = typename ForwardRange::iterator;

    iterator begin() const { return forwardrange.rbegin(); }
    iterator end() const { return forwardrange.rend(); }
    const_iterator cbegin() const { return forwardrange.crbegin(); }
    const_iterator cend() const { return forwardrange.crend(); }

    reverse_iterator rbegin() const { return forwardrange.begin(); }
    reverse_iterator rend() const { return forwardrange.end(); }
    const_reverse_iterator crbegin() const { return forwardrange.cbegin(); }
    const_reverse_iterator crend() const { return forwardrange.cend(); }

  private:
    ForwardRange forwardrange;
};

} // namespace alpaqa