#pragma once

namespace alpaqa::util {

/// Deleter for `std::unique_ptr` that just destructs the object, without
/// deallocating.
template <class T>
struct noop_delete {
    constexpr noop_delete() noexcept = default;
    template <class U>
    constexpr noop_delete(const noop_delete<U> &) noexcept {}
    constexpr void operator()(T *t) const noexcept { t->~T(); }
};

} // namespace alpaqa::util