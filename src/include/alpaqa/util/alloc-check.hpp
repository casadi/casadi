#pragma once

#include <Eigen/Core>

namespace alpaqa {

#ifdef EIGEN_RUNTIME_NO_MALLOC
template <bool Allow>
struct ScopedMallocChecker {
    ScopedMallocChecker() : prev(Eigen::internal::is_malloc_allowed()) {
        Eigen::internal::set_is_malloc_allowed(Allow);
    }
    ~ScopedMallocChecker() { Eigen::internal::set_is_malloc_allowed(prev); }
    bool prev;
};
#else
template <bool>
struct [[maybe_unused]] ScopedMallocChecker {};
#endif
struct [[maybe_unused]] ScopedMallocBlocker : ScopedMallocChecker<false> {};
struct [[maybe_unused]] ScopedMallocAllower : ScopedMallocChecker<true> {};

} // namespace alpaqa