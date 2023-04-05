#pragma once

#include <chrono>

namespace alpaqa::util {
template <class T>
struct Timed {
    Timed(T &time) : time(time) {
        time -= std::chrono::steady_clock::now().time_since_epoch();
    }
    ~Timed() { time += std::chrono::steady_clock::now().time_since_epoch(); }
    Timed(const Timed &)            = delete;
    Timed(Timed &&)                 = delete;
    Timed &operator=(const Timed &) = delete;
    Timed &operator=(Timed &&)      = delete;
    T &time;
};
#ifndef DOXYGEN
template <class T>
Timed(T &) -> Timed<T>;
#endif
} // namespace alpaqa::util