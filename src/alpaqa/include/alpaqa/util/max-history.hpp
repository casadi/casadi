#pragma once

#include <algorithm>
#include <cstddef>
#include <vector>

namespace alpaqa {

/// Keep track of the maximum value over a specified horizon length.
template <class T>
class MaxHistory {
  public:
    MaxHistory(size_t memory) : buffer(std::vector<T>(memory)) {}

    void add(T newt) {
        if (full) {
            T oldt = std::move(*it);
            *it    = std::move(newt);
            if (*it > max_)
                max_ = *it;
            else if (oldt == max_)
                max_ = *std::max_element(buffer.begin(), buffer.end());
            ++it;
            if (it == buffer.end())
                it = buffer.begin();
        } else {
            if (it == buffer.begin() || newt > max_)
                max_ = newt;
            *it = std::move(newt);
            ++it;
            if (it == buffer.end()) {
                it   = buffer.begin();
                full = true;
            }
        }
    }

    const T &max() const { return max_; }

  private:
    std::vector<T> buffer;
    bool full                   = false;
    decltype(buffer.begin()) it = buffer.begin();
    T max_{};
};

} // namespace alpaqa