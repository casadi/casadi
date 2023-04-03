#pragma once

#include <alpaqa/util/demangled-typename.hpp>
#include <optional>
#include <set>
#include <stdexcept>

template <class T>
class ThreadChecker {
    using set_t      = std::set<const T *>;
    using iterator_t = typename set_t::iterator;
    static set_t set;
    std::optional<iterator_t> iterator;

  public:
    ThreadChecker(const T *ptr) {
        auto [iter, inserted] = set.insert(ptr);
        if (!inserted) {
            std::string name = "instance of type " + demangled_typename(typeid(T));
            if constexpr (requires { ptr->get_name(); })
                name = "instance of " + std::string(ptr->get_name());
            throw std::runtime_error("Same " + name +
                                     " used in multiple threads (consider making a copy)");
        }
        iterator = iter;
    }
    ~ThreadChecker() {
        if (iterator)
            set.erase(*iterator);
    }
    ThreadChecker(const ThreadChecker &)            = delete;
    ThreadChecker &operator=(const ThreadChecker &) = delete;
    ThreadChecker(ThreadChecker &&o) noexcept { std::swap(this->iterator, o.iterator); }
    ThreadChecker &operator=(ThreadChecker &&o) noexcept {
        this->iterator = std::move(o.iterator);
        o.iterator.reset();
        return *this;
    }
};

template <class T>
typename ThreadChecker<T>::set_t ThreadChecker<T>::set;
