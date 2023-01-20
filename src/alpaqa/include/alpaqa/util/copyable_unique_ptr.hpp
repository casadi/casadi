#pragma once

#include <memory>

namespace alpaqa::util {

template <class T>
struct copyable_unique_ptr {
    copyable_unique_ptr(std::unique_ptr<T> ptr) : ptr{std::move(ptr)} {}
    copyable_unique_ptr() = default;
    copyable_unique_ptr(const copyable_unique_ptr &o)
        : ptr{o.ptr ? std::make_unique<T>(*o.ptr) : nullptr} {}
    copyable_unique_ptr &operator=(const copyable_unique_ptr &o) {
        this->ptr = o.ptr ? std::make_unique<T>(*o.ptr) : nullptr;
        return *this;
    }
    copyable_unique_ptr(copyable_unique_ptr &&) noexcept            = default;
    copyable_unique_ptr &operator=(copyable_unique_ptr &&) noexcept = default;

    operator std::unique_ptr<T> &() & { return ptr; }
    operator const std::unique_ptr<T> &() const & { return ptr; }
    operator std::unique_ptr<T> &&() && { return std::move(ptr); }

    std::unique_ptr<T> &operator->() { return ptr; }
    const std::unique_ptr<T> &operator->() const { return ptr; }
    auto &operator*() { return *ptr; }
    auto &operator*() const { return *ptr; }

    std::unique_ptr<T> ptr;
};

} // namespace alpaqa::util