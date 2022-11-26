#pragma once

#include <pybind11/iostream.h>
#include <utility>

template <class T>
class StreamReplacer {
    pybind11::detail::pythonbuf buffer{pybind11::module_::import("sys").attr("stdout")};
    std::ostream os{&buffer};
    T *ptr                    = nullptr;
    std::ostream *original_os = nullptr;

  public:
    StreamReplacer(T *ptr) : ptr(ptr) { original_os = std::exchange(ptr->os, &os); }
    ~StreamReplacer() { ptr->os = original_os; }
    StreamReplacer(const StreamReplacer &)            = delete;
    StreamReplacer &operator=(const StreamReplacer &) = delete;
    StreamReplacer(const StreamReplacer &&)           = delete;
    StreamReplacer &operator=(StreamReplacer &&)      = delete;
};
