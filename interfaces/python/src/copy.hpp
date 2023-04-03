#pragma once

#include <pybind11/pybind11.h>
namespace py = pybind11;

template <class T, class... Args>
void default_copy(py::class_<T, Args...> &cls) {
    using namespace py::literals;
    cls.def("__copy__", [](const T &self) { return T{self}; });
}

template <class T, class... Args>
void default_deepcopy(py::class_<T, Args...> &cls) {
    using namespace py::literals;
    cls.def(
        "__deepcopy__", [](const T &self, py::dict) { return T{self}; }, "memo"_a);
}

template <class T, class... Args>
void default_copy_ctor(py::class_<T, Args...> &cls) {
    using namespace py::literals;
    cls.def(py::init<const T &>(), "other"_a, "Create a copy");
}

template <class T, class... Args>
void default_copy_methods(py::class_<T, Args...> &cls) {
    default_copy_ctor(cls);
    default_copy(cls);
    default_deepcopy(cls);
}