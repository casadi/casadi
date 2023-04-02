#pragma once

#include <pybind11/pybind11.h>
namespace py = pybind11;

template <auto M>
auto member_ref() {
    return []<class T, class A>(A T::*) {
        return py::cpp_function([](T &self) -> A & { return self.*M; },
                                py::return_value_policy::reference_internal);
    }(M);
}

template <auto M>
auto member_ptr() {
    return []<class T, class A>(A *T::*) {
        return py::cpp_function([](T &self) -> A * { return self.*M; },
                                py::return_value_policy::reference_internal);
    }(M);
}