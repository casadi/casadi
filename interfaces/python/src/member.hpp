#pragma once

#include <pybind11/pybind11.h>
namespace py = pybind11;

/// Exposes a member @p M by reference, with
/// `py::return_value_policy::reference_internal`.
template <auto M>
auto member_ref() {
    return []<class T, class A>(A T::*) {
        return py::cpp_function([](T &self) -> A & { return self.*M; },
                                py::return_value_policy::reference_internal);
    }(M);
}

/// Exposes a member @p M by pointer, with
/// `py::return_value_policy::reference_internal`.
template <auto M>
auto member_ptr() {
    return []<class T, class A>(A *T::*) {
        return py::cpp_function([](T &self) -> A * { return self.*M; },
                                py::return_value_policy::reference_internal);
    }(M);
}

/// Exposes a member @p M by reference as `rvec`, with
/// `py::return_value_policy::reference_internal`.
template <auto M>
auto vector_getter() {
    return []<class T, class A>(A T::*) {
        return py::cpp_function([](T &self) -> typename T::rvec { return self.*M; },
                                py::return_value_policy::reference_internal);
    }(M);
}

/// Returns a setter for a member @p M, ensuring that the size of the new vector
/// matches the size of the old one.
/// @throws std::invalid_argument
///         If the sizes don't match.
/// @param  name
///         Name to use in the exception's error message.
template <auto M>
auto vector_setter(std::string_view name) {
    return [name]<class T, class A>(A T::*) {
        return py::cpp_function([name](T &self, typename T::crvec value) {
            if (value.size() != (self.*M).size())
                throw std::invalid_argument("Invalid dimension for '" + std::string(name) +
                                            "': got " + std::to_string(value.size()) +
                                            ", should be " + std::to_string((self.*M).size()) +
                                            ".");
            self.*M = value;
        });
    }(M);
}
