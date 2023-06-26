#pragma once

#include <alpaqa/export.hpp>
#include <alpaqa/util/float.hpp>

#include <iosfwd>
#include <limits>
#include <string>
#include <string_view>

#include <Eigen/Core>

namespace alpaqa {

template <std::floating_point F>
ALPAQA_EXPORT std::string
float_to_str(F value, int precision = std::numeric_limits<F>::max_digits10);

template <class T>
ALPAQA_EXPORT std::ostream &
print_csv_impl(std::ostream &os, const T &M, std::string_view sep = ",",
               std::string_view begin = "", std::string_view end = "\n");

template <class T>
ALPAQA_EXPORT std::ostream &print_matlab_impl(std::ostream &os, const T &M,
                                              std::string_view end = ";\n");

template <class T>
ALPAQA_EXPORT std::ostream &print_python_impl(std::ostream &os, const T &M,
                                              std::string_view end = "\n");

#define ALPAQA_PRINT_CSV_OVL_IMPL(type)                                        \
    inline std::ostream &print_csv(                                            \
        std::ostream &os, const Eigen::Ref<const Eigen::MatrixX<type>> &M,     \
        std::string_view sep, std::string_view begin, std::string_view end) {  \
        return print_csv_impl(os, M, sep, begin, end);                         \
    }                                                                          \
    inline std::ostream &print_csv(                                            \
        std::ostream &os, const Eigen::Ref<const Eigen::MatrixX<type>> &M,     \
        std::string_view sep, std::string_view begin) {                        \
        return print_csv_impl(os, M, sep, begin);                              \
    }                                                                          \
    inline std::ostream &print_csv(                                            \
        std::ostream &os, const Eigen::Ref<const Eigen::MatrixX<type>> &M,     \
        std::string_view sep) {                                                \
        return print_csv_impl(os, M, sep);                                     \
    }                                                                          \
    inline std::ostream &print_csv(                                            \
        std::ostream &os, const Eigen::Ref<const Eigen::MatrixX<type>> &M) {   \
        return print_csv_impl(os, M);                                          \
    }
#define ALPAQA_PRINT_OVL_IMPL(name, type)                                      \
    inline std::ostream &print_##name(                                         \
        std::ostream &os, const Eigen::Ref<const Eigen::MatrixX<type>> &M,     \
        std::string_view end) {                                                \
        return print_##name##_impl(os, M, end);                                \
    }                                                                          \
    inline std::ostream &print_##name(                                         \
        std::ostream &os, const Eigen::Ref<const Eigen::MatrixX<type>> &M) {   \
        return print_##name##_impl(os, M);                                     \
    }
#define ALPAQA_PRINT_OVL(type)                                                 \
    ALPAQA_PRINT_CSV_OVL_IMPL(type)                                            \
    ALPAQA_PRINT_OVL_IMPL(matlab, type)                                        \
    ALPAQA_PRINT_OVL_IMPL(python, type)

ALPAQA_PRINT_OVL(float)
ALPAQA_PRINT_OVL(double)
ALPAQA_PRINT_OVL(long double)
#ifdef ALPAQA_WITH_QUAD_PRECISION
ALPAQA_PRINT_OVL(__float128)
#endif

#undef ALPAQA_PRINT_OVL
#undef ALPAQA_PRINT_OVL_IMPL
#undef ALPAQA_PRINT_CSV_OVL_IMPL

} // namespace alpaqa