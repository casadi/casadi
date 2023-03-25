#pragma once

#include <alpaqa/config/config.hpp>
#include <alpaqa/export.h>
#include <iosfwd>
#include <stdexcept>
#include <vector>

namespace alpaqa::csv {

struct ALPAQA_EXPORT read_error : std::runtime_error {
    using std::runtime_error::runtime_error;
};

template <std::floating_point F>
void ALPAQA_EXPORT read_row_impl(std::istream &is,
                                 Eigen::Ref<Eigen::VectorX<F>> v,
                                 char sep = ',');

template <std::floating_point F>
std::vector<F> ALPAQA_EXPORT read_row_std_vector(std::istream &is,
                                                 char sep = ',');

#define ALPAQA_READ_ROW_OVL(type)                                              \
    inline void read_row(std::istream &is, Eigen::Ref<Eigen::VectorX<type>> v, \
                         char sep) {                                           \
        return read_row_impl<type>(is, v, sep);                                \
    }                                                                          \
    inline void read_row(std::istream &is,                                     \
                         Eigen::Ref<Eigen::VectorX<type>> v) {                 \
        return read_row_impl<type>(is, v);                                     \
    }

ALPAQA_READ_ROW_OVL(float)
ALPAQA_READ_ROW_OVL(double)
ALPAQA_READ_ROW_OVL(long double)
// #ifdef ALPAQA_WITH_QUAD_PRECISION
// ALPAQA_READ_ROW_OVL(__float128)
// #endif

#undef ALPAQA_READ_ROW_OVL

} // namespace alpaqa::csv