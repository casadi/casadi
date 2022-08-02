#pragma once

#include <alpaqa/util/print.hpp>

#include <charconv>
#include <cmath>
#include <limits>
#include <ostream>
#include <string_view>

namespace alpaqa {

#if __cpp_lib_to_chars
std::string_view float_to_str_vw(
    auto &buf, std::floating_point auto value,
    int precision = std::numeric_limits<decltype(value)>::max_digits10) {
    auto begin = buf.data();
    if (!std::signbit(value))
        *begin++ = '+';
    auto [end, _] = std::to_chars(begin, buf.data() + buf.size(), value,
                                  std::chars_format::scientific, precision);
    return std::string_view{buf.data(), end};
}
#else
#pragma message "Using snprintf as a fallback to replace std::to_chars"

inline std::string_view float_to_str_vw_snprintf(auto &buf,
                                                 std::floating_point auto value,
                                                 int prec, const char *fmt) {
    int n = std::snprintf(buf.data(), buf.size(), fmt, prec, value);
    assert((size_t)n < buf.size());
    return {buf.data(), buf.data() + n};
}

inline std::string_view float_to_str_vw(
    auto &buf, double value,
    int precision = std::numeric_limits<decltype(value)>::max_digits10) {
    return float_to_str_vw_snprintf(buf, value, precision, "%+-#.*e");
}
inline std::string_view float_to_str_vw(
    auto &buf, float value,
    int precision = std::numeric_limits<decltype(value)>::max_digits10) {
    return float_to_str_vw(buf, static_cast<double>(value), precision);
}
inline std::string_view float_to_str_vw(
    auto &buf, long double value,
    int precision = std::numeric_limits<decltype(value)>::max_digits10) {
    return float_to_str_vw_snprintf(buf, value, precision, "%+-#.*Le");
}
#endif

#ifdef ALPAQA_WITH_QUAD_PRECISION
std::string_view float_to_str_vw(auto &buf, __float128 value) {
    int prec = std::numeric_limits<decltype(value)>::max_digits10;
    int n = quadmath_snprintf(buf.data(), buf.size(), "%+-#.*Qe", prec, value);
    assert((size_t)n < buf.size());
    return {buf.data(), buf.data() + n};
}
#endif

std::string float_to_str(std::floating_point auto value) {
    std::array<char, 64> buf;
    return std::string{float_to_str_vw(buf, value)};
}

template <std::floating_point F>
void print_elem_matlab(auto &buf, F value, std::ostream &os) {
    os << float_to_str_vw(buf, value);
}

template <std::floating_point F>
void print_elem_matlab(auto &buf, std::complex<F> value, std::ostream &os) {
    os << float_to_str_vw(buf, value.real()) << " + "
       << float_to_str_vw(buf, value.imag()) << 'i';
}

template <float_or_complex_float F>
std::ostream &print_matlab(std::ostream &os, const Eigen::MatrixX<F> &M) {
    os << '[';
    std::array<char, 64> buf;
    for (decltype(M.rows()) r{}; r < M.rows(); ++r) {
        for (decltype(M.cols()) c{}; c < M.cols(); ++c) {
            print_elem_matlab(buf, M(r, c), os);
            os << ' ';
        }
        if (r != M.rows() - 1)
            os << ";\n ";
    }
    return os << "];\n";
}

template <float_or_complex_float F>
std::ostream &print_matlab(std::ostream &os, const Eigen::VectorX<F> &M) {
    os << '[';
    std::array<char, 64> buf{};
    for (decltype(M.rows()) r{}; r < M.rows(); ++r) {
        print_elem_matlab(buf, M(r), os);
        if (r != M.rows() - 1)
            os << "; ";
    }
    return os << "];\n";
}

} // namespace alpaqa