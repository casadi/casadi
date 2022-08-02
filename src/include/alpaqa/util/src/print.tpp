#pragma once

#include <alpaqa/util/print.hpp>

#include <charconv>
#include <cmath>
#include <limits>
#include <ostream>
#include <string_view>

namespace alpaqa {

inline std::string_view float_to_str_vw_snprintf(auto &&print, auto &buf,
                                                 std::floating_point auto value,
                                                 int precision,
                                                 const char *fmt) {
    int n = print(buf.data(), buf.size(), fmt, precision, value);
    assert((size_t)n < buf.size());
    return {buf.data(), buf.data() + n};
}

#if __cpp_lib_to_chars
template <std::floating_point F>
std::string_view
float_to_str_vw(auto &buf, F value,
                int precision = std::numeric_limits<F>::max_digits10) {
    auto begin = buf.data();
    if (!std::signbit(value))
        *begin++ = '+';
    auto [end, _] = std::to_chars(begin, buf.data() + buf.size(), value,
                                  std::chars_format::scientific, precision);
    return std::string_view{buf.data(), end};
}
#else
#pragma message "Using snprintf as a fallback to replace std::to_chars"

inline std::string_view
float_to_str_vw(auto &buf, double value,
                int precision = std::numeric_limits<double>::max_digits10) {
    return float_to_str_vw_snprintf(std::snprintf, buf, value, precision,
                                    "%+-#.*e");
}
inline std::string_view
float_to_str_vw(auto &buf, float value,
                int precision = std::numeric_limits<float>::max_digits10) {
    return float_to_str_vw(buf, static_cast<double>(value), precision);
}
inline std::string_view float_to_str_vw(
    auto &buf, long double value,
    int precision = std::numeric_limits<long double>::max_digits10) {
    return float_to_str_vw_snprintf(std::snprintf, buf, value, precision,
                                    "%+-#.*Le");
}
#endif

#ifdef ALPAQA_WITH_QUAD_PRECISION
std::string_view
float_to_str_vw(auto &buf, __float128 value,
                int precision = std::numeric_limits<__float128>::max_digits10) {
    return float_to_str_vw_snprintf(quadmath_snprintf, buf, value, precision,
                                    "%+-#.*Qe");
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