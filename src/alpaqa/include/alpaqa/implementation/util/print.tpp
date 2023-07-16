#pragma once

#include <alpaqa/util/print.hpp>

#include <cassert>
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
    return {buf.data(), (size_t)n};
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
#pragma message "Using std::snprintf as a fallback to replace std::to_chars"

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

template <std::floating_point F>
std::string float_to_str(F value, int precision) {
    std::array<char, 64> buf;
    return std::string{float_to_str_vw(buf, value, precision)};
}

template <std::floating_point F>
void print_elem(auto &buf, F value, std::ostream &os) {
    os << float_to_str_vw(buf, value);
}

template <std::floating_point F>
void print_elem(auto &buf, std::complex<F> value, std::ostream &os) {
    os << float_to_str_vw(buf, value.real()) << " + "
       << float_to_str_vw(buf, value.imag()) << 'j';
}

template <class T>
std::ostream &print_csv_impl(std::ostream &os, const T &M, std::string_view sep,
                             std::string_view begin, std::string_view end) {
    std::array<char, 64> buf;
    if (M.cols() == 1) {
        os << begin;
        for (decltype(M.rows()) r{}; r < M.rows(); ++r) {
            print_elem(buf, M(r, 0), os);
            if (r != M.rows() - 1)
                os << sep;
        }
        return os << end;
    } else {
        for (decltype(M.rows()) r{}; r < M.rows(); ++r) {
            os << begin;
            for (decltype(M.cols()) c{}; c < M.cols(); ++c) {
                print_elem(buf, M(r, c), os);
                if (c != M.cols() - 1)
                    os << sep;
            }
            os << end;
        }
        return os;
    }
}

template <class T>
std::ostream &print_matlab_impl(std::ostream &os, const T &M,
                                std::string_view end) {
    if (M.cols() == 1) {
        return print_csv_impl<T>(os, M, " ", "[", "]") << end;
    } else {
        os << '[';
        std::array<char, 64> buf;
        for (decltype(M.rows()) r{}; r < M.rows(); ++r) {
            for (decltype(M.cols()) c{}; c < M.cols(); ++c) {
                print_elem(buf, M(r, c), os);
                if (c != M.cols() - 1)
                    os << ' ';
            }
            if (r != M.rows() - 1)
                os << ";\n ";
        }
        return os << ']' << end;
    }
}

template <class T>
std::ostream &print_python_impl(std::ostream &os, const T &M,
                                std::string_view end) {
    if (M.cols() == 1) {
        return print_csv_impl<T>(os, M, ", ", "[", "]") << end;
    } else {
        os << "[[";
        std::array<char, 64> buf;
        for (decltype(M.rows()) r{}; r < M.rows(); ++r) {
            for (decltype(M.cols()) c{}; c < M.cols(); ++c) {
                print_elem(buf, M(r, c), os);
                if (c != M.cols() - 1)
                    os << ", ";
            }
            if (r != M.rows() - 1)
                os << "],\n [";
        }
        return os << "]]" << end;
    }
}

} // namespace alpaqa