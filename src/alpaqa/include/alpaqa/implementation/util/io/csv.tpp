#include <alpaqa/util/io/csv.hpp>

#include <algorithm>
#include <cassert>
#include <charconv>
#include <ios>
#include <iostream>

#if !__cpp_lib_to_chars
#include <cerrno>
#include <cstdlib> // strtod
#endif

namespace alpaqa::csv {

template <std::floating_point F>
struct CSVReader {
    static constexpr std::streamsize bufmaxsize = 64;
    std::array<char, bufmaxsize + 1> s;
    std::streamsize bufidx    = 0;
    bool keep_reading         = true;
    static constexpr char end = '\n';

    [[nodiscard]] F read(std::istream &is, char sep) {
        // Get some characters to process
        if (keep_reading) {
            if (!is.get(s.data() + bufidx, bufmaxsize - bufidx, end))
                throw read_error("csv::read_row extraction failed: " +
                                 std::to_string(is.bad()) + " " +
                                 std::to_string(is.fail()) + " " +
                                 std::to_string(is.eof()));
            bufidx += is.gcount();
            keep_reading = is.peek() != end && !is.eof();
            assert(bufidx < bufmaxsize);
        }
        // Parse a number
        F v;
        char *bufend    = s.data() + bufidx;
        const char *ptr = read_single(s.data(), bufend, v);
        // Check separator
        if (ptr != bufend && *ptr != sep)
            throw read_error("csv::read_row unexpected character '" +
                             std::string{*ptr} + "'");
        // Shift the buffer over
        if (ptr != bufend) {
            std::copy(ptr + 1, static_cast<const char *>(bufend), s.data());
            bufidx -= ptr + 1 - s.data();
        } else {
            bufidx = 0;
        }
        return v;
    }

#if __cpp_lib_to_chars
    static const char *read_single(const char *bufbegin, const char *bufend,
                                   F &v) {
        if (bufbegin != bufend && *bufbegin == '+')
            ++bufbegin;
        const auto [ptr, ec] = std::from_chars(bufbegin, bufend, v);
        const auto bufvw     = std::string_view(bufbegin, bufend);
        if (ec != std::errc{})
            throw read_error("csv::read_row conversion failed '" +
                             std::string(bufvw) +
                             "': " + std::make_error_code(ec).message());
        return ptr;
    }
#else
    static void strtod_ovl(const char *str, char **str_end, float &v) {
        v = std::strtof(str, str_end);
    }
    static void strtod_ovl(const char *str, char **str_end, double &v) {
        v = std::strtod(str, str_end);
    }
    static void strtod_ovl(const char *str, char **str_end, long double &v) {
        v = std::strtold(str, str_end);
    }
    static const char *read_single(const char *bufbegin, char *bufend, F &v) {
        *bufend = '\0';
        char *ptr;
        errno = 0;
        strtod_ovl(bufbegin, &ptr, v);
        if (errno || ptr == bufbegin)
            throw read_error("csv::read_row conversion failed '" +
                             std::string(bufbegin) +
                             "': " + std::to_string(errno));
        return ptr;
    }
#endif

    void check_end(std::istream &is) const {
        if (bufidx > 0 || (is.get() != end && is))
            throw read_error("csv::read_row line not fully consumed");
    }

    [[nodiscard]] bool done(std::istream &is) const {
        bool keep_reading = is.peek() != end && !is.eof();
        return bufidx == 0 && !keep_reading;
    }
};

template <std::floating_point F>
void read_row_impl(std::istream &is, Eigen::Ref<Eigen::VectorX<F>> v,
                   char sep) {
    CSVReader<F> reader;
    for (auto &vv : v)
        vv = reader.read(is, sep);
    reader.check_end(is);
}

template <std::floating_point F>
std::vector<F> read_row_std_vector(std::istream &is, char sep) {
    CSVReader<F> reader;
    std::vector<F> v;
    while (!reader.done(is))
        v.push_back(reader.read(is, sep));
    reader.check_end(is);
    return v;
}

} // namespace alpaqa::csv