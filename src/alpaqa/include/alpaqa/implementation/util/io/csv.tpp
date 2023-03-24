#include <alpaqa/util/io/csv.hpp>

#include <algorithm>
#include <charconv>
#include <ios>
#include <iostream>

namespace alpaqa::csv {

template <std::floating_point F>
struct CSVReader {
    static constexpr std::streamsize bufmaxsize = 64;
    std::array<char, bufmaxsize> s;
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
        }
        // Parse a number
        F v;
        const char *bufend   = s.data() + bufidx;
        const auto [ptr, ec] = std::from_chars(s.data(), bufend, v);
        const auto bufvw     = std::string_view(s.data(), bufend);
        if (ec != std::errc{})
            throw read_error("csv::read_row conversion failed '" +
                             std::string(bufvw) + "'");
        // Check separator
        if (ptr != bufend && *ptr != sep)
            throw read_error("csv::read_row unexpected character '" +
                             std::string{*ptr} + "'");
        // Shift the buffer over
        if (ptr != bufend) {
            std::copy(ptr + 1, bufend, s.data());
            bufidx -= ptr + 1 - s.data();
        } else {
            bufidx = 0;
        }
        return v;
    }

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