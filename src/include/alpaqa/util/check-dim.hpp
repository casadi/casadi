#pragma once

#include <alpaqa/config/config.hpp>
#include <stdexcept>
#include <string>

namespace alpaqa::util {

template <Config Conf>
void check_dim_msg(crvec<Conf> v, auto sz, std::string msg) {
    if (v.size() != sz) {
        msg += "\n(should be ";
        msg += std::to_string(sz);
        msg += ", got ";
        msg += std::to_string(v.size());
        msg += ")";
        throw std::invalid_argument(msg);
    }
}

template <Config Conf>
void check_dim(std::string name, crvec<Conf> v, auto sz) {
    name += ": dimension mismatch";
    check_dim_msg<Conf>(v, sz, name);
}

template <Config Conf>
void check_dim_msg(crmat<Conf> m, auto rows, auto cols, std::string msg) {
    if (m.cols() != cols || m.rows() != rows) {
        msg += "\n(should be ";
        msg += std::to_string(rows);
        msg += "×";
        msg += std::to_string(cols);
        msg += ", got ";
        msg += std::to_string(m.rows());
        msg += "×";
        msg += std::to_string(m.cols());
        msg += ")";
        throw std::invalid_argument(msg);
    }
}

template <Config Conf>
void check_dim(std::string name, crmat<Conf> m, auto rows, auto cols) {
    name += ": dimension mismatch";
    check_dim_msg<Conf>(m, rows, cols, name);
}

} // namespace alpaqa::util