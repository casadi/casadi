#pragma once

#include <stdexcept>
#include <string>

void check_dim(const char *name, const auto &vec, auto sz) {
    if (vec.size() != sz)
        throw std::invalid_argument(std::string(name) + ": dimension mismatch");
}

void check_dim_msg(const auto &vec, auto sz, std::string msg) {
    if (vec.size() != sz)
        throw std::invalid_argument(std::move(msg));
}