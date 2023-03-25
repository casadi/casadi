#pragma once

#include <alpaqa/config/config.hpp>
#include <alpaqa/util/print.hpp>
#include <concepts>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <span>
#include <string_view>

namespace fs = std::filesystem;

// Fundamental types
template <std::floating_point F>
void python_literal(std::ostream &os, F v) {
    os << alpaqa::float_to_str(v);
}
template <std::integral I>
void python_literal(std::ostream &os, I v) {
    os << v;
}
void python_literal(std::ostream &os, bool v) { os << (v ? "True" : "False"); }
// Array of doubles
void python_literal(std::ostream &os,
                    const alpaqa::vec<alpaqa::DefaultConfig> &v) {
    alpaqa::print_python(os << "np.array(", v, ")");
}
void python_literal(std::ostream &os,
                    const alpaqa::mat<alpaqa::DefaultConfig> &v) {
    alpaqa::print_python(os << "np.array(", v, ")");
}
void python_literal(std::ostream &os, std::span<const double> s) {
    USING_ALPAQA_CONFIG(alpaqa::DefaultConfig);
    vec v = cmvec{s.data(), static_cast<length_t>(s.size())};
    python_literal(os, v);
}
// Strings
void python_literal(std::ostream &os, std::string_view v) {
    os << std::quoted(v);
}
void python_literal(std::ostream &os, const std::string &v) {
    python_literal(os, std::string_view{v});
}
void python_literal(std::ostream &os, const char *v) {
    python_literal(os, std::string_view{v});
}
void python_literal(std::ostream &os, const fs::path &v) { os << v; }
// List of strings
void python_literal(std::ostream &os, std::span<const std::string_view> v) {
    os << "[";
    for (const auto &el : v)
        os << std::quoted(el) << ", ";
    os << "]";
}
