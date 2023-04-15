#pragma once

#include <iterator>
#include <numeric>
#include <string>

std::string format_string_list(
    const auto &container,
    const auto &proj = [](const auto &x) -> decltype(auto) { return x; }) {
    if (container.empty())
        return std::string{};
    auto penult       = std::prev(container.end());
    auto quote_concat = [&](std::string &&a, const auto &b) {
        return a + "'" + std::string(proj(b)) + "', ";
    };
    return std::accumulate(container.begin(), penult, std::string{},
                           quote_concat) +
           "'" + std::string(proj(*penult)) + "'";
}