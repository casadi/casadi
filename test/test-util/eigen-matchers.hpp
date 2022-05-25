#pragma once

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <alpaqa/util/print.hpp>

/// @file
/// @see https://google.github.io/googletest/reference/matchers.html#defining-matchers

MATCHER_P(EigenEqual, expect, "") {
    auto diffnorm = (arg - expect).template lpNorm<Eigen::Infinity>();
    if (auto *os = result_listener->stream()) {
        *os << "\nactual = ...\n";
        alpaqa::print_matlab(*os, arg.eval());
        *os << "and expected = ...\n";
        alpaqa::print_matlab(*os, expect.eval());
        *os << "with difference = ...\n";
        alpaqa::print_matlab(*os, (arg - expect).eval());
        *os << "which has infinity norm                       = " << diffnorm;
    }
    return diffnorm == 0;
}

MATCHER_P2(EigenAlmostEqual, expect, atol, "") {
    auto diffnorm = (arg - expect).template lpNorm<Eigen::Infinity>();
    if (auto *os = result_listener->stream()) {
        *os << "\nactual = ...\n";
        alpaqa::print_matlab(*os, arg.eval());
        *os << "and expected = ...\n";
        alpaqa::print_matlab(*os, expect.eval());
        *os << "with difference = ...\n";
        alpaqa::print_matlab(*os, (arg - expect).eval());
        *os << "which has infinity norm                      = " << diffnorm;
        *os << ",\nwhich is greater than the absolute tolerance = " << atol;
    }
    return diffnorm <= atol;
}

MATCHER_P2(EigenAlmostEqualRel, expect, rtol, "") {
    auto diffnorm =
        (arg - expect).cwiseQuotient(expect).template lpNorm<Eigen::Infinity>();
    if (auto *os = result_listener->stream()) {
        *os << "\nactual = ...\n";
        alpaqa::print_matlab(*os, arg.eval());
        *os << "and expected = ...\n";
        alpaqa::print_matlab(*os, expect.eval());
        *os << "with difference = ...\n";
        alpaqa::print_matlab(*os, (arg - expect).eval());
        *os << "which has relative infinity norm              = " << diffnorm;
        *os << ",\nwhich is greater than the relative tolerance = " << rtol;
    }
    return diffnorm <= rtol;
}
