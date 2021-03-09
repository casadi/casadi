#pragma once

#include <Eigen/Eigen>
#include <gmock/gmock-matchers.h>
#include <gmock/gmock.h>
#include <gtest/gtest.h>
#include <iomanip>

/// @file
/// https://stackoverflow.com/questions/25146997/teach-google-test-how-to-print-eigen-matrix

MATCHER_P(EigenEqual, expect, testing::PrintToString(expect)) {
    return arg == expect;
}

MATCHER_P2(EigenAlmostEqual, expect, atol, testing::PrintToString(expect)) {
    return (arg - expect).template lpNorm<Eigen::Infinity>() <= atol;
}

template <class Base>
class EigenPrintWrap : public Base {
    friend void PrintTo(const EigenPrintWrap &m, std::ostream *o) {
        if (m.cols() == 1)
            *o << std::scientific << std::setprecision(17) << m.transpose();
        else
            *o << "\n" << std::scientific << std::setprecision(17) << m;
    }
};

template <class Base>
const EigenPrintWrap<Base> &print_wrap(const Base &base) {
    return static_cast<const EigenPrintWrap<Base> &>(base);
}
