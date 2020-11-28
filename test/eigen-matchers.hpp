#pragma once

#include <gmock/gmock-matchers.h>
#include <gmock/gmock.h>
#include <gtest/gtest.h>

/// @file
/// https://stackoverflow.com/questions/25146997/teach-google-test-how-to-print-eigen-matrix

MATCHER_P(EigenEqual, expect,
          std::string(negation ? "isn't" : "is") + " equal to" +
              testing::PrintToString(expect)) {
    return arg == expect;
}

template <class Base>
class EigenPrintWrap : public Base {
    friend void PrintTo(const EigenPrintWrap &m, std::ostream *o) {
        *o << "\n" << m;
    }
};

template <class Base>
const EigenPrintWrap<Base> &print_wrap(const Base &base) {
    return static_cast<const EigenPrintWrap<Base> &>(base);
}
