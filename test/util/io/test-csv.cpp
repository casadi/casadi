#include <alpaqa/config/config.hpp>
#include <alpaqa/util/io/csv.hpp>
#include <gtest/gtest.h>
#include <test-util/eigen-matchers.hpp>
#include <sstream>

TEST(csv, read) {
    USING_ALPAQA_CONFIG(alpaqa::DefaultConfig);
    std::istringstream is{"1.00000000000000000,"
                          "+2.00000000000000000,"
                          "-3.00000000000000000,"
                          "4.00000000000000000,"
                          "5.00000000000000000\n"};
    vec v(5);
    alpaqa::csv::read_row(is, rvec{v});
    vec expected(5);
    expected << 1, +2, -3, 4, 5;
    EXPECT_THAT(v, EigenEqual(expected));
}

TEST(csv, readEOF) {
    USING_ALPAQA_CONFIG(alpaqa::DefaultConfig);
    std::istringstream is{"1.00000000000000000,"
                          "2.00000000000000000,"
                          "3.00000000000000000,"
                          "4.00000000000000000,"
                          "5.00000000000000000"};
    vec v(5);
    alpaqa::csv::read_row(is, rvec{v});
    vec expected(5);
    expected << 1, 2, 3, 4, 5;
    EXPECT_THAT(v, EigenEqual(expected));
}

TEST(csv, readEndWithSep) {
    USING_ALPAQA_CONFIG(alpaqa::DefaultConfig);
    std::istringstream is{"1.00000000000000000,"
                          "2.00000000000000000,"
                          "3.00000000000000000,"
                          "4.00000000000000000,"
                          "5.00000000000000000,\n"};
    vec v(5);
    alpaqa::csv::read_row(is, rvec{v});
    vec expected(5);
    expected << 1, 2, 3, 4, 5;
    EXPECT_THAT(v, EigenEqual(expected));
}

TEST(csv, readEndWithSepEOF) {
    USING_ALPAQA_CONFIG(alpaqa::DefaultConfig);
    std::istringstream is{"1.00000000000000000,"
                          "2.00000000000000000,"
                          "3.00000000000000000,"
                          "4.00000000000000000,"
                          "5.00000000000000000,"};
    vec v(5);
    alpaqa::csv::read_row(is, rvec{v});
    vec expected(5);
    expected << 1, 2, 3, 4, 5;
    EXPECT_THAT(v, EigenEqual(expected));
}

TEST(csv, readLong) {
    USING_ALPAQA_CONFIG(alpaqa::DefaultConfig);
    std::istringstream is{"1.000000000000000000000000000000000000000000000000,"
                          "2.000000000000000000000000000000000000000000000000,"
                          "3.000000000000000000000000000000000000000000000000,"
                          "4.000000000000000000000000000000000000000000000000,"
                          "5.000000000000000000000000000000000000000000000000"};
    vec v(5);
    alpaqa::csv::read_row(is, rvec{v});
    vec expected(5);
    expected << 1, 2, 3, 4, 5;
    EXPECT_THAT(v, EigenEqual(expected));
}

TEST(csv, readShort) {
    USING_ALPAQA_CONFIG(alpaqa::DefaultConfig);
    std::istringstream is{"1.0,"
                          "2.0,"
                          "3.0,"
                          "4.0,"
                          "5.0"};
    vec v(5);
    alpaqa::csv::read_row(is, rvec{v});
    vec expected(5);
    expected << 1, 2, 3, 4, 5;
    EXPECT_THAT(v, EigenEqual(expected));
}

TEST(csv, readNaNInf) {
    USING_ALPAQA_CONFIG(alpaqa::DefaultConfig);
    std::istringstream is{"inf,"
                          "+inf,"
                          "-inf,"
                          "nan,"
                          "-nan"};
    vec v(5);
    alpaqa::csv::read_row(is, rvec{v});
    vec expected(5);
    expected << alpaqa::inf<config_t>, +alpaqa::inf<config_t>,
        -alpaqa::inf<config_t>, alpaqa::NaN<config_t>, -alpaqa::NaN<config_t>;
    EXPECT_EQ(std::memcmp(v.data(), expected.data(), 5 * sizeof(*v.data())), 0);
}

TEST(csv, readEmpty) {
    USING_ALPAQA_CONFIG(alpaqa::DefaultConfig);
    std::istringstream is{"\n"};
    vec v(0);
    alpaqa::csv::read_row(is, rvec{v});
    vec expected(0);
    EXPECT_THAT(v, EigenEqual(expected));
}

TEST(csv, readEmptyEOF) {
    USING_ALPAQA_CONFIG(alpaqa::DefaultConfig);
    std::istringstream is{""};
    vec v(0);
    alpaqa::csv::read_row(is, rvec{v});
    vec expected(0);
    EXPECT_THAT(v, EigenEqual(expected));
}

TEST(csv, readTooMany) {
    USING_ALPAQA_CONFIG(alpaqa::DefaultConfig);
    std::istringstream is{"1.0,"
                          "2.0,"
                          "3.0,"
                          "4.0,"
                          "5.0,"
                          "6.0"};
    vec v(5);
    EXPECT_THROW(alpaqa::csv::read_row(is, rvec{v}), alpaqa::csv::read_error);
}

TEST(csv, readInvalidCharPrefix) {
    USING_ALPAQA_CONFIG(alpaqa::DefaultConfig);
    std::istringstream is{"1.0,"
                          "2.0,"
                          "k3.0,"
                          "4.0,"
                          "5.0"};
    vec v(5);
    EXPECT_THROW(alpaqa::csv::read_row(is, rvec{v}), alpaqa::csv::read_error);
}

TEST(csv, readInvalidCharInfix) {
    USING_ALPAQA_CONFIG(alpaqa::DefaultConfig);
    std::istringstream is{"1.0,"
                          "2.0,"
                          "3k.0,"
                          "4.0,"
                          "5.0"};
    vec v(5);
    EXPECT_THROW(alpaqa::csv::read_row(is, rvec{v}), alpaqa::csv::read_error);
}

TEST(csv, readInvalidCharSuffix) {
    USING_ALPAQA_CONFIG(alpaqa::DefaultConfig);
    std::istringstream is{"1.0,"
                          "2.0,"
                          "3.0k,"
                          "4.0,"
                          "5.0"};
    vec v(5);
    EXPECT_THROW(alpaqa::csv::read_row(is, rvec{v}), alpaqa::csv::read_error);
}

TEST(csv, readTooFew) {
    USING_ALPAQA_CONFIG(alpaqa::DefaultConfig);
    std::istringstream is{"1.0,"
                          "2.0,"
                          "3.0,"
                          "4.0"};
    vec v(5);
    EXPECT_THROW(alpaqa::csv::read_row(is, rvec{v}), alpaqa::csv::read_error);
}

// ---

TEST(csv, stdvecRead) {
    USING_ALPAQA_CONFIG(alpaqa::DefaultConfig);
    std::istringstream is{"1.00000000000000000,"
                          "2.00000000000000000,"
                          "3.00000000000000000,"
                          "4.00000000000000000,"
                          "5.00000000000000000\n"};
    auto v = alpaqa::csv::read_row_std_vector<real_t>(is);
    std::vector<real_t> expected{1, 2, 3, 4, 5};
    EXPECT_EQ(v, expected);
}

TEST(csv, stdvecReadEOF) {
    USING_ALPAQA_CONFIG(alpaqa::DefaultConfig);
    std::istringstream is{"1.00000000000000000,"
                          "2.00000000000000000,"
                          "3.00000000000000000,"
                          "4.00000000000000000,"
                          "5.00000000000000000"};
    auto v = alpaqa::csv::read_row_std_vector<real_t>(is);
    std::vector<real_t> expected{1, 2, 3, 4, 5};
    EXPECT_EQ(v, expected);
}

TEST(csv, stdvecReadEndWithSep) {
    USING_ALPAQA_CONFIG(alpaqa::DefaultConfig);
    std::istringstream is{"1.00000000000000000,"
                          "2.00000000000000000,"
                          "3.00000000000000000,"
                          "4.00000000000000000,"
                          "5.00000000000000000,\n"};
    auto v = alpaqa::csv::read_row_std_vector<real_t>(is);
    std::vector<real_t> expected{1, 2, 3, 4, 5};
    EXPECT_EQ(v, expected);
}

TEST(csv, stdvecReadEndWithSepEOF) {
    USING_ALPAQA_CONFIG(alpaqa::DefaultConfig);
    std::istringstream is{"1.00000000000000000,"
                          "2.00000000000000000,"
                          "3.00000000000000000,"
                          "4.00000000000000000,"
                          "5.00000000000000000,"};
    auto v = alpaqa::csv::read_row_std_vector<real_t>(is);
    std::vector<real_t> expected{1, 2, 3, 4, 5};
    EXPECT_EQ(v, expected);
}

TEST(csv, stdvecReadLong) {
    USING_ALPAQA_CONFIG(alpaqa::DefaultConfig);
    std::istringstream is{"1.000000000000000000000000000000000000000000000000,"
                          "2.000000000000000000000000000000000000000000000000,"
                          "3.000000000000000000000000000000000000000000000000,"
                          "4.000000000000000000000000000000000000000000000000,"
                          "5.000000000000000000000000000000000000000000000000"};
    auto v = alpaqa::csv::read_row_std_vector<real_t>(is);
    std::vector<real_t> expected{1, 2, 3, 4, 5};
    EXPECT_EQ(v, expected);
}

TEST(csv, stdvecReadShort) {
    USING_ALPAQA_CONFIG(alpaqa::DefaultConfig);
    std::istringstream is{"1.0,"
                          "2.0,"
                          "3.0,"
                          "4.0,"
                          "5.0"};
    auto v = alpaqa::csv::read_row_std_vector<real_t>(is);
    std::vector<real_t> expected{1, 2, 3, 4, 5};
    EXPECT_EQ(v, expected);
}

TEST(csv, stdvecReadNaNInf) {
    USING_ALPAQA_CONFIG(alpaqa::DefaultConfig);
    std::istringstream is{"inf,"
                          "-inf,"
                          "nan,"
                          "-nan"};
    auto v = alpaqa::csv::read_row_std_vector<real_t>(is);
    std::vector<real_t> expected{
        alpaqa::inf<config_t>,
        -alpaqa::inf<config_t>,
        alpaqa::NaN<config_t>,
        -alpaqa::NaN<config_t>,
    };
    EXPECT_EQ(std::memcmp(v.data(), expected.data(), 4 * sizeof(*v.data())), 0);
}

TEST(csv, stdvecReadEmpty) {
    USING_ALPAQA_CONFIG(alpaqa::DefaultConfig);
    std::istringstream is{"\n"};
    auto v = alpaqa::csv::read_row_std_vector<real_t>(is);
    std::vector<real_t> expected;
    EXPECT_EQ(v, expected);
}

TEST(csv, stdvecReadEmptyEOF) {
    USING_ALPAQA_CONFIG(alpaqa::DefaultConfig);
    std::istringstream is{""};
    auto v = alpaqa::csv::read_row_std_vector<real_t>(is);
    std::vector<real_t> expected;
    EXPECT_EQ(v, expected);
}

TEST(csv, stdvecReadInvalidCharPrefix) {
    USING_ALPAQA_CONFIG(alpaqa::DefaultConfig);
    std::istringstream is{"1.0,"
                          "2.0,"
                          "k3.0,"
                          "4.0,"
                          "5.0"};
    EXPECT_THROW(alpaqa::csv::read_row_std_vector<real_t>(is),
                 alpaqa::csv::read_error);
}

TEST(csv, stdvecReadInvalidCharInfix) {
    USING_ALPAQA_CONFIG(alpaqa::DefaultConfig);
    std::istringstream is{"1.0,"
                          "2.0,"
                          "3k.0,"
                          "4.0,"
                          "5.0"};
    EXPECT_THROW(alpaqa::csv::read_row_std_vector<real_t>(is),
                 alpaqa::csv::read_error);
}

TEST(csv, stdvecReadInvalidCharSuffix) {
    USING_ALPAQA_CONFIG(alpaqa::DefaultConfig);
    std::istringstream is{"1.0,"
                          "2.0,"
                          "3.0k,"
                          "4.0,"
                          "5.0"};
    EXPECT_THROW(alpaqa::csv::read_row_std_vector<real_t>(is),
                 alpaqa::csv::read_error);
}
