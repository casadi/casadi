// 
//     MINOTAUR -- It's only 1/2 bull
// 
//     (C)opyright 2009 - 2014 The MINOTAUR Team.
// 

#ifndef LAPACKUT_H
#define LAPACKUT_H

#include <vector>

#include <cppunit/TestCase.h>
#include <cppunit/TestCaller.h>
#include <cppunit/TestSuite.h>
#include <cppunit/TestResult.h>
#include <cppunit/extensions/HelperMacros.h>

#include "Types.h"

using namespace Minotaur;

class LapackTest : public CppUnit::TestCase {

  public:
    LapackTest(std::string name) : TestCase(name) {}
    LapackTest() {}

    void setUp() {};
    void tearDown() {};

    CPPUNIT_TEST_SUITE(LapackTest);
    CPPUNIT_TEST(testWrapper);
    CPPUNIT_TEST(testEigenValues);
    CPPUNIT_TEST(testEigenValues2);
    CPPUNIT_TEST_SUITE_END();

    void testWrapper();
    void testEigenValues();
    void testEigenValues2();

  private:

};

#endif

// Local Variables: 
// mode: c++ 
// eval: (c-set-style "k&r") 
// eval: (c-set-offset 'innamespace 0) 
// eval: (setq c-basic-offset 2) 
// eval: (setq fill-column 78) 
// eval: (auto-fill-mode 1) 
// eval: (setq column-number-mode 1) 
// eval: (setq indent-tabs-mode nil) 
// End:
