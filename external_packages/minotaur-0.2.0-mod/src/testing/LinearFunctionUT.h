// 
//     MINOTAUR -- It's only 1/2 bull
// 
//     (C)opyright 2009 - 2014 The MINOTAUR Team.
// 

#ifndef LINEARFUNCTIONUT_H
#define LINEARFUNCTIONUT_H

#include <string>

#include <cppunit/TestCase.h>
#include <cppunit/TestCaller.h>
#include <cppunit/TestSuite.h>
#include <cppunit/TestResult.h>
#include <cppunit/extensions/HelperMacros.h>

#include <Problem.h>
#include <LinearFunction.h>

using namespace Minotaur;

class LinearFunctionTest : public CppUnit::TestCase {

public:
  LinearFunctionTest(std::string name) : TestCase(name) {}
  LinearFunctionTest() {}

  void setUp();
  void tearDown();

  CPPUNIT_TEST_SUITE(LinearFunctionTest);
  CPPUNIT_TEST(testGetCoeffs);
  CPPUNIT_TEST(testGetObj);
  CPPUNIT_TEST(testOperations);
  CPPUNIT_TEST(testFix);
  CPPUNIT_TEST_SUITE_END();

  void testGetCoeffs();
  void testGetObj();
  void testOperations();
  void testFix();

private:
  ProblemPtr instance_;
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
