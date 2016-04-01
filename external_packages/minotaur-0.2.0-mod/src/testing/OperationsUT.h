// 
//     MINOTAUR -- It's only 1/2 bull
// 
//     (C)opyright 2009 - 2014 The MINOTAUR Team.
// 

#ifndef OPERATIONSUT_H
#define OPERATIONSUT_H

#include <vector>

#include <cppunit/TestCase.h>
#include <cppunit/TestCaller.h>
#include <cppunit/TestSuite.h>
#include <cppunit/TestResult.h>
#include <cppunit/extensions/HelperMacros.h>

#include "Types.h"

using namespace Minotaur;

class OperationsTest : public CppUnit::TestCase {

public:
  OperationsTest(std::string name) : TestCase(name) {}
  OperationsTest() {}

  CPPUNIT_TEST_SUITE(OperationsTest);
  CPPUNIT_TEST(testGcd);
  CPPUNIT_TEST(testToLower);
  CPPUNIT_TEST(testSortVarX);
  CPPUNIT_TEST_SUITE_END();

  void testGcd();
  void testToLower();
  void testSortVarX();
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
