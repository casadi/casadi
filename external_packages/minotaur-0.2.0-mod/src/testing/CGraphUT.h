// 
//     MINOTAUR -- It's only 1/2 bull
// 
//     (C)opyright 2009 - 2014 The MINOTAUR Team.
// 

#ifndef CGRAPHUT_H
#define CGRAPHUT_H

#include <cppunit/TestCase.h>
#include <cppunit/TestCaller.h>
#include <cppunit/TestSuite.h>
#include <cppunit/TestResult.h>
#include <cppunit/extensions/HelperMacros.h>

#include <Types.h>
using namespace Minotaur;

class CGraphUT : public CppUnit::TestCase {

public:
  CGraphUT(std::string name) : TestCase(name) {}
  CGraphUT() {}

  void setUp() { }      // need not implement
  void tearDown() { }   // need not implement
  void testIdentical();
  void testLin();
  void testQuad();

  CPPUNIT_TEST_SUITE(CGraphUT);
  CPPUNIT_TEST(testIdentical);
  CPPUNIT_TEST(testLin);
  CPPUNIT_TEST(testQuad);
  CPPUNIT_TEST_SUITE_END();

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
