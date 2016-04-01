// 
//     MINOTAUR -- It's only 1/2 bull
// 
//     (C)opyright 2009 - 2014 The MINOTAUR Team.
// 

#ifndef AMPLOSIUT_H
#define AMPLOSIUT_H

#include <cppunit/TestCase.h>
#include <cppunit/TestCaller.h>
#include <cppunit/TestSuite.h>
#include <cppunit/TestResult.h>
#include <cppunit/extensions/HelperMacros.h>

#include "Types.h"
#include "OsiLPEngine.h"
#include "AMPLInterface.h"

using namespace MINOTAUR_AMPL;

class AMPLOsiUT : public CppUnit::TestCase {

public:
  AMPLOsiUT(std::string name) : TestCase(name) {}
  AMPLOsiUT() {}

  void testOsiLP();
  void testOsiLP2();
  void testOsiWarmStart();
  void testOsiBnB();
  void setUp();
  void tearDown();

  CPPUNIT_TEST_SUITE(AMPLOsiUT);
  CPPUNIT_TEST(testOsiLP);
  CPPUNIT_TEST(testOsiLP2);
  CPPUNIT_TEST(testOsiWarmStart);
  CPPUNIT_TEST(testOsiBnB);
  CPPUNIT_TEST_SUITE_END();

private:
  AMPLInterfacePtr iface_;
  Minotaur::OsiLPEnginePtr engine_ptr_;
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
