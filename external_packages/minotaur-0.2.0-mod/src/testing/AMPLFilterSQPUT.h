// 
//     MINOTAUR -- It's only 1/2 bull
// 
//     (C)opyright 2009 - 2014 The MINOTAUR Team.
// 

#ifndef AMPLFILTERSQPUT_H
#define AMPLFILTERSQPUT_H

#include <string>

#include <cppunit/TestCase.h>
#include <cppunit/TestCaller.h>
#include <cppunit/TestSuite.h>
#include <cppunit/TestResult.h>
#include <cppunit/extensions/HelperMacros.h>

#include <AMPLInterface.h>
#include <Problem.h>

using namespace MINOTAUR_AMPL;

class AMPLFilterSQPUT : public CppUnit::TestCase {

public:
  AMPLFilterSQPUT(std::string name) : TestCase(name) {}
  AMPLFilterSQPUT() {}

  void testNLP();
  void testLP();
  void setUp() { }      // need not implement
  void tearDown() { }   // need not implement

  CPPUNIT_TEST_SUITE(AMPLFilterSQPUT);
  CPPUNIT_TEST(testNLP);
  CPPUNIT_TEST(testLP);
  CPPUNIT_TEST_SUITE_END();

private:
  AMPLInterfacePtr iface_;
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
