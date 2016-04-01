/* 
 *     MINOTAUR -- It's only 1/2 bull
 *
 *     (C)opyright 2009 - 2014 The MINOTAUR Team.
 */

#ifndef AMPLINSTANCEUT_H
#define AMPLINSTANCEUT_H

#include <string>

#include <cppunit/TestCase.h>
#include <cppunit/TestCaller.h>
#include <cppunit/TestSuite.h>
#include <cppunit/TestResult.h>
#include <cppunit/extensions/HelperMacros.h>

#include <AMPLInterface.h>
#include <Problem.h>

using namespace MINOTAUR_AMPL;


// read instance using ampl and test:
// Number of variables and their types,
// Number of constraints and their types,
// Function evaluations,
// Gradient evaluations,
// Hessian evaluations.

class AMPLInstanceUT : public CppUnit::TestCase {

public:
  AMPLInstanceUT(std::string name) : TestCase(name) {}
  AMPLInstanceUT() {}

  void setUp();         
  void tearDown();   // need not implement

  void testSize();
  void testVariables();
  void testConstraints();
  void testObjective();
  void testJacobian();
  void testObjectiveGradient();
  

  CPPUNIT_TEST_SUITE(AMPLInstanceUT);
  CPPUNIT_TEST(testSize);
  CPPUNIT_TEST(testVariables);
  CPPUNIT_TEST(testConstraints);
  CPPUNIT_TEST(testJacobian);
  CPPUNIT_TEST(testObjective);
  CPPUNIT_TEST(testObjectiveGradient);

  CPPUNIT_TEST_SUITE_END();

private:
  AMPLInterfacePtr iface_;
  Minotaur::ProblemPtr inst_;
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
