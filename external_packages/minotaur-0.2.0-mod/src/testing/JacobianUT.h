// 
//     MINOTAUR -- It's only 1/2 bull
// 
//     (C)opyright 2009 - 2014 The MINOTAUR Team.
// 

#ifndef JACOBIANUT_H
#define JACOBIANUT_H

#include <string>
#include <cppunit/TestCase.h>
#include <cppunit/TestCaller.h>
#include <cppunit/TestSuite.h>
#include <cppunit/TestResult.h>
#include <cppunit/extensions/HelperMacros.h>

#include "Function.h"
#include "Jacobian.h"
#include "LinearFunction.h"
#include "Problem.h"

using namespace Minotaur;

class JacobianUT : public CppUnit::TestCase {

public:
  JacobianUT(std::string name) : TestCase(name) {}
  JacobianUT() {}

  void setUp();
  void tearDown();

  CPPUNIT_TEST_SUITE(JacobianUT);
  CPPUNIT_TEST(testLinearEval);
  CPPUNIT_TEST(testQuadEval);
  CPPUNIT_TEST_SUITE_END();

  void testLinearEval();
  void testQuadEval();

private:
  LinearFunctionPtr lf_;
  QuadraticFunctionPtr qf_;
  ProblemPtr instance_;
  std::vector<ConstraintPtr> cons_;
  std::vector<VariablePtr> vars_;
  FunctionPtr f_;
  
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
