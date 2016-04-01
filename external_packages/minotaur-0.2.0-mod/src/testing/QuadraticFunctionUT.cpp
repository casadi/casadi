// 
//    MINOTAUR -- It's only 1/2 bull
// 
//    (C)opyright 2009 - 2014 The MINOTAUR Team.
// 


#include <cmath>

#include "MinotaurConfig.h"
#include "LinearFunction.h"
#include "QuadraticFunctionUT.h"
#include "Variable.h"

CPPUNIT_TEST_SUITE_REGISTRATION(QuadraticFunctionTest);
CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(QuadraticFunctionTest, "QuadraticFunctionUT");

using namespace Minotaur;


void QuadraticFunctionTest::setUp()
{
  std::string vname = "common_var_name";
  VariablePtr v;
  v = (VariablePtr) new Variable(0, 0, 0.0, 1.0, Integer, vname);
  vars_.push_back(v);

  v = (VariablePtr) new Variable(1, 1, 9.0, 1.0, Integer, vname);
  vars_.push_back(v);

  v = (VariablePtr) new Variable(2, 2, -1.0, 1.0, Integer, vname);
  vars_.push_back(v);

  v = (VariablePtr) new Variable(3, 3, -100.0, 100.0, Integer, vname);
  vars_.push_back(v);

  q_ = (QuadraticFunctionPtr) new QuadraticFunction();
  q1_ = (QuadraticFunctionPtr) new QuadraticFunction();

  // x0^2 + 2x0.x1 + x1^2
  q_->addTerm(vars_[0], vars_[1], 2.0);
  q_->addTerm(vars_[0], vars_[0], 1.0);
  q_->addTerm(vars_[1], vars_[1], 1.0);

  // -2x0.x1 - 2x1^2
  q1_->addTerm(vars_[0], vars_[1], -2.0);
  q1_->addTerm(vars_[1], vars_[1], -2.0);
}


void QuadraticFunctionTest::tearDown()
{
  vars_.clear();
}


void QuadraticFunctionTest::testGetCoeffs()
{
  VarCountConstMap * vmap = q_->getVarMap();
  CPPUNIT_ASSERT(q_->getNumTerms() == 3);
  CPPUNIT_ASSERT(q_->getNumVars() == 2);
  CPPUNIT_ASSERT(q_->getWeight(vars_[0], vars_[0]) == 1.0);
  CPPUNIT_ASSERT(q_->getWeight(vars_[0], vars_[1]) == 2.0);
  CPPUNIT_ASSERT(q_->getWeight(vars_[1], vars_[0]) == 2.0);
  CPPUNIT_ASSERT(q_->getWeight(vars_[1], vars_[1]) == 1.0);

  CPPUNIT_ASSERT((vmap->find(vars_[0]))->second==3);
  CPPUNIT_ASSERT((vmap->find(vars_[1]))->second==3);

  vmap = q1_->getVarMap();
  CPPUNIT_ASSERT((vmap->find(vars_[0]))->second==1);
  CPPUNIT_ASSERT((vmap->find(vars_[1]))->second==3);
}


void QuadraticFunctionTest::testEvaluate()
{
  std::vector<double> x(2,0);
  std::vector<double> dq(2,0);
  CPPUNIT_ASSERT(q_->eval(x) == 0);
  CPPUNIT_ASSERT(q1_->eval(x) == 0);

  q_->evalGradient(x, dq);
  CPPUNIT_ASSERT(dq[0] == 0);
  CPPUNIT_ASSERT(dq[1] == 0);

  x[0] = 12.0;
  CPPUNIT_ASSERT(q_->eval(x) == 144);
  CPPUNIT_ASSERT(q1_->eval(x) == 0);

  q_->evalGradient(x, dq);
  CPPUNIT_ASSERT(dq[0] == 24);
  CPPUNIT_ASSERT(dq[1] == 24);

  dq = std::vector<double>(2,0);
  q1_->evalGradient(x, dq);
  CPPUNIT_ASSERT(dq[0] == 0);
  CPPUNIT_ASSERT(dq[1] == -24);

  x[1] = -1.0;
  CPPUNIT_ASSERT(q_->eval(x) == 121);
  CPPUNIT_ASSERT(q1_->eval(x) == 22);

  dq = std::vector<double>(2,0);
  q_->evalGradient(x, dq);
  CPPUNIT_ASSERT(dq[0] == 22);
  CPPUNIT_ASSERT(dq[1] == 22);

  dq = std::vector<double>(2,0);
  q1_->evalGradient(x, dq);
  CPPUNIT_ASSERT(dq[0] == 2);
  CPPUNIT_ASSERT(dq[1] == -20);
}


void QuadraticFunctionTest::testOperations()
{
  q2_ = q_ + q_;
  CPPUNIT_ASSERT(q2_->getNumTerms() == 3);
  CPPUNIT_ASSERT(q2_->getNumVars() == 2);
  CPPUNIT_ASSERT(q2_->getWeight(vars_[0], vars_[0]) == 2.0);
  CPPUNIT_ASSERT(q2_->getWeight(vars_[0], vars_[1]) == 4.0);
  CPPUNIT_ASSERT(q2_->getWeight(vars_[1], vars_[0]) == 4.0);
  CPPUNIT_ASSERT(q2_->getWeight(vars_[1], vars_[1]) == 2.0);

  q2_ = -6.0*q_;
  CPPUNIT_ASSERT(q2_->getNumTerms() == 3);
  CPPUNIT_ASSERT(q2_->getNumVars() == 2);
  CPPUNIT_ASSERT(q2_->getWeight(vars_[0], vars_[0]) == -6.0);
  CPPUNIT_ASSERT(q2_->getWeight(vars_[0], vars_[1]) == -12.0);
  CPPUNIT_ASSERT(q2_->getWeight(vars_[1], vars_[0]) == -12.0);
  CPPUNIT_ASSERT(q2_->getWeight(vars_[1], vars_[1]) == -6.0);

  q2_ = q1_ + q2_;
  CPPUNIT_ASSERT(q2_->getNumTerms() == 3);
  CPPUNIT_ASSERT(q2_->getNumVars() == 2);
  CPPUNIT_ASSERT(q2_->getWeight(vars_[0], vars_[0]) == -6.0);
  CPPUNIT_ASSERT(q2_->getWeight(vars_[0], vars_[1]) == -14.0);
  CPPUNIT_ASSERT(q2_->getWeight(vars_[1], vars_[0]) == -14.0);
  CPPUNIT_ASSERT(q2_->getWeight(vars_[1], vars_[1]) == -8.0);

  q2_ = q1_ - q2_;
  CPPUNIT_ASSERT(q2_->getNumTerms() == 3);
  CPPUNIT_ASSERT(q2_->getNumVars() == 2);
  CPPUNIT_ASSERT(q2_->getWeight(vars_[0], vars_[0]) ==  6.0);
  CPPUNIT_ASSERT(q2_->getWeight(vars_[0], vars_[1]) ==  12.0);
  CPPUNIT_ASSERT(q2_->getWeight(vars_[1], vars_[0]) ==  12.0);
  CPPUNIT_ASSERT(q2_->getWeight(vars_[1], vars_[1]) ==  6.0);

  (*q2_) += q1_;
  CPPUNIT_ASSERT(q2_->getNumTerms() == 3);
  CPPUNIT_ASSERT(q2_->getNumVars() == 2);
  CPPUNIT_ASSERT(q2_->getWeight(vars_[0], vars_[0]) ==  6.0);
  CPPUNIT_ASSERT(q2_->getWeight(vars_[0], vars_[1]) ==  10.0);
  CPPUNIT_ASSERT(q2_->getWeight(vars_[1], vars_[0]) ==  10.0);
  CPPUNIT_ASSERT(q2_->getWeight(vars_[1], vars_[1]) ==  4.0);

  q2_ = q_->clone();
  CPPUNIT_ASSERT(q2_->getNumTerms() == 3);
  CPPUNIT_ASSERT(q2_->getNumVars() == 2);
  CPPUNIT_ASSERT(q2_->getWeight(vars_[0], vars_[0]) == 1.0);
  CPPUNIT_ASSERT(q2_->getWeight(vars_[0], vars_[1]) == 2.0);
  CPPUNIT_ASSERT(q2_->getWeight(vars_[1], vars_[0]) == 2.0);
  CPPUNIT_ASSERT(q2_->getWeight(vars_[1], vars_[1]) == 1.0);

  LinearFunctionPtr l1 = (LinearFunctionPtr) new LinearFunction();
  LinearFunctionPtr l2;

  l1->addTerm(vars_[0], 2);
  q2_ = l1 * l2;
  // q2_ should be NULL
  CPPUNIT_ASSERT(!q2_);

  l2 = (LinearFunctionPtr) new LinearFunction();
  l2->addTerm(vars_[0], -1);
  l2->addTerm(vars_[1],  2);
  q2_ = l1 * l2;
  CPPUNIT_ASSERT(q2_->getNumTerms() == 2);
  CPPUNIT_ASSERT(q2_->getNumVars() == 2);
  CPPUNIT_ASSERT(q2_->getWeight(vars_[0], vars_[0]) == -2.0);
  CPPUNIT_ASSERT(q2_->getWeight(vars_[0], vars_[1]) == 4.0);
  CPPUNIT_ASSERT(q2_->getWeight(vars_[1], vars_[0]) == 4.0);
}


void QuadraticFunctionTest::testEigen()
{
  EigenCalculator *e_cal = new EigenCalculator();
  EigenPtr ePtr = e_cal->findValues(q_);

  CPPUNIT_ASSERT(ePtr->numNegative() == 0);
  CPPUNIT_ASSERT(ePtr->numZero() == 1);
  CPPUNIT_ASSERT(ePtr->numPositive() == 1);

  EigenPairConstIterator it=ePtr->begin();
  CPPUNIT_ASSERT(it->first == 0.0);
  ++it;
  CPPUNIT_ASSERT(it->first == 2.0);

  ePtr = e_cal->findValues(q1_);
  CPPUNIT_ASSERT(ePtr->numNegative() == 1);
  CPPUNIT_ASSERT(ePtr->numZero() == 0);
  CPPUNIT_ASSERT(ePtr->numPositive() == 1);

  it=ePtr->begin();
  CPPUNIT_ASSERT(fabs(it->first + 2.414214) < 1e-5);
  ++it;
  CPPUNIT_ASSERT(fabs(it->first - 0.414214) < 1e-5);

  delete e_cal;
}


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
