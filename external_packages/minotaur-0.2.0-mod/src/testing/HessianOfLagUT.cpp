// 
//    MINOTAUR -- It's only 1/2 bull
// 
//    (C)opyright 2009 - 2014 The MINOTAUR Team.
// 

#include <cmath>

#include "MinotaurConfig.h"
#include "Function.h"
#include "HessianOfLagUT.h"
#include "LinearFunction.h"
#include "QuadraticFunction.h"

CPPUNIT_TEST_SUITE_REGISTRATION(HessianOfLagUT);
CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(HessianOfLagUT, "HessianOfLagUT");

using namespace Minotaur;


void HessianOfLagUT::setUp()
{
  instance_ = (ProblemPtr) new Problem();
  // 6 variables
  vars_.push_back(instance_->newVariable(0.0, INFINITY, Integer));
  vars_.push_back(instance_->newVariable(0.0, 1.0, Integer));
  vars_.push_back(instance_->newVariable(0.0, 2.0, Integer));
  vars_.push_back(instance_->newVariable(0.0, 3.0, Integer));
  vars_.push_back(instance_->newVariable(0.0, 4.0, Integer));
  vars_.push_back(instance_->newVariable(0.0, 5.0, Integer));

  // 7x_0 - 2x_1 <= 14
  lf_ = (LinearFunctionPtr) new LinearFunction();
  lf_->addTerm(vars_[0], 7.0);
  lf_->addTerm(vars_[1], -2.0);
  f_ = (FunctionPtr) new Function(lf_);
  instance_->newConstraint(f_, -INFINITY, 14.0, "cons0");

  // x_0 - x_2 <= 14
  lf_ = (LinearFunctionPtr) new LinearFunction();
  lf_->addTerm(vars_[0], 1.0);
  lf_->addTerm(vars_[2], -1.0);
  f_ = (FunctionPtr) new Function(lf_);
  instance_->newConstraint(f_, -INFINITY, 14.0, "cons1");

  f_ = FunctionPtr();
  instance_->newObjective(f_, 0.0, Minimize);
}


void HessianOfLagUT::tearDown()
{
  vars_.clear();
  lf_.reset();
  qf_.reset();
}


void HessianOfLagUT::testEmpty()
{
  HessianOfLagPtr hess;

  qf_ = (QuadraticFunctionPtr) new QuadraticFunction();
  lf_ = LinearFunctionPtr();
  f_ = (FunctionPtr) new Function(lf_, qf_);
  instance_->newConstraint(f_, -INFINITY, 10.0, "cons2");
  instance_->setNativeDer();

  // test size
  hess = instance_->getHessian();
  CPPUNIT_ASSERT(hess->getNumNz() == 0);
}


void HessianOfLagUT::testLinearEval()
{
  HessianOfLagPtr hess;
  instance_->setNativeDer();

  // test size
  hess = instance_->getHessian();
  CPPUNIT_ASSERT(hess->getNumNz() == 0);
}
 
 
void HessianOfLagUT::testQuadEval()
{
  HessianOfLagPtr hess;
  int error = 0;
  double mults[] = {1.0, 1.0, 4.0, 7.0, -1.0};

  // add quadratics to instance
  // x_2^2 + 2x_3^2 <= 10
  // x_3^2 + 8x_3 + x_4_x_5 <= 10
  // x_5^2 + x_1x_5 + x_1x_2 + x_1 <= 10
  qf_ = (QuadraticFunctionPtr) new QuadraticFunction();
  lf_ = LinearFunctionPtr();
  qf_->addTerm(vars_[2], vars_[2], 1.0);
  qf_->addTerm(vars_[3], vars_[3], 2.0);
  f_ = (FunctionPtr) new Function(lf_, qf_);
  instance_->newConstraint(f_, -INFINITY, 10.0, "cons2");

  qf_ = (QuadraticFunctionPtr) new QuadraticFunction();
  lf_ = (LinearFunctionPtr) new LinearFunction();
  qf_->addTerm(vars_[3], vars_[3], 1.0);
  lf_->addTerm(vars_[3], 8.0);
  qf_->addTerm(vars_[4], vars_[5], 1.0);
  f_ = (FunctionPtr) new Function(lf_, qf_);
  instance_->newConstraint(f_, -INFINITY, 10.0, "cons3");

  qf_ = (QuadraticFunctionPtr) new QuadraticFunction();
  lf_ = (LinearFunctionPtr) new LinearFunction();
  qf_->addTerm(vars_[5], vars_[5], 1.0);
  qf_->addTerm(vars_[1], vars_[5], 1.0);
  qf_->addTerm(vars_[1], vars_[2], 1.0);
  lf_->addTerm(vars_[1], 1.0);
  f_ = (FunctionPtr) new Function(lf_, qf_);
  instance_->newConstraint(f_, -INFINITY, 10.0, "cons4");

  instance_->setNativeDer();

  // test size
  hess = instance_->getHessian();
  CPPUNIT_ASSERT(hess->getNumNz() == 6);

  // test indices
  UInt iRow[6];
  UInt jCol[6];
  hess->fillRowColIndices(&iRow[0], &jCol[0]);
  CPPUNIT_ASSERT(iRow[0] == 2);
  CPPUNIT_ASSERT(jCol[0] == 1);

  CPPUNIT_ASSERT(iRow[1] == 2);
  CPPUNIT_ASSERT(jCol[1] == 2);

  CPPUNIT_ASSERT(iRow[2] == 3);
  CPPUNIT_ASSERT(jCol[2] == 3);

  CPPUNIT_ASSERT(iRow[3] == 5);
  CPPUNIT_ASSERT(jCol[3] == 1);

  CPPUNIT_ASSERT(iRow[4] == 5);
  CPPUNIT_ASSERT(jCol[4] == 4);

  CPPUNIT_ASSERT(iRow[5] == 5);
  CPPUNIT_ASSERT(jCol[5] == 5);

  // test values
  // x_2^2 + 2x_3^2 <= 10                    (lambda[2] = 4)
  // x_3^2 + 8x_3 + x_4_x_5 <= 10            (lambda[3] = 7)
  // x_5^2 + x_1x_5 + x_1x_2 + x_1 <= 10     (lambda[4] = -1)
  
  double values[6];
  double x[6];
  x[0] = 0.0;
  x[1] = 1.0;
  x[2] = 2.0;
  x[3] = -3.0;
  x[4] = -3.0;
  x[5] = 10.0;
  std::fill(&values[0], &values[0]+6, 0);
  hess->fillRowColValues(x, 1, mults, values, &error);
  CPPUNIT_ASSERT(0==error);
  CPPUNIT_ASSERT(values[0] == -1.0);
  CPPUNIT_ASSERT(values[1] == 8.0);
  CPPUNIT_ASSERT(values[2] == 30.0);
  CPPUNIT_ASSERT(values[3] == -1.0);
  CPPUNIT_ASSERT(values[4] == 7.0);
  CPPUNIT_ASSERT(values[5] == -2.0);
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
