// 
//    MINOTAUR -- It's only 1/2 bull
// 
//    (C)opyright 2009 - 2014 The MINOTAUR Team.
// 

#include <cmath>

#include "MinotaurConfig.h"
#include "Jacobian.h"
#include "JacobianUT.h"
#include "QuadraticFunction.h"

CPPUNIT_TEST_SUITE_REGISTRATION(JacobianUT);
CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(JacobianUT, "JacobianUT");

using namespace Minotaur;


void JacobianUT::setUp()
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

}


void JacobianUT::tearDown()
{
  vars_.clear();
  lf_.reset();
  qf_.reset();
}


void JacobianUT::testLinearEval()
{
  JacobianPtr jac;
  int error = 0;
  instance_->setNativeDer();

  // test size
  jac = instance_->getJacobian();
  CPPUNIT_ASSERT(jac->getNumNz() == 4);

  // test indices
  UInt iRow[4];
  UInt jCol[4];
  jac->fillRowColIndices(&iRow[0], &jCol[0]);
  CPPUNIT_ASSERT(iRow[0] == 0);
  CPPUNIT_ASSERT(iRow[1] == 0);
  CPPUNIT_ASSERT(iRow[2] == 1);
  CPPUNIT_ASSERT(iRow[3] == 1);
  CPPUNIT_ASSERT(jCol[0] == 0);
  CPPUNIT_ASSERT(jCol[1] == 1);
  CPPUNIT_ASSERT(jCol[2] == 0);
  CPPUNIT_ASSERT(jCol[3] == 2);

  // test values
  double values[4];
  double x[4];
  x[0] = 0.0;
  x[1] = 1.0;
  x[2] = 2.0;
  x[3] = -3.0;
  values[0] = values[1] = values[2] = values[3] = 0;
  jac->fillRowColValues(&x[0], &values[0], &error);
  CPPUNIT_ASSERT(0==error);
  CPPUNIT_ASSERT(values[0] == 7.0);
  CPPUNIT_ASSERT(values[1] == -2.0);
  CPPUNIT_ASSERT(values[2] == 1.0);
  CPPUNIT_ASSERT(values[3] == -1.0);
}
 
 
void JacobianUT::testQuadEval()
{
  JacobianPtr jac;
  int error = 0;

  // add quadratics to instance
  // x_2^2 + 2x_3^2 <= 10
  // x_3^2 + 8x_3 + x_4_x_5 <= 10
  // x_5^2 + x_5 + x_1x_2 + x_1 <= 10
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
  lf_->addTerm(vars_[5], 1.0);
  qf_->addTerm(vars_[1], vars_[2], 1.0);
  lf_->addTerm(vars_[1], 1.0);
  f_ = (FunctionPtr) new Function(lf_, qf_);
  instance_->newConstraint(f_, -INFINITY, 10.0, "cons4");

  instance_->setNativeDer();

  // test size
  jac = instance_->getJacobian();
  CPPUNIT_ASSERT(jac->getNumNz() == 12);

  // test indices
  UInt iRow[12];
  UInt jCol[12];
  jac->fillRowColIndices(&iRow[0], &jCol[0]);
  // first 4 already tested
  CPPUNIT_ASSERT(iRow[4] == 2);
  CPPUNIT_ASSERT(iRow[5] == 2);
  CPPUNIT_ASSERT(iRow[6] == 3);
  CPPUNIT_ASSERT(iRow[7] == 3);
  CPPUNIT_ASSERT(iRow[8] == 3);
  CPPUNIT_ASSERT(iRow[9] == 4);
  CPPUNIT_ASSERT(iRow[10] == 4);
  CPPUNIT_ASSERT(iRow[11] == 4);
  CPPUNIT_ASSERT(jCol[4] == 2);
  CPPUNIT_ASSERT(jCol[5] == 3);
  CPPUNIT_ASSERT(jCol[6] == 3);
  CPPUNIT_ASSERT(jCol[7] == 4);
  CPPUNIT_ASSERT(jCol[8] == 5);
  CPPUNIT_ASSERT(jCol[9] == 1);
  CPPUNIT_ASSERT(jCol[10] == 2);
  CPPUNIT_ASSERT(jCol[11] == 5);

  // test values
  // x_2^2 + 2x_3^2 <= 10
  // x_3^2 + 8x_3 + x_4_x_5 <= 10
  // x_5^2 + x_5 + x_1x_2 + x_1 <= 10
  double values[12];
  double x[6];
  x[0] = 0.0;
  x[1] = 1.0;
  x[2] = 2.0;
  x[3] = -3.0;
  x[4] = -3.0;
  x[5] = 10.0;
  std::fill(&values[0], &values[0]+12, 0);
  jac->fillRowColValues(&x[0], &values[0], &error);
  CPPUNIT_ASSERT(values[4] == 4.0);
  CPPUNIT_ASSERT(values[5] == -12.0);
  CPPUNIT_ASSERT(values[6] == 2.0);
  CPPUNIT_ASSERT(values[7] == 10.0);
  CPPUNIT_ASSERT(values[8] == -3.0);
  CPPUNIT_ASSERT(values[9] == 3.0);
  CPPUNIT_ASSERT(values[10] == 1.0);
  CPPUNIT_ASSERT(values[11] == 21.0);
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
