// 
//     MINOTAUR -- It's only 1/2 bull
// 
//     (C)opyright 2009 - 2014 The MINOTAUR Team.
// 



#include <cmath>
#include <fstream>
#include <iostream>

#include "MinotaurConfig.h"
#include "FunctionUT.h"
#include "LinearFunction.h"
#include "QuadraticFunction.h"
#include "Variable.h"

CPPUNIT_TEST_SUITE_REGISTRATION(FunctionTest);
CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(FunctionTest, "FunctionUT");

using namespace boost;
using namespace Minotaur;
using namespace std;


void FunctionTest::setUp()
{

  std::string vname = "common_var_name";
  for(UInt i = 0; i < 10; ++i) {
    // Create some variables
    vars_.push_back(VariablePtr(new Variable(i, i, -INFINITY, INFINITY, 
            Continuous, vname)));
  }

  // 7x1 - 2x2
  lf1_ = LinearFunctionPtr(new LinearFunction());
  lf1_->addTerm(vars_[0], 7.0);
  lf1_->addTerm(vars_[1], -2.0);
  
  // 2 x1*x2
  qf1_ = QuadraticFunctionPtr(new QuadraticFunction());
  qf1_->addTerm(VariablePair(vars_[0],vars_[1]), 2.0);
  
  // 3 x1^2 
  qf2_ = QuadraticFunctionPtr(new QuadraticFunction());
  qf2_->addTerm(VariablePair(vars_[0],vars_[0]), 3.0);
}


void FunctionTest::testEval()
{
  int error = 0;
  double x[2] = {1.0, 1.0};
  FunctionPtr f = (FunctionPtr) new Function(lf1_);
  CPPUNIT_ASSERT(f->eval(&x[0], &error)==5.0);
  CPPUNIT_ASSERT(0==error);

  double y[2] = {0.0, 1.0};
  CPPUNIT_ASSERT(f->eval(&y[0], &error)==-2.0);
  CPPUNIT_ASSERT(0==error);

  double z[4] = {1.0, 0.0, 0.0, 1.0};
  CPPUNIT_ASSERT(f->eval(&z[0], &error)==7.0);
  CPPUNIT_ASSERT(0==error);

  CPPUNIT_ASSERT(f->getNumVars()==2);
  CPPUNIT_ASSERT(lf1_->getNumTerms() == 2);

  f = (FunctionPtr) new Function(lf1_, qf1_);
  CPPUNIT_ASSERT(f->getNumVars()==2);
  CPPUNIT_ASSERT(f->eval(&x[0], &error)==7.0);
  CPPUNIT_ASSERT(0==error);
  CPPUNIT_ASSERT(f->eval(&y[0], &error)==-2.0);
  CPPUNIT_ASSERT(0==error);
  CPPUNIT_ASSERT(f->eval(&z[0], &error)==7.0);
  CPPUNIT_ASSERT(0==error);
}

void FunctionTest::testBilinearRecognize()
{
  FunctionPtr f1 = FunctionPtr(new Function(lf1_));
  CPPUNIT_ASSERT(f1->getType() != Bilinear);

  FunctionPtr f2 = FunctionPtr(new Function(lf1_, qf1_));
  CPPUNIT_ASSERT(f2->getType() == Bilinear);
}

void FunctionTest::testGetFixVarOffset() 
{
  FunctionPtr f1 = FunctionPtr(new Function(lf1_));
  CPPUNIT_ASSERT(f1->getFixVarOffset(vars_[0], 1.0) == 7.0); 

  FunctionPtr f2 = FunctionPtr(new Function(lf1_, qf2_));
  CPPUNIT_ASSERT(f2->getFixVarOffset(vars_[0], 2.0) == 26.0);
}

void FunctionTest::testEvalGradient()
{
  double y[] = {1.0, 1.0}; 
  double grad_f[] = {0.0, 0.0}; 
  double* x = &y[0];
  int error = 0;

  FunctionPtr f1 = FunctionPtr(new Function(lf1_)); 
  f1->evalGradient(x, grad_f, &error);
  CPPUNIT_ASSERT(grad_f[0] == 7.0); 
  CPPUNIT_ASSERT(grad_f[1] == -2.0);
  CPPUNIT_ASSERT(error == 0);
  
  grad_f[0] = 0.0;
  grad_f[1] = 0.0;
  FunctionPtr f2 = FunctionPtr(new Function(lf1_,qf1_));
  f2->evalGradient(x, grad_f, &error);    
  CPPUNIT_ASSERT(grad_f[0] == 9.0); 
  CPPUNIT_ASSERT(grad_f[1] == 0.0);
  CPPUNIT_ASSERT(error == 0);
}

void FunctionTest::tearDown()
{
  vars_.clear();
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
