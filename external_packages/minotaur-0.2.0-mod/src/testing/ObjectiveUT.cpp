// 
//     MINOTAUR -- It's only 1/2 bull
// 
//     (C)opyright 2009 - 2014 The MINOTAUR Team.
// 

#include "MinotaurConfig.h"
#include "Function.h"
#include "LinearFunction.h"
#include "Objective.h"
#include "ObjectiveUT.h"
#include "Variable.h"


CPPUNIT_TEST_SUITE_REGISTRATION(ObjectiveUT);
CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(ObjectiveUT, "ObjectiveUT");

using namespace Minotaur;

void ObjectiveUT::setUp()
{
  // XXX: add a quadratic
  // 3x0 + 10x1 - x3 + 6x0x1x2 + x3^3 + 2.0
  std::string vname = "common_var_name";

  // Define variables
  v0_ = (VariablePtr) new Variable (0, 0, -1.0, 1.0, Continuous, vname);
  v1_ = (VariablePtr) new Variable (1, 1, -1.0, 1.0, Continuous, vname);
  v2_ = (VariablePtr) new Variable (2, 2, -1.0, 1.0, Continuous, vname);
  v3_ = (VariablePtr) new Variable (3, 3, -1.0, 1.0, Continuous, vname);

  LinearFunctionPtr lPtr = (LinearFunctionPtr) new LinearFunction();
  lPtr->addTerm(v0_, 3.0);
  lPtr->addTerm(v1_, 10.0);
  lPtr->addTerm(v3_, -1.0);

  QuadraticFunctionPtr qPtr = QuadraticFunctionPtr(); //NULL
  myNLFun2Ptr nlPtr = (myNLFun2Ptr) new myNLFun2(); 

  FunctionPtr fPtr = (FunctionPtr) new Function(lPtr, qPtr, nlPtr);

  objective_ = (ObjectivePtr) new Objective(fPtr, 2.0, Minimize);
}


void ObjectiveUT::testGetVal()
{
  double x[4] = {1.0, 1.0, 1.0, 1.0};
  int error = 0;

  CPPUNIT_ASSERT(objective_->getFunction()->eval(x, &error) == 19.0);
  CPPUNIT_ASSERT(0==error);
}

// ------------------------------------------------------------------------ //
// ------------------------------------------------------------------------ //
// our own nonlinear function
double myNLFun2::eval(const double *x, int *error) 
{
  *error = 0;
  return 6*x[0]*x[1]*x[2] + x[3]*x[3]*x[3];
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
