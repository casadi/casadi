// 
//    MINOTAUR -- It's only 1/2 bull
// 
//    (C)opyright 2009 - 2014 The MINOTAUR Team.
// 



#include <cmath>
#include <fstream>
#include <iostream>

#include "MinotaurConfig.h"
#include "Constraint.h"
#include "Function.h"
#include "LinearFunctionUT.h"
#include "Objective.h"
#include "Variable.h"

CPPUNIT_TEST_SUITE_REGISTRATION(LinearFunctionTest);
CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(LinearFunctionTest, "LinearFunctionUT");

using namespace Minotaur;


void LinearFunctionTest::setUp()
{
  // This is instance on page 95 in Wolsey.
  //  max 4x1 - x2, 7x1 - 2 x2 <= 14, x2 <= 3, 2 x1 - 2 x2 <= 3
  //

  FunctionPtr f;
  ConstraintPtr cPtr;
  ObjectivePtr oPtr;
  LinearFunctionPtr lo, lf;

  instance_ = (ProblemPtr) new Problem();
  std::vector<VariablePtr> vars;

  vars.push_back(instance_->newVariable(0.0, INFINITY, Integer));
  vars.push_back(instance_->newVariable(0.0, 3.0, Integer));
  
  // 7x1 - 2x2 <= 14
  lf = (LinearFunctionPtr) new LinearFunction();
  lf->addTerm(vars[0], 7.0);
  lf->addTerm(vars[1], -2.0);

  f = (FunctionPtr) new Function(lf);
  cPtr = instance_->newConstraint(f, -INFINITY, 14.0, "cons0");
  
  // 2x1 - 2 x2
  lf = (LinearFunctionPtr) new LinearFunction();
  lf->addTerm(vars[0], 2.0);
  lf->addTerm(vars[1], -2.0);

  f = (FunctionPtr) new Function(lf);
  cPtr = instance_->newConstraint(f, -INFINITY, 3.0, "cons1");

  // max 4x1 - x2
  lo = (LinearFunctionPtr) new LinearFunction();
  lo->addTerm(vars[0], 4.0);
  lo->addTerm(vars[1], -1.0);
  f = (FunctionPtr) new Function(lo);


  oPtr = instance_->newObjective(f, 0.0, Maximize);
}


void LinearFunctionTest::tearDown()
{
  instance_->clear();
}


void LinearFunctionTest::testGetCoeffs()
{
  // get coefficients in each constraint
  ConstraintPtr cPtr;
  LinearFunctionPtr lf;
  VariableGroupConstIterator v_g_iter;
  ConstVariablePtr vPtr;


  // first constraint
  cPtr = instance_->getConstraint(0);
  lf = cPtr->getLinearFunction();

  v_g_iter = lf->termsBegin();
  vPtr = v_g_iter->first;
  CPPUNIT_ASSERT(vPtr->getId()==0);
  CPPUNIT_ASSERT(v_g_iter->second==7.0);

  ++v_g_iter;
  vPtr = v_g_iter->first;
  CPPUNIT_ASSERT(vPtr->getId()==1);
  CPPUNIT_ASSERT(v_g_iter->second==-2.0);

  // second constraint
  cPtr = instance_->getConstraint(1);
  lf = cPtr->getLinearFunction();

  v_g_iter = lf->termsBegin();
  vPtr = v_g_iter->first;
  CPPUNIT_ASSERT(vPtr->getId()==0);
  CPPUNIT_ASSERT(v_g_iter->second==2.0);

  ++v_g_iter;
  vPtr = v_g_iter->first;
  CPPUNIT_ASSERT(vPtr->getId()==1);
  CPPUNIT_ASSERT(v_g_iter->second==-2.0);
}


// test if linear parts of objective are correctly reported.
void LinearFunctionTest::testGetObj()
{
  ObjectivePtr oPtr;
  LinearFunctionPtr lf;
  VariableGroupConstIterator v_g_iter;
  ConstVariablePtr vPtr;

  oPtr = instance_->getObjective();

  CPPUNIT_ASSERT(oPtr->getObjectiveType()==Maximize);

  lf = oPtr->getLinearFunction();
  v_g_iter = lf->termsBegin();
  vPtr = v_g_iter->first;
  CPPUNIT_ASSERT(vPtr->getId()==0);
  CPPUNIT_ASSERT(v_g_iter->second==4.0);

  ++v_g_iter;
  vPtr = v_g_iter->first;
  CPPUNIT_ASSERT(vPtr->getId()==1);
  CPPUNIT_ASSERT(v_g_iter->second==-1.0);
}


void LinearFunctionTest::testOperations()
{
  std::vector<VariablePtr> vars;
  VariablePtr var;
  ConstVariablePtr vPtr;
  LinearFunctionPtr lf = (LinearFunctionPtr) new LinearFunction();
  LinearFunctionPtr lf1 = (LinearFunctionPtr) new LinearFunction();
  LinearFunctionPtr lf2;
  std::string vname = "common_var_name";

  var = (VariablePtr) new Variable(0, 0, 0.0, 1.0, Integer, vname);
  vars.push_back(var);

  var = (VariablePtr) new Variable(1, 1, 9.0, 1.0, Integer, vname);
  vars.push_back(var);

  var = (VariablePtr) new Variable(2, 2, -1.0, 1.0, Integer, vname);
  vars.push_back(var);

  var = (VariablePtr) new Variable(3, 3, -100.0, 100.0, Integer, vname);
  vars.push_back(var);

  // 2x0 - 6x1 + x2
  lf->addTerm(vars[0], 2.0);
  lf->addTerm(vars[1],-6.0);
  lf->addTerm(vars[2], 1.0);

  // 4x0 + 13x1 + 24.5x3
  lf1->addTerm(vars[0], 4.0);
  lf1->addTerm(vars[1], 13.0);
  lf1->addTerm(vars[3], 24.5);

  lf2 = lf + lf1;
  CPPUNIT_ASSERT(lf2->getWeight(vars[0]) == 6.0);
  CPPUNIT_ASSERT(lf2->getWeight(vars[1]) == 7.0);
  CPPUNIT_ASSERT(lf2->getWeight(vars[2]) == 1.0);
  CPPUNIT_ASSERT(lf2->getWeight(vars[3]) == 24.5);

  lf2 = lf - lf1;
  CPPUNIT_ASSERT(lf2->getWeight(vars[0]) == -2.0);
  CPPUNIT_ASSERT(lf2->getWeight(vars[1]) == -19.0);
  CPPUNIT_ASSERT(lf2->getWeight(vars[2]) == 1.0);
  CPPUNIT_ASSERT(lf2->getWeight(vars[3]) == -24.5);

  lf2 = lf - lf;
  CPPUNIT_ASSERT(lf2->getWeight(vars[0]) == 0.0); 
  CPPUNIT_ASSERT(lf2->getWeight(vars[3]) == 0.0); 

  // multiply by scalar.
  lf2 = 1.5*lf;
  CPPUNIT_ASSERT(lf2->getWeight(vars[0]) == 3.0); 
  CPPUNIT_ASSERT(lf2->getWeight(vars[1]) == -9.0); 
  CPPUNIT_ASSERT(lf2->getWeight(vars[3]) == 0.0); 

  // increment
  (*lf) += lf1;
  CPPUNIT_ASSERT(lf->getWeight(vars[0]) == 6.0); 
  CPPUNIT_ASSERT(lf->getWeight(vars[1]) == 7.0); 
  CPPUNIT_ASSERT(lf->getWeight(vars[2]) == 1.0);
  CPPUNIT_ASSERT(lf->getWeight(vars[3]) == 24.5);

  ConstLinearFunctionPtr lf6 = lf->clone();
  for (VariableGroupConstIterator it = lf6->termsBegin(); 
      it != lf6->termsEnd(); ++it) {
    CPPUNIT_ASSERT(lf6->getWeight(vars[0]) == 6.0); 
    CPPUNIT_ASSERT(lf6->getWeight(vars[1]) == 7.0); 
    CPPUNIT_ASSERT(lf6->getWeight(vars[2]) == 1.0);
    CPPUNIT_ASSERT(lf6->getWeight(vars[3]) == 24.5);
  }

}


void LinearFunctionTest::testFix()
{
  VariablePtr x0, x1, x2, x3;
  LinearFunctionPtr lf = (LinearFunctionPtr) new LinearFunction();
  double x[4] = {1.0, 5.0, -1.0, 11.0};

  x0 = (VariablePtr) new Variable(0, 0, 0.0, 1.0, Integer, "x0");

  x1 = (VariablePtr) new Variable(1, 1, 9.0, 1.0, Integer, "x1");

  x2 = (VariablePtr) new Variable(2, 2, -1.0, 1.0, Integer, "x2");

  x3 = (VariablePtr) new Variable(3, 3, -100.0, 100.0, Integer, "x3");

  // 2x0 - 6x1 + x2 + 9x3
  lf->addTerm(x0, 2.0);
  lf->addTerm(x1,-6.0);
  lf->addTerm(x2, 1.0);
  lf->addTerm(x3, 9.0);

  CPPUNIT_ASSERT(lf->eval(x) == (2.0*1.0 - 6.0*5.0 - 1.0*1.0 + 9.0*11.0));

  CPPUNIT_ASSERT(lf->getFixVarOffset(x0,0) == 0.0);
  CPPUNIT_ASSERT(lf->getFixVarOffset(x0,1.0) == 2.0);
  CPPUNIT_ASSERT(lf->getFixVarOffset(x1,3.0) == -18.0);
  CPPUNIT_ASSERT(lf->getFixVarOffset(x2,3.0) ==  3.0);
  CPPUNIT_ASSERT(lf->getFixVarOffset(x3,0.5) ==  4.5);

  lf->removeVar(x2, 0.0);
  CPPUNIT_ASSERT(lf->eval(x) ==  (2.0*1.0 - 6.0*5.0 + 9.0*11.0));
  lf->removeVar(x3, 0.0);
  CPPUNIT_ASSERT(lf->eval(x) ==  (2.0*1.0 - 6.0*5.0));

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
