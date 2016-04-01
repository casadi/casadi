// 
//     MINOTAUR -- It's only 1/2 bull
// 
//     (C)opyright 2009 - 2014 The MINOTAUR Team.
// 

#ifndef OBJECTIVEUT_H
#define OBJECTIVEUT_H

#include <cppunit/TestCase.h>
#include <cppunit/TestCaller.h>
#include <cppunit/TestSuite.h>
#include <cppunit/TestResult.h>
#include <cppunit/extensions/HelperMacros.h>

#include <Objective.h>
#include "NonlinearFunction.h"

using namespace Minotaur;


class ObjectiveUT : public CppUnit::TestCase {
  public:
    ObjectiveUT(std::string name) : TestCase(name) {}
    ObjectiveUT() {}

    void testGetVal();
    void setUp();
    void tearDown() { };

    CPPUNIT_TEST_SUITE(ObjectiveUT);
    CPPUNIT_TEST(testGetVal);
    CPPUNIT_TEST_SUITE_END();

  private:
    ObjectivePtr objective_;
    VariablePtr v0_, v1_, v2_, v3_;
};

// ------------------------------------------------------------------------- //
// ------------------------------------------------------------------------- //
// our own nonlinear function
class myNLFun2 : public NonlinearFunction {
  public:
    myNLFun2() {}
    NonlinearFunctionPtr cloneWithVars(VariableConstIterator, int *err) const 
    {*err = 1; return NonlinearFunctionPtr();};

    double eval(const double *x, int *error) ;
    void evalGradient(const double *, double *, int *) {};
    void evalHessian(const double, const double *, 
                     const LTHessStor *, double *, 
                     int *) {};
    void fillHessStor(LTHessStor *) {};
    void fillJac(const double*, double *, int*) {};
    void finalHessStor(const LTHessStor *) {};
    void getVars(VariableSet *) {};
    void multiply(const double) {};
    void prepJac(VarSetConstIter, VarSetConstIter) {};
};
typedef boost::shared_ptr<myNLFun2> myNLFun2Ptr;

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
