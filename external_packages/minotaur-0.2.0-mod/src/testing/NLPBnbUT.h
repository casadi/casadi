// 
//     MINOTAUR -- It's only 1/2 bull
// 
//     (C)opyright 2009 - 2014 The MINOTAUR Team.
// 

#ifndef NLPBNBUT_H
#define NLPBNBUT_H

#include <cppunit/TestCase.h>
#include <cppunit/TestCaller.h>
#include <cppunit/TestSuite.h>
#include <cppunit/TestResult.h>
#include <cppunit/extensions/HelperMacros.h>

#include "Types.h"
#include "HessianOfLag.h"
#include "Jacobian.h"
#include "NonlinearFunction.h"
#include "Problem.h"

using namespace Minotaur;


class NLPBnbUT : public CppUnit::TestCase {

public:
  /// Setup for cpp-unit-test.
  NLPBnbUT(std::string name) : TestCase(name) {}

  /// Default constructor.
  NLPBnbUT() {}

  /// setUp necessary for cpp-unit-test.
  void setUp() { };         

  // /**
  // tearDown necessary for cpp-unit-test.
  // */
  void tearDown() { };   

  /**
   * solve the following MINLP:
   * 
   * (Simple modification of Nocedal & Wright, Chapter 12, page 360)
   * min x0x1
   * s.t.
   * x0^2 + x1^2 = 2
   * 
   * x0 \in {-1, 1}
   * x1 \in {-1, 1}
   * 
   */
  void testNLPBnb();

  /**
   * solve the following MINLP:
   * 
   * min x0.x0 + x1.x1 - x0 - 1.5y0 
   * 
   * s.t. -2x0 + 2x1 <= 1
   *     
   * x0 \in {0, 1}
   * x1 \in {0, 1}
   */
  void testNLPBnb1();

  CPPUNIT_TEST_SUITE(NLPBnbUT);
  CPPUNIT_TEST(testNLPBnb);
  CPPUNIT_TEST(testNLPBnb1);
  CPPUNIT_TEST_SUITE_END();

private:
  // /**
  // Create the problem we want to solve by adding constraints, variables,
  // nonlinear functions, jacobian of the constraint matrix and hessian of the
  // lagrangean.
  // */
  ProblemPtr createInstance_();
  ProblemPtr createInstance1_();
};

// ------------------------------------------------------------------------- //
// ------------------------------------------------------------------------- //
// objective function
class myNLFun3 : public NonlinearFunction {
  public:
    myNLFun3() {}
    NonlinearFunctionPtr cloneWithVars(VariableConstIterator,
                                       int *err) const 
    {*err = 1; return NonlinearFunctionPtr();};
    
    
    double eval(const double *x, int *error) ;
    void evalGradient(const double *x, double *grad_f, int *error) ;
    void evalHessian(const double, const double *, 
                     const LTHessStor *, double *, 
                     int *) {};
    void fillHessStor(LTHessStor *) {};
    void fillJac(const double*, double *, int*) {};
    void finalHessStor(const LTHessStor *) {};
    UInt getGradientNzCount();
    UInt getHessianNzCount();
    void getVars(VariableSet *) {};
    bool isZero() const { return true; }
    void multiply(const double) {};
    void prepJac(VarSetConstIter, VarSetConstIter) {};
};

// ------------------------------------------------------------------------- //
// ------------------------------------------------------------------------- //
// function in constraint 
class myNLFun4 : public NonlinearFunction {
  public:
    myNLFun4() {}
    NonlinearFunctionPtr cloneWithVars(VariableConstIterator,
                                       int *err) const 
    {*err = 1; return NonlinearFunctionPtr();};
    
    double eval(const double *x, int *error) ;
    void evalGradient(const double *, double *, int *) {}
    void evalHessian(const double, const double *, 
                     const LTHessStor *, double *, 
                     int *) {};
    void fillHessStor(LTHessStor *) {};
    void fillJac(const double*, double *, int*) {};
    void finalHessStor(const LTHessStor *) {};
    UInt getGradientNzCount();
    UInt getHessianNzCount();
    void getVars(VariableSet *) {};
    bool isZero() const { return true; }
    void multiply(const double) {};
    void prepJac(VarSetConstIter, VarSetConstIter) {};
};

// ------------------------------------------------------------------------- //
// ------------------------------------------------------------------------- //
// jacobian 
class myJac2 : public Jacobian {
  public:
    myJac2();
    UInt getNumNz();
    void fillRowColIndices(UInt *iRow, UInt *jCol);
    void fillRowColValues(const double *x, double *values, int *error);


  private:
    UInt nrows_;
    UInt ncols_;

};

// ------------------------------------------------------------------------- //
// ------------------------------------------------------------------------- //
// hessian 
class myHess2 : public HessianOfLag {
  public:
    myHess2();
    UInt getNumNz() const;
    void fillRowColIndices(UInt *iRow, UInt *jCol);
    void fillRowColValues(const double *x, double obj_mult,
                          const double *con_mult, double *values,  int *error);

};

// ------------------------------------------------------------------------- //
// ------------------------------------------------------------------------- //

// objective function
class myNLFun5 : public NonlinearFunction {
  public:
    myNLFun5() {}
    NonlinearFunctionPtr cloneWithVars(VariableConstIterator,
                                       int *err) const 
    {*err = 1; return NonlinearFunctionPtr();};
    
    
    double eval(const double *x, int *error) ;
    void evalGradient(const double *x, double *grad_f, int *error) ;
    void evalHessian(const double, const double *, 
                     const LTHessStor *, double *, 
                     int *) {};
    void fillHessStor(LTHessStor *) {};
    void fillJac(const double*, double *, int*) {};
    void finalHessStor(const LTHessStor *) {};
    UInt getGradientNzCount();
    UInt getHessianNzCount();
    void getVars(VariableSet *) {};
    bool isZero() const { return true; }
    void prepJac(VarSetConstIter, VarSetConstIter) {};
    void multiply(const double) {};
};

// ------------------------------------------------------------------------- //
// ------------------------------------------------------------------------- //
// function in constraint 
class myNLFun6 : public NonlinearFunction {
  public:
    myNLFun6() {}
    NonlinearFunctionPtr cloneWithVars(VariableConstIterator,
                                       int *err) const 
    {*err = 1; return NonlinearFunctionPtr();};
    
    double eval(const double *x, int *error) ;
    void evalGradient(const double *, double *, int *) {}
    void evalHessian(const double, const double *, 
                     const LTHessStor *, double *, 
                     int *) {};
    void fillHessStor(LTHessStor *) {};
    void fillJac(const double*, double *, int*) {};
    void finalHessStor(const LTHessStor *) {};
    UInt getGradientNzCount();
    UInt getHessianNzCount();
    void getVars(VariableSet *) {};
    bool isZero() const { return true; }
    void multiply(const double) {};
    void prepJac(VarSetConstIter, VarSetConstIter) {};
};

// ------------------------------------------------------------------------- //
// ------------------------------------------------------------------------- //
// jacobian 
class myJac3 : public Jacobian {
  public:
    myJac3();
    UInt getNumNz();
    void fillRowColIndices(UInt *iRow, UInt *jCol);
    void fillRowColValues(const double *x, double *values, int *error);


  private:
    UInt nrows_;
    UInt ncols_;

};

// ------------------------------------------------------------------------- //
// ------------------------------------------------------------------------- //
// hessian 
class myHess3 : public HessianOfLag {
  public:
    myHess3();
    UInt getNumNz() const;
    void fillRowColIndices(UInt *iRow, UInt *jCol);
    void fillRowColValues(const double *x, double obj_mult, 
                          const double *con_mult,  double *values, int *error);

};


typedef boost::shared_ptr<myNLFun3> myNLFun3Ptr;
typedef boost::shared_ptr<myNLFun4> myNLFun4Ptr;
typedef boost::shared_ptr<myJac2> myJac2Ptr;
typedef boost::shared_ptr<myHess2> myHess2Ptr;
typedef boost::shared_ptr<myNLFun5> myNLFun5Ptr;
typedef boost::shared_ptr<myNLFun6> myNLFun6Ptr;
typedef boost::shared_ptr<myJac3> myJac3Ptr;
typedef boost::shared_ptr<myHess3> myHess3Ptr;

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
