// 
//     MINOTAUR -- It's only 1/2 bull
// 
//     (C)opyright 2009 - 2014 The MINOTAUR Team.
// 

#include <cmath>

// define minotaur specific definitions
#include "MinotaurConfig.h"
#include "Environment.h"
#include "Function.h"
#include "IpoptEngineUT.h"
#include "IpoptEngine.h"
#include "Logger.h"
#include "NonlinearFunction.h"

CPPUNIT_TEST_SUITE_REGISTRATION(IpoptEngineUT);
CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(IpoptEngineUT, "IpoptEngineUT");

using namespace Minotaur;

void IpoptEngineUT::setUp()
{

}

void IpoptEngineUT::tearDown()
{

}

void IpoptEngineUT::createInstance_()
{
  /* 
   * (Nocedal & Wright, Chapter 12, page 360)
   * min x0x1
   * s.t.
   * x0^2 + x1^2 = 1
   */
  VariablePtr x0, x1;
  FunctionPtr fPtr;
  ConstraintPtr cPtr;
  ObjectivePtr oPtr;
  myNLFun0Ptr obj_f;
  myNLFun1Ptr con0_f;

  /* create instance and add variables */
  instance_ = (ProblemPtr) new Problem();

  /* x0 */
  x0 = instance_->newVariable(-10, 10, Continuous);

  /* x1 */
  x1 = instance_->newVariable(-10, 10, Continuous);

  /* objective */
  obj_f = (myNLFun0Ptr) new myNLFun0();
  fPtr = (FunctionPtr) new Function(obj_f);
  oPtr = instance_->newObjective(fPtr, 0.0, Minimize);
  

  /* constraint */
  con0_f = (myNLFun1Ptr) new myNLFun1();
  fPtr = (FunctionPtr) new Function(con0_f);
  cPtr = instance_->newConstraint(fPtr, 1.0, 1.0, "cons0");
  

}

void IpoptEngineUT::testGetObjVal()
{
  EnvPtr env = (EnvPtr) new Environment();
  double initial_pt[2] = {1, -0};
  createInstance_();

  //create a new engine
  IpoptEnginePtr ipopt_e = (IpoptEnginePtr) new IpoptEngine(env);

  // create my own jacobian
  myJacPtr jPtr = (myJacPtr) new myJac();

  // create my own hessian
  myHessPtr hPtr = (myHessPtr) new myHess();

  instance_->setJacobian(jPtr);
  instance_->setHessian(hPtr);

  // set an initial point
  instance_->setInitialPoint(initial_pt);

  LoggerPtr logger = (LoggerPtr) new Logger(LogDebug);
  instance_->setLogger(logger);

  //load the problem
  ipopt_e->load(instance_);

  //solve
  ipopt_e->solve();

  // get status
  EngineStatus status = ipopt_e->getStatus();
  CPPUNIT_ASSERT(status==ProvenLocalOptimal);

  // get objvalue
  double value = ipopt_e->getSolutionValue();
  CPPUNIT_ASSERT(fabs(value+0.5) < 1e-7);

  // get solution
  const double *x = ipopt_e->getSolution()->getPrimal();
  CPPUNIT_ASSERT(fabs(x[0]-0.7071067812) < 1e-6);
  CPPUNIT_ASSERT(fabs(x[1]+0.7071067812) < 1e-6);
}

// ------------------------------------------------------------------------- //
// ------------------------------------------------------------------------- //

UInt myNLFun0::getGradientNzCount()
{
  // gradient of our objective is
  // [x1 x0]
  return 2;
}

UInt myNLFun0::getHessianNzCount()
{
  // x0x1 
  // 0  .   
  // 1  0  
  return 1;
}

double myNLFun0::eval(const double *x, int *error) 
{
  //std::cout << "\n\n obj = " << x[0]*x[1] << "\n\n";
  *error = 0;
  return x[0]*x[1];
}

void myNLFun0::evalGradient(const double *x, double *grad_f, int *error)
{
  *error = 0;
  grad_f[0] = x[1];
  grad_f[1] = x[0];
}


// ------------------------------------------------------------------------- //
// ------------------------------------------------------------------------- //

UInt myNLFun1::getGradientNzCount()
{
  // gradient of our function is
  // [2x0 2x1]
  return 2;
}


UInt myNLFun1::getHessianNzCount()
{
  // x0^2 + x1^2
  // 2  .   
  // 0  2  
  return 2;
}

double myNLFun1::eval(const double *x, int *error)
{
  *error = 0;
  return x[0]*x[0] + x[1]*x[1];
}

// ------------------------------------------------------------------------- //
// ------------------------------------------------------------------------- //
// implement my own jacobian, since i know what it looks like.
myJac::myJac()
{
  nrows_ = 1;
  ncols_ = 2;
}


UInt myJac::getNumNz() 
{
  // There is only one constraint.
  // x0^2 + x1^2 = 1
  // jacobian is
  // [2x0 2x1]
  return 2;
}


void myJac::fillRowColIndices(UInt *iRow, UInt *jCol)
{
  iRow[0] = 0;
  jCol[0] = 0;

  iRow[1] = 0;
  jCol[1] = 1;
}


void myJac::fillRowColValues(const double *x, double *values, int *error)
{
  *error = 0;
  values[0] = 2*x[0];
  values[1] = 2*x[1];
}


// ------------------------------------------------------------------------- //
// ------------------------------------------------------------------------- //
// implement my own hessian, since i know what it looks like.

myHess::myHess()
{
  // hessian of x0x1 = [0 0 0] (in lower triangular form)
  // hessian of x0^2 + x1^2 = [2 0 2] (in lower triangular form)
}

UInt myHess::getNumNz() const
{
  return 2;
}

void myHess::fillRowColIndices(UInt *iRow, UInt *jCol)
{
  iRow[0] = 0;
  jCol[0] = 0;

  iRow[1] = 1;
  jCol[1] = 1;
}

void myHess::fillRowColValues(const double *, double obj_mult,
                              const double *con_mult, double *values, int *error)
{
  *error = 0;
  values[0] = 2*obj_mult;
  values[1] = 2*con_mult[0];
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
