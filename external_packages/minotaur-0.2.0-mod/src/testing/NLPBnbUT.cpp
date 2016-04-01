// 
//     MINOTAUR -- It's only 1/2 bull
// 
//     (C)opyright 2009 - 2014 The MINOTAUR Team.
// 

#include <cmath>

#include "MinotaurConfig.h"
#include "BranchAndBound.h"
#include "Environment.h"
#include "EngineFactory.h"
#include "Function.h"
#include "IntVarHandler.h"
#include "LinearHandler.h"
#include "LPEngine.h"
#include "LPProcessor.h"
#include "NLPBnbUT.h"
#include "NLPEngine.h"
#include "NodeIncRelaxer.h"
#include "Option.h"
#include "Relaxation.h"
#include "ReliabilityBrancher.h"

CPPUNIT_TEST_SUITE_REGISTRATION(NLPBnbUT);
CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(NLPBnbUT, "NLPBnbUT");

using namespace Minotaur;

ProblemPtr NLPBnbUT::createInstance_()
{
  // 
  // (Minor changes to Nocedal & Wright, Chapter 12, page 360)
  // min x0x1
  // s.t.
  // x0^2 + x1^2 = 2
  // 

  ProblemPtr instance;
  VariablePtr x0, x1;
  FunctionPtr fPtr;
  ConstraintPtr cPtr;
  ObjectivePtr oPtr;
  myNLFun3Ptr obj_f;
  myNLFun4Ptr con0_f;
  double initial_pt[2] = {0.5, 0.0};

  // create instance and add variables 
  instance = (ProblemPtr) new Problem();

  // variable x0.
  x0 = instance->newVariable(-1, 1, Integer);

  /* x1 */
  x1 = instance->newVariable(-1, 1, Integer);

  /* objective */
  obj_f = (myNLFun3Ptr) new myNLFun3();
  fPtr = (FunctionPtr) new Function(obj_f);
  oPtr = instance->newObjective(fPtr, 0.0, Minimize);
  

  /* constraint */
  con0_f = (myNLFun4Ptr) new myNLFun4();
  fPtr = (FunctionPtr) new Function(con0_f);
  cPtr = instance->newConstraint(fPtr, 2.0, 2.0);

  // create my own jacobian
  myJac2Ptr jPtr = (myJac2Ptr) new myJac2();

  // create my own hessian
  myHess2Ptr hPtr = (myHess2Ptr) new myHess2();

  instance->setJacobian(jPtr);
  instance->setHessian(hPtr);

  instance->setInitialPoint(initial_pt);

  return instance;
}

void NLPBnbUT::testNLPBnb()
{
  EnvPtr env = (EnvPtr) new Environment();
  env->getOptions()->findBool("modify_rel_only")->setValue(true);

  // create the instance.
  ProblemPtr instance = createInstance_();

  // create a NULL LP engine.
  LPEnginePtr lp_e = LPEnginePtr(); //NULL

  // create an IpoptEnginePtr.
  //IpoptEnginePtr nlp_e = (IpoptEnginePtr) new IpoptEngine(env);

  //BabInitPtr bit = (BabInitPtr) new BabInit(env, instance, lp_e, nlp_e, qp_e);
  //BranchAndBoundPtr bab = bit->createBab();
  //bab->setLogLevel(LogDebug2);
  //bab->solve();
  //CPPUNIT_ASSERT(fabs(bab->getUb() + 1.0) < 1e-6);

  instance->clear();

}

// ------------------------------------------------------------------------- //
// ------------------------------------------------------------------------- //

UInt myNLFun3::getGradientNzCount()
{
  // gradient of our objective is
  // [x1 x0]
  return 2;
}

UInt myNLFun3::getHessianNzCount()
{
  // x0x1 
  // 0  .   
  // 1  0  
  return 1;
}

double myNLFun3::eval(const double *x, int *error) 
{
  *error = 0;
  return x[0]*x[1];
}

void myNLFun3::evalGradient(const double *x, double *grad_f, int *error) 
{
  *error = 0;
  grad_f[0] = x[1];
  grad_f[1] = x[0];
}


// ------------------------------------------------------------------------- //
// ------------------------------------------------------------------------- //

UInt myNLFun4::getGradientNzCount()
{
  // gradient of our function is
  // [2x0 2x1]
  return 2;
}


UInt myNLFun4::getHessianNzCount()
{
  // x0^2 + x1^2
  // 2  .   
  // 0  2  
  return 2;
}

double myNLFun4::eval(const double *x, int *error) 
{
  *error = 0;
  return x[0]*x[0] + x[1]*x[1];
}

// ------------------------------------------------------------------------- //
// ------------------------------------------------------------------------- //
// implement my own jacobian, since i know what it looks like.
myJac2::myJac2()
{
  nrows_ = 1;
  ncols_ = 2;
}

UInt myJac2::getNumNz()
{
  // There is only one constraint.
  // x0^2 + x1^2 = 1
  // jacobian is
  // [2x0 2x1]
  return 2;
}

void myJac2::fillRowColIndices(UInt *iRow, UInt *jCol)
{
  iRow[0] = 0;
  jCol[0] = 0;

  iRow[1] = 0;
  jCol[1] = 1;
}

void myJac2::fillRowColValues(const double *x, double *values, int *error)
{
  *error = 0;
  values[0] = 2*x[0];
  values[1] = 2*x[1];
}


// ------------------------------------------------------------------------- //
// ------------------------------------------------------------------------- //
// implement my own hessian, since i know what it looks like.

myHess2::myHess2()
{
  // hessian of x0x1 = [0 1 0] (in lower triangular form)
  // hessian of x0^2 + x1^2 = [2 0 2] (in lower triangular form)
}

UInt myHess2::getNumNz() const
{
  return 3;
}

void myHess2::fillRowColIndices(UInt *iRow, UInt *jCol)
{
  iRow[0] = 0;
  jCol[0] = 0;

  iRow[1] = 1;
  jCol[1] = 0;

  iRow[2] = 1;
  jCol[2] = 1;
}

void myHess2::fillRowColValues(const double *, double obj_mult, 
                               const double *con_mult, double *values, 
                               int *error)
{
  *error = 0;
  values[0] = 2*con_mult[0];
  values[1] = 1*obj_mult;
  values[2] = 2*con_mult[0];
}

// ------------------------------------------------------------------------- //
// ------------------------------------------------------------------------- //
ProblemPtr NLPBnbUT::createInstance1_()
{
  // solve the following MINLP:
  //
  // min x0.x0 + x1.x1 - x0 - 1.5x1 
  //
  // s.t. -2x0 + 2x1 <= 1
  //     
  // x0 \in {0, 1}
  // x1 \in {0, 1}

  ProblemPtr instance;
  VariablePtr x0, x1;
  FunctionPtr fPtr;
  ConstraintPtr cPtr;
  ObjectivePtr oPtr;
  myNLFun5Ptr obj_f;
  myNLFun6Ptr con0_f;
  double initial_pt[2] = {0.0, 0.0};

  // create instance and add variables 
  instance = (ProblemPtr) new Problem();

  // variable x0.
  x0 = instance->newBinaryVariable();

  /* x1 */
  x1 = instance->newBinaryVariable();

  /* objective */
  obj_f = (myNLFun5Ptr) new myNLFun5();
  fPtr = (FunctionPtr) new Function(obj_f);
  oPtr = instance->newObjective(fPtr, 0.0, Minimize);
  

  /* constraint */
  con0_f = (myNLFun6Ptr) new myNLFun6();
  fPtr = (FunctionPtr) new Function(con0_f);
  cPtr = instance->newConstraint(fPtr, -INFINITY, 1.0);

  // create my own jacobian
  myJac3Ptr jPtr = (myJac3Ptr) new myJac3();

  // create my own hessian
  myHess3Ptr hPtr = (myHess3Ptr) new myHess3();

  instance->setJacobian(jPtr);
  instance->setHessian(hPtr);

  instance->setInitialPoint(initial_pt);

  return instance;
}

void NLPBnbUT::testNLPBnb1()
{
  EnvPtr env = (EnvPtr) new Environment();
  env->getOptions()->findBool("modify_rel_only")->setValue(true);
  env->getOptions()->findString("nlp_engine")->setValue("IPOPT");
  ProblemPtr p = createInstance1_();
  HandlerVector handlers;
  ReliabilityBrancherPtr br;
  RelaxationPtr rel;
  EnginePtr e;
  int err = 0;

  env->startTimer(err);
  BranchAndBound *bab = new BranchAndBound(env, p);

  IntVarHandlerPtr v_hand = (IntVarHandlerPtr) new IntVarHandler(env, p);
  LinearHandlerPtr l_hand = (LinearHandlerPtr) new LinearHandler(env, p);
  handlers.push_back(v_hand);
  handlers.push_back(l_hand);
  v_hand->setModFlags(false, true);
  l_hand->setModFlags(false, true);

  EngineFactory efac(env);
  e = efac.getNLPEngine();

  LPProcessorPtr nproc = (LPProcessorPtr) new LPProcessor(env, e, handlers);
  br= (ReliabilityBrancherPtr) new ReliabilityBrancher(env, handlers);
  br->setEngine(e);
  nproc->setBrancher(br);
  bab->setNodeProcessor(nproc);

  NodeIncRelaxerPtr nr = (NodeIncRelaxerPtr) new NodeIncRelaxer(env, handlers);
  bab->setNodeRelaxer(nr);

  rel = (RelaxationPtr) new Relaxation(p);
  rel->calculateSize();
  rel->setJacobian(p->getJacobian());
  rel->setHessian(p->getHessian());
  nr->setRelaxation(rel);
  nr->setEngine(e);
  nr->setModFlag(false);

  bab->shouldCreateRoot(false);
  bab->setLogLevel(LogNone);
  bab->solve();
  CPPUNIT_ASSERT(fabs(bab->getUb() + 0.5) < 1e-6);

  p->clear();
  delete bab;

  //CPPUNIT_ASSERT(!"implement me!");
}
// ------------------------------------------------------------------------- //
// ------------------------------------------------------------------------- //

UInt myNLFun5::getGradientNzCount()
{
  // gradient of our objective is
  // [2x0 - 1, 2x1 - 1.5]
  return 2;
}

UInt myNLFun5::getHessianNzCount()
{
  //   2    .   
  //   0    2
  return 2;
}

double myNLFun5::eval(const double *x, int *error) 
{
  *error = 0;
  return x[0]*x[0] + x[1]*x[1] - x[0] - 1.5*x[1];
}

void myNLFun5::evalGradient(const double *x, double *grad_f, int *error) 
{
  *error = 0;
  grad_f[0] = 2*x[0] - 1;
  grad_f[1] = 2*x[1] - 1.5;
}


// ------------------------------------------------------------------------- //
// ------------------------------------------------------------------------- //

UInt myNLFun6::getGradientNzCount()
{
  // gradient of our function is
  // [-2 2]
  return 2;
}


UInt myNLFun6::getHessianNzCount()
{
  // -2x0 + 2x1 <= 1
  // 0  .   
  // 0  0  
  return 0;
}

double myNLFun6::eval(const double *x, int *error) 
{
  *error = 0;
  return -2*x[0] + 2*x[1];
}

// ------------------------------------------------------------------------- //
// ------------------------------------------------------------------------- //
// implement my own jacobian, since i know what it looks like.
myJac3::myJac3()
{
  nrows_ = 1;
  ncols_ = 2;
}

UInt myJac3::getNumNz()
{
  // There is only one constraint.
  // -2x0 + 2x1 <= 1
  // jacobian is
  // [-2 2]
  return 2;
}

void myJac3::fillRowColIndices(UInt *iRow, UInt *jCol)
{
  iRow[0] = 0;
  jCol[0] = 0;

  iRow[1] = 0;
  jCol[1] = 1;
}

void myJac3::fillRowColValues(const double *, double *values, int *error)
{
  *error = 0;
  values[0] = -2;
  values[1] =  2;
}


// ------------------------------------------------------------------------- //
// ------------------------------------------------------------------------- //
// implement my own hessian, since i know what it looks like.

myHess3::myHess3()
{
}

UInt myHess3::getNumNz() const
{
  return 2;
}

void myHess3::fillRowColIndices(UInt *iRow, UInt *jCol)
{
  iRow[0] = 0;
  jCol[0] = 0;

  iRow[1] = 1;
  jCol[1] = 1;
}

void myHess3::fillRowColValues(const double *, 
    double obj_mult, const double *, double *values, int *error)
{
  *error = 0;
  values[0] = 2*obj_mult;
  values[1] = 2*obj_mult;
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
