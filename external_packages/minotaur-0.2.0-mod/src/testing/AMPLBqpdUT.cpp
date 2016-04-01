// 
//     MINOTAUR -- It's only 1/2 bull
// 
//     (C)opyright 2009 - 2014 The MINOTAUR Team.
// 

#include "MinotaurConfig.h"
#include "AMPLBqpdUT.h"
#include "AMPLJacobian.h"
#include "AMPLHessian.h"
#include "BranchAndBound.h"
#include "BqpdEngine.h"
#include "Environment.h"
#include "Option.h"

CPPUNIT_TEST_SUITE_REGISTRATION(AMPLBqpdUT);
CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(AMPLBqpdUT, "AMPLBqpdUT");
using namespace boost;
using namespace MINOTAUR_AMPL;
using namespace std;

void AMPLBqpdUT::testNLP()
{
  Minotaur::EnvPtr env = (Minotaur::EnvPtr) new Minotaur::Environment();
  double solval = 0.0;
  env->getOptions()->findInt("ampl_log_level")->setValue(Minotaur::LogNone);
  iface_ = (AMPLInterfacePtr) new AMPLInterface(env);

  // read an instance
  Minotaur::ProblemPtr inst = iface_->readInstance("instances/hs021");

  // set starting point if any
  inst->setInitialPoint(iface_->getInitialPoint());

  inst->calculateSize();

  if (inst->isQuadratic() || inst->isLinear()) { 
    inst->setNativeDer();
  } else {
    // create the jacobian
    Minotaur::JacobianPtr jPtr = (AMPLJacobianPtr) new AMPLJacobian(iface_);
    inst->setJacobian(jPtr);

    // create the hessian
    Minotaur::HessianOfLagPtr hPtr = (AMPLHessianPtr) new AMPLHessian(iface_);
    inst->setHessian(hPtr);
  }

  //create a new engine
  Minotaur::BqpdEnginePtr bqpd_e = (Minotaur::BqpdEnginePtr) 
    new Minotaur::BqpdEngine(env);

  //load the problem
  bqpd_e->load(inst);

  //solve
  bqpd_e->solve();

  // get status
  Minotaur::EngineStatus status = bqpd_e->getStatus();
  CPPUNIT_ASSERT(status==Minotaur::ProvenLocalOptimal);

  // get objvalue
  solval = bqpd_e->getSolutionValue();
  CPPUNIT_ASSERT(fabs(solval+99.96) < 1e-7);

  inst->clear();
  delete iface_;
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
