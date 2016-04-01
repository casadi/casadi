// 
//     MINOTAUR -- It's only 1/2 bull
// 
//     (C)opyright 2009 - 2014 The MINOTAUR Team.
// 

#include "MinotaurConfig.h"
#include "AMPLIpoptUT.h"

#include "AMPLHessian.h"
#include "AMPLJacobian.h"
#include "BranchAndBound.h"
#include "Environment.h"
#include "IpoptEngine.h"
#include "Option.h"
#undef F77_FUNC_

CPPUNIT_TEST_SUITE_REGISTRATION(AMPLIpoptUT);
CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(AMPLIpoptUT, "AMPLIpoptUT");

using namespace MINOTAUR_AMPL;

void AMPLIpoptUT::testNLP()
{
  Minotaur::EnvPtr env = (Minotaur::EnvPtr) new Minotaur::Environment();
  env->getOptions()->findInt("ampl_log_level")->setValue(Minotaur::LogNone);
  iface_ = (AMPLInterfacePtr) new AMPLInterface(env);

  // read an instance
  Minotaur::ProblemPtr inst = iface_->readInstance("instances/3pk");

  // set starting point if any
  inst->setInitialPoint(iface_->getInitialPoint());

  inst->calculateSize();

  if (inst->isQuadratic()) { 
    inst->setNativeDer();
  } else {
    // create the jacobian
    AMPLJacobianPtr jPtr = (AMPLJacobianPtr) new AMPLJacobian(iface_);
    inst->setJacobian(jPtr);

    // create the hessian
    AMPLHessianPtr hPtr = (AMPLHessianPtr) new AMPLHessian(iface_);
    inst->setHessian(hPtr);
  }

  //create a new engine
  Minotaur::IpoptEnginePtr ipopt_e = (Minotaur::IpoptEnginePtr) 
    new Minotaur::IpoptEngine(env);

  //load the problem
  ipopt_e->load(inst);

  //solve
  ipopt_e->solve();

  // get status
  Minotaur::EngineStatus status = ipopt_e->getStatus();
  CPPUNIT_ASSERT(status==Minotaur::ProvenLocalOptimal);

  // get objvalue
  double value = ipopt_e->getSolutionValue();
  CPPUNIT_ASSERT(fabs(value-1.7201185) < 1e-7);

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
