// 
//     MINOTAUR -- It's only 1/2 bull
//
//     (C)opyright 2009 - 2014 The MINOTAUR Team.
//

#include "MinotaurConfig.h"
#include "TransformerUT.h"

#include "Environment.h"
#include "Handler.h"
#include "Option.h"
#include "Transformer.h"
#include "Types.h"

CPPUNIT_TEST_SUITE_REGISTRATION(TransformerUT);
CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(TransformerUT, "TransformerUT");

using namespace boost;
using namespace Minotaur;
using namespace std;

// minlp_eg0.mod
// # simple example of a minlp
// 
// var x0 integer;
// var x1 binary;
// var x2;
// var x3 >=4, <=10;
// var x4 >=0;
// 
// minimize obj: x0*x3 + x1*x2 + x4;
// 
// subject to cons0:
//         x0*x0 + x1*x1 + x2*x2 = 1;
// 
// subject to cons1:
//         x0^3 + x0^2 <= 100;
// 
// subject to cons2:
//         x0 + x1 - x2 >= 0;
// 
// subject to cons3:
//         x0 + x1 + x2 <= 3;
// 
// subject to cons4:
//         x3 + x4 <= 10;

// keep in mind that ampl changes the order of variables and constraints.

void TransformerUT::setUp()
{
  env_ = (EnvPtr) new Environment();
  env_->getOptions()->findInt("ampl_log_level")->setValue(Minotaur::LogNone);
  env_->setLogLevel(LogError);
  iface_ = (MINOTAUR_AMPL::AMPLInterfacePtr)
    new MINOTAUR_AMPL::AMPLInterface(env_);
  env_->getOptions()->findBool("use_native_cgraph")->setValue(true);
  inst_ = iface_->readInstance("instances/minlp_eg0");
}


void TransformerUT::tearDown()
{
  delete iface_;
  inst_->clear();
}


void TransformerUT::testSize()
{
  CPPUNIT_ASSERT(inst_->getNumCons() == 5);
  CPPUNIT_ASSERT(inst_->getNumVars() == 5);
}


void TransformerUT::testMakeObjectiveLinear()
{
  HandlerVector handlers;
  ProblemPtr newp = (ProblemPtr) new Problem();
  TransPoly t(env_, inst_);


  //std::cout << "Before: " << std::endl;
  //inst_->write(std::cout);

  //t.makeObjectiveLinear();

  //std::cout << "After: " << std::endl;
  //inst_->write(std::cout);

  CPPUNIT_ASSERT(inst_->getNumCons() == 5);
  CPPUNIT_ASSERT(inst_->getNumVars() == 5);
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
