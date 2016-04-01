// 
//     MINOTAUR -- It's only 1/2 bull
//
//     (C)opyright 2009 - 2014 The MINOTAUR Team.
//

#include "MinotaurConfig.h"
#include "AMPLCGraphUT.h"
#include "AMPLHessian.h"
#include "AMPLJacobian.h"
#include "CGraph.h"
#include "Constraint.h"
#include "Environment.h"
#include "Function.h"
#include "LinearFunction.h"
#include "Objective.h"
#include "Option.h"
#include "QuadraticFunction.h"
#include "Variable.h"

CPPUNIT_TEST_SUITE_REGISTRATION(AMPLCGraphUT);
CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(AMPLCGraphUT, "AMPLCGraphUT");
using namespace boost;
using namespace MINOTAUR_AMPL;
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

void AMPLCGraphUT::setUp()
{
  Minotaur::EnvPtr env = (Minotaur::EnvPtr) new Minotaur::Environment();
  env->getOptions()->findInt("ampl_log_level")->setValue(Minotaur::LogNone);
  env->setLogLevel(Minotaur::LogError);
  iface_ = (AMPLInterfacePtr) new AMPLInterface(env);
  env->getOptions()->findBool("use_native_cgraph")->setValue(true);
}


void AMPLCGraphUT::tearDown()
{
  delete iface_;
  inst_->clear();
}


void AMPLCGraphUT::testAllFuns()
{
  Minotaur::JacobianPtr jac;
  Minotaur::HessianOfLagPtr hess;

  inst_ = iface_->readInstance("instances/allfuns");
  //inst_->write(std::cout);
  CPPUNIT_ASSERT(inst_->getNumVars() == 74);
  inst_->setNativeDer();
  jac = inst_->getJacobian();
  hess = inst_->getHessian();
  CPPUNIT_ASSERT(hess->getNumNz() == 115);

}


void AMPLCGraphUT::testSize()
{
  inst_ = iface_->readInstance("instances/minlp_eg0");
  CPPUNIT_ASSERT(inst_->getNumCons() == 5);
  CPPUNIT_ASSERT(inst_->getNumVars() == 5);
}


void AMPLCGraphUT::testVariables()
{
  Minotaur::VariablePtr vPtr;
  std::string vName;

  inst_ = iface_->readInstance("instances/minlp_eg0");
  // check Ids
  for (Minotaur::UInt i=0; i<inst_->getNumVars(); ++i) {
    vPtr = inst_->getVariable(i);
    CPPUNIT_ASSERT(vPtr->getId() == i);

    // ampl reorders the variables, so we need to check each variable's name
    // and then perform further tests.
    vName = vPtr->getName();
    if (vName=="x0") {
       CPPUNIT_ASSERT(vPtr->getType() == Minotaur::Integer);
       CPPUNIT_ASSERT(vPtr->getLb() == -INFINITY);
       CPPUNIT_ASSERT(vPtr->getUb() ==  INFINITY);
       CPPUNIT_ASSERT(vPtr->getState() ==  Minotaur::NormalVar);
    } else if (vName=="x1") {
       // ampl converts binary to integer for nonlinear variables.
       // XXX: may be we should check bounds while reading the problem.
       CPPUNIT_ASSERT(vPtr->getType() == Minotaur::Integer);
       CPPUNIT_ASSERT(vPtr->getLb() ==  0);
       CPPUNIT_ASSERT(vPtr->getUb() ==  1);
       CPPUNIT_ASSERT(vPtr->getState() ==  Minotaur::NormalVar);
    } else if (vName=="x2") {
       CPPUNIT_ASSERT(vPtr->getType() == Minotaur::Continuous);
       CPPUNIT_ASSERT(vPtr->getLb() == -INFINITY);
       CPPUNIT_ASSERT(vPtr->getUb() ==  INFINITY);
       CPPUNIT_ASSERT(vPtr->getState() ==  Minotaur::NormalVar);
    } else if (vName=="x3") {
       CPPUNIT_ASSERT(vPtr->getType() == Minotaur::Continuous);
       CPPUNIT_ASSERT(vPtr->getLb() ==  4);
       CPPUNIT_ASSERT(vPtr->getUb() ==  10);
       CPPUNIT_ASSERT(vPtr->getState() ==  Minotaur::NormalVar);
    } else if (vName=="x4") {
       CPPUNIT_ASSERT(vPtr->getType() == Minotaur::Continuous);
       CPPUNIT_ASSERT(vPtr->getLb() ==  0);
       CPPUNIT_ASSERT(vPtr->getUb() ==  INFINITY);
       CPPUNIT_ASSERT(vPtr->getState() ==  Minotaur::NormalVar);
    } else {
       CPPUNIT_ASSERT_MESSAGE("variable has an unexpected name", false);
    }
  }
}


void AMPLCGraphUT::testConstraints()
{
  Minotaur::ConstraintPtr cPtr;
  std::string cName;
  Minotaur::LinearFunctionPtr lfPtr;
  Minotaur::QuadraticFunctionPtr qfPtr;
  Minotaur::NonlinearFunctionPtr nlfPtr;
                       // x2,x0,x1,x3,x4   (remember ampl changes order)
  double x[5] = {0, 0, 0, 0, 0};
  double y[5] = {7, 0, 2, 9, 6};
  double z[5] = {1,-1, 2,-4, 5};
  int error = 0;

  inst_ = iface_->readInstance("instances/minlp_eg0");
  for (Minotaur::UInt i=0; i<inst_->getNumCons(); ++i) {
    cPtr = inst_->getConstraint(i);
    CPPUNIT_ASSERT(cPtr->getId() == i);

    // ampl reorders the constraints, so we need to check each constraint's name
    // and then perform further tests.
    cName = cPtr->getName();
    if (cName=="cons0") {
      lfPtr = cPtr->getLinearFunction();
      if (lfPtr) {
        CPPUNIT_ASSERT(!"cons 0 should not have a linear function!");
      }

      qfPtr = cPtr->getQuadraticFunction();
      if (qfPtr) {
        CPPUNIT_ASSERT(!"cons 0 should not have a quadratic function!");
      } else {
        CPPUNIT_ASSERT(cPtr->getFunction()->eval(x, &error) == 0);
        CPPUNIT_ASSERT(cPtr->getFunction()->eval(y, &error) == 53);
        CPPUNIT_ASSERT(cPtr->getFunction()->eval(z, &error) == 6);
      }

      nlfPtr = cPtr->getNonlinearFunction();
      if (!nlfPtr) {
        CPPUNIT_ASSERT(!"cons 0 should have a nonlinear function!");
      }  

    } else if (cName=="cons1") {
      lfPtr = cPtr->getLinearFunction();
      if (lfPtr) {
        CPPUNIT_ASSERT(!"cons 1 should not have a linear function!");
      }

      CPPUNIT_ASSERT(cPtr->getFunction()->eval(x, &error) == 0);
      CPPUNIT_ASSERT(cPtr->getFunction()->eval(y, &error) == 0);
      CPPUNIT_ASSERT(cPtr->getFunction()->eval(z, &error) == 0);

      nlfPtr = cPtr->getNonlinearFunction();
      if (!nlfPtr) {
        CPPUNIT_ASSERT(!"cons 1 should have a nonlinear function!");
      } else {
        CPPUNIT_ASSERT(nlfPtr->eval(x, &error) == 0);
        CPPUNIT_ASSERT(0==error);
        CPPUNIT_ASSERT(nlfPtr->eval(y, &error) == 0);
        CPPUNIT_ASSERT(0==error);
        CPPUNIT_ASSERT(nlfPtr->eval(z, &error) == 0);
        CPPUNIT_ASSERT(0==error);
      } 

    } else if (cName=="cons2") {
      lfPtr = cPtr->getLinearFunction();
      if (!lfPtr) {
        CPPUNIT_ASSERT(!"cons 2 should have a linear function!");
      } else {
        CPPUNIT_ASSERT(lfPtr->eval(x) == 0);
        CPPUNIT_ASSERT(lfPtr->eval(y) == -5);
        CPPUNIT_ASSERT(lfPtr->eval(z) == 0);
      }

      qfPtr = cPtr->getQuadraticFunction();
      if (qfPtr) {
        CPPUNIT_ASSERT(!"cons 2 should not have a quadratic function!");
      }
      nlfPtr = cPtr->getNonlinearFunction();
      if (nlfPtr) {
        CPPUNIT_ASSERT(!"cons 2 should not have a nonlinear function!");
      } 

    } else if (cName=="cons3") {
      lfPtr = cPtr->getLinearFunction();
      if (!lfPtr) {
        CPPUNIT_ASSERT(!"cons 3 should have a linear function!");
      } else {
        CPPUNIT_ASSERT(lfPtr->eval(x) == 0);
        CPPUNIT_ASSERT(lfPtr->eval(y) == 9);
        CPPUNIT_ASSERT(lfPtr->eval(z) == 2);
      }
      qfPtr = cPtr->getQuadraticFunction();
      if (qfPtr) {
        CPPUNIT_ASSERT(!"cons 3 should not have a quadratic function!");
      }
      nlfPtr = cPtr->getNonlinearFunction();
      if (nlfPtr) {
        CPPUNIT_ASSERT(!"cons 3 should not have a nonlinear function!");
      } 

    } else if (cName=="cons4") {
      lfPtr = cPtr->getLinearFunction();
      if (!lfPtr) {
        CPPUNIT_ASSERT(!"cons 4 should have a linear function!");
      } else {
        CPPUNIT_ASSERT(lfPtr->eval(x) == 0);
        CPPUNIT_ASSERT(lfPtr->eval(y) == 15);
        CPPUNIT_ASSERT(lfPtr->eval(z) == 1);
      }

      qfPtr = cPtr->getQuadraticFunction();
      if (qfPtr) {
        CPPUNIT_ASSERT(!"cons 4 should not have a quadratic function!");
      }
      nlfPtr = cPtr->getNonlinearFunction();
      if (nlfPtr) {
        CPPUNIT_ASSERT(!"cons 4 should not have a nonlinear function!");
      } 

    } else {
       CPPUNIT_ASSERT_MESSAGE("constraint has an unexpected name", false);
    }
  }
}


void AMPLCGraphUT::testJacobian()
{
  //
  // order of variables in the instance is x2, x0, x1, x3, x4
  // order of constraints in the instance is cons0, 1, 2, 3, 4
  //
  // order of nonzeros is jacobian is column wise.
  //

  Minotaur::UInt iRow[12];
  Minotaur::UInt jCol[12];
  double values[12];
  int error = 0;

  // x2 appears in 3 rows, x0 in 4, x1 in 3, x3 in 1, x4 in 1.
  //Minotaur::UInt correctCol[12] = {0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 3, 4};
  //Minotaur::UInt correctRow[12] = {0, 2, 3, 0, 1, 2, 3, 0, 2, 3, 4, 4};

  Minotaur::UInt correctCol[12] = {0, 1, 2, 1, 0, 1, 2, 0, 1, 2, 3, 4};
  Minotaur::UInt correctRow[12] = {0, 0, 0, 1, 2, 2, 2, 3, 3, 3, 4, 4};

  double x0[5] = {0, 0, 0, 0, 0};
  double x1[5] = {1, 2, 3, 4, 5};

  //Minotaur::UInt correctCol[12] = {0, 1, 2,  1,  0, 1, 2, 0, 1, 2, 3, 4};
  //Minotaur::UInt correctRow[12] = {0, 0, 0,  1,  2, 2, 2, 3, 3, 3, 4, 4};
  double correctValues0[12] = {0, 0, 0,  0, -1, 1, 1, 1, 1, 1, 1, 1};
  double correctValues1[12] = {2, 4, 6, 16, -1, 1, 1, 1, 1, 1, 1, 1};

  inst_ = iface_->readInstance("instances/minlp_eg0");
  // set jacobian
  inst_->setNativeDer();
  Minotaur::JacobianPtr jacPtr = inst_->getJacobian();

  // check jacobian non-zeros
  CPPUNIT_ASSERT(jacPtr->getNumNz() == 12);

  // check indices of non-zeros
  jacPtr->fillRowColIndices(iRow, jCol);

  for (Minotaur::UInt i=0; i<jacPtr->getNumNz(); ++i) {
    CPPUNIT_ASSERT(iRow[i] == correctRow[i]);
    CPPUNIT_ASSERT(jCol[i] == correctCol[i]);
  }

  // get jacobian at [0 0 0 0 0]
  jacPtr->fillRowColValues(x0, values, &error);
  CPPUNIT_ASSERT(0==error);
  for (Minotaur::UInt i=0; i<jacPtr->getNumNz(); ++i) {
    CPPUNIT_ASSERT(fabs(values[i] - correctValues0[i]) < 1e-7);
  }

  // get jacobian at [1 2 3 4 5]
  jacPtr->fillRowColValues(x1, values, &error);
  CPPUNIT_ASSERT(0==error);
  for (Minotaur::UInt i=0; i<jacPtr->getNumNz(); ++i) {
    CPPUNIT_ASSERT(fabs(values[i] - correctValues1[i]) < 1e-7);
  }

}


void AMPLCGraphUT::testHessian()
{
  // variables are ordered as x2, x0, x1, x3, x4
  Minotaur::HessianOfLagPtr h;
  Minotaur::UInt row[5];
  Minotaur::UInt col[5];
  Minotaur::UInt rowexp[5] = {0, 1, 2, 2, 3};
  Minotaur::UInt colexp[5] = {0, 1, 0, 2, 1};
  double x[5] = {0.0, 0.0, 0.0, 0.0, 0.0};
  double hval[5] = {0.0, 0.0, 0.0, 0.0, 0.0};
  double con_mult[5] = {0.0, 0.0, 0.0, 0.0, 0.0};
  double exph[5];
  int error = 0;

  inst_ = iface_->readInstance("instances/minlp_eg0");
  inst_->setNativeDer();

  h = inst_->getHessian();
  CPPUNIT_ASSERT(h->getNumNz() == 5);
  h->fillRowColIndices(row,col);
  for (int i=0; i<5; ++i) {
    CPPUNIT_ASSERT(row[i] == rowexp[i]);
    CPPUNIT_ASSERT(col[i] == colexp[i]);
  }

  h->fillRowColValues(x, 1.0, con_mult, hval, &error);
  exph[0] =0.0; exph[1] = 0.0; exph[2] = 1.0; exph[3] = 0.0; exph[4] = 1.0;
  for (int i=0; i<5; ++i) {
    CPPUNIT_ASSERT(fabs(hval[i] - exph[i])<1e-10);
  }

  memset(hval, 0, sizeof(double));
  con_mult[0] = -1.0;
  exph[0] = -2.0; exph[1] = -2.0; exph[2] =  0.0; exph[3] =-2.0; exph[4] = 0.0;
  h->fillRowColValues(x, 0.0, con_mult, hval, &error);
  for (int i=0; i<5; ++i) {
    CPPUNIT_ASSERT(fabs(hval[i] - exph[i])<1e-10);
  }

  memset(hval, 0, sizeof(double));
  con_mult[0] =  0.0; con_mult[1] =  5.5;
  exph[0] =  0.0; exph[1] =  11.0; exph[2] =  0.0; exph[3] = 0.0; exph[4] = 0.0;
  h->fillRowColValues(x, 0.0, con_mult, hval, &error);
  for (int i=0; i<5; ++i) {
    CPPUNIT_ASSERT(fabs(hval[i] - exph[i])<1e-10);
  }

  memset(hval, 0, sizeof(double));
  con_mult[0] =  1.0; con_mult[1] =  5.5;
  h->fillRowColValues(x, -3.0, con_mult, hval, &error);
  exph[0] =  2.0; exph[1] =  13.0; exph[2] =  -3.0; exph[3] = 2.0; exph[4] = -3.0;
  for (int i=0; i<5; ++i) {
    CPPUNIT_ASSERT(fabs(hval[i] - exph[i])<1e-10);
  }


  memset(hval, 0, sizeof(double));
  x[1] = 7.0;
  h->fillRowColValues(x, -3.0, con_mult, hval, &error);
  exph[0] =  2.0; exph[1] =  244.0; exph[2] =  -3.0; exph[3] = 2.0; exph[4] = -3.0;
  for (int i=0; i<5; ++i) {
    //std::cout << hval[i] << " " << exph[i] << std::endl;
    CPPUNIT_ASSERT(fabs(hval[i] - exph[i])<1e-10);
  }
}


void AMPLCGraphUT::testObjective()
{
  Minotaur::ObjectivePtr oPtr;
  std::string oName;
  Minotaur::NonlinearFunctionPtr nlfPtr;
  Minotaur::QuadraticFunctionPtr qfPtr;

                       // x2,x0,x1,x3,x4   (remember ampl changes order)
  double y[5] = {7, 0, 2, 9, 6};
  double z[5] = {1,-1, 2,-4, 5};
  int error = 0;
  inst_ = iface_->readInstance("instances/minlp_eg0");

  oPtr = inst_->getObjective();
  if (!oPtr) {
    CPPUNIT_ASSERT(!"objective is missing!");
  } else {
    oName = oPtr->getName();
    CPPUNIT_ASSERT(oName == "obj");

    // test eval
    nlfPtr = oPtr->getNonlinearFunction();
    //CPPUNIT_ASSERT(!nlfPtr); // nlfPtr is NULL
    //qfPtr = oPtr->getQuadraticFunction();
    CPPUNIT_ASSERT(oPtr->getFunction()->eval(y, &error) == 20);
    CPPUNIT_ASSERT(oPtr->getFunction()->eval(z, &error) == 11);
  }
}


void AMPLCGraphUT::testObjectiveGradient()
{
  Minotaur::ObjectivePtr oPtr;
  int error;

  // test points
                       // x2,x0,x1,x3,x4   (remember ampl changes order)
  double x[5] = {0, 0, 0, 0, 0};
  double y[5] = {7, 0, 2, 9, 6};
  double z[5] = {1,-1, 2,-4, 5};

  // expected results for the whole function.
  double gx[5] = {0, 0, 0, 0, 1};
  double gy[5] = {2, 9, 7, 0, 1};
  double gz[5] = {2, -4, 1, -1, 1};
  //double gz[5] = {1, 2, 3, 4, 0};

  // acually calculated values for just the quadratic parts
  double gradient[5] = {0, 0, 0, 0, 0};
  inst_ = iface_->readInstance("instances/minlp_eg0");
  oPtr = inst_->getObjective();
  if (!oPtr) {
    CPPUNIT_ASSERT(!"objective is missing!");
    return;
  }

  // check values for the whole function.
  std::fill(gradient, gradient+5, 0);
  oPtr->eval(x, &error);
  oPtr->evalGradient(x, gradient, &error);
  for (Minotaur::UInt i=0; i<5; ++i) {
    CPPUNIT_ASSERT(fabs(gradient[i] - gx[i]) < 1e-7);
  }
  std::fill(gradient, gradient+5, 0);

  oPtr->eval(y, &error);
  oPtr->evalGradient(y, gradient, &error);
  for (Minotaur::UInt i=0; i<5; ++i) {
    CPPUNIT_ASSERT(fabs(gradient[i] - gy[i]) < 1e-7);
  }

  std::fill(gradient, gradient+5, 0);
  oPtr->eval(z, &error);
  oPtr->evalGradient(z, gradient, &error);
  for (Minotaur::UInt i=0; i<5; ++i) {
    CPPUNIT_ASSERT(fabs(gradient[i] - gz[i]) < 1e-7);
  }
}


void AMPLCGraphUT::testNl()
{
  Minotaur::EnvPtr env = (Minotaur::EnvPtr) new Minotaur::Environment();
  env->getOptions()->findBool("use_native_cgraph")->setValue(true);
  env->setLogLevel(Minotaur::LogError);
  AMPLInterface iface(env);
  double x[5] = {1.0,5,1.0,5,1.0};
  double h[4] = {0.0,0.0,0.0,0.0};
  Minotaur::UInt rows[4] = {0,0,0,0};
  Minotaur::UInt cols[4] = {0,0,0,0};
  int err = 0;
  double m[2] = {1.0, 1.0};


  inst_ = iface.readInstance("instances/hess");
  inst_->setNativeDer();
  inst_->getHessian()->fillRowColIndices(rows, cols);
  h[0] = h[1] = h[2] = h[3] = 0.0;
  inst_->getHessian()->fillRowColValues(x, 0.0, m, h, &err);
  CPPUNIT_ASSERT(0==err);
  CPPUNIT_ASSERT(3==inst_->getHessian()->getNumNz());
  CPPUNIT_ASSERT(fabs(h[0]+0.138704)<1e-7);
  CPPUNIT_ASSERT(fabs(h[1]-0.0276854)<1e-7);
  CPPUNIT_ASSERT(fabs(h[2]+0.00552603)<1e-7);

  inst_->subst(inst_->getVariable(1), inst_->getVariable(3), 1.0);
  inst_->subst(inst_->getVariable(2), inst_->getVariable(4), 1.0);
  
  inst_->setNativeDer();
  inst_->getHessian()->fillRowColIndices(rows, cols);
  h[0] = h[1] = h[2] = h[3] = 0.0;
  inst_->getHessian()->fillRowColValues(x, 0.0, m, h, &err);
  CPPUNIT_ASSERT(0==err);
  CPPUNIT_ASSERT(3==inst_->getHessian()->getNumNz());
  CPPUNIT_ASSERT(fabs(h[0]+0.138704)<1e-7);
  CPPUNIT_ASSERT(fabs(h[1]-0.0276854)<1e-7);
  CPPUNIT_ASSERT(fabs(h[2]+0.00552603)<1e-7);
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
