// 
//     MINOTAUR -- It's only 1/2 bull
// 
//     (C)opyright 2009 - 2014 The MINOTAUR Team.
// 

#include <cmath>

#include "MinotaurConfig.h"
#include "CGraphUT.h"
#include "CGraph.h"
#include "CNode.h"
#include "Problem.h"
#include "Variable.h"

CPPUNIT_TEST_SUITE_REGISTRATION(CGraphUT);
CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(CGraphUT, "CGraphUT");
using namespace Minotaur;

void CGraphUT::testIdentical()
{
  VariablePtr v0 = (VariablePtr) new Variable(0, 0, 0.0, 10.0, Continuous, "x0");

  CNode *n0, *n1;
  CGraphPtr cg1 = (CGraphPtr) new CGraph();
  CGraphPtr cg2 = (CGraphPtr) new CGraph();

  // x^2 = x^2
  n0 = cg1->newNode(v0);
  n1 = cg1->newNode(OpSqr, n0, 0);
  cg1->setOut(n1);
  cg1->finalize();

  n0 = cg2->newNode(v0);
  n1 = cg2->newNode(OpSqr, n0, 0);
  cg2->setOut(n1);
  cg2->finalize();

  CPPUNIT_ASSERT(cg2->isIdenticalTo(cg1)); 
  CPPUNIT_ASSERT(cg1->isIdenticalTo(cg2)); 
  CPPUNIT_ASSERT(cg1->isIdenticalTo(cg1)); 

  // x^2 != x^3
  cg1 = (CGraphPtr) new CGraph();
  cg2 = (CGraphPtr) new CGraph();

  n0 = cg1->newNode(v0);
  n1 = cg1->newNode(OpSqr, n0, 0);
  cg1->setOut(n1);
  cg1->finalize();

  n0 = cg2->newNode(v0);
  n1 = cg2->newNode(3);
  n1 = cg2->newNode(OpPowK, n0, n1);
  cg2->setOut(n1);
  cg2->finalize();

  CPPUNIT_ASSERT(cg2->isIdenticalTo(cg2)); 
  CPPUNIT_ASSERT(false == cg1->isIdenticalTo(cg2)); 

  // x^3 == x^3
  cg1 = (CGraphPtr) new CGraph();
  n0 = cg1->newNode(v0);
  n1 = cg1->newNode(3);
  n1 = cg1->newNode(OpPowK, n0, n1);
  cg1->setOut(n1);
  cg1->finalize();
  CPPUNIT_ASSERT(cg1->isIdenticalTo(cg2)); 
}


void CGraphUT::testLin()
{
  CNode *n0, *n1, *n2, *n3;
  CGraph cgraph;
  int error = 0;

  VariablePtr v0 = (VariablePtr) new Variable(0, 0, 0.0, 10.0, Continuous, "x0");
  VariablePtr v1 = (VariablePtr) new Variable(1, 1, 0.0, 10.0, Continuous, "x1");
  VariablePtr v2 = (VariablePtr) new Variable(2, 2, 0.0, 10.0, Continuous, "x2");
  VariablePtr v3 = (VariablePtr) new Variable(3, 3, 0.0, 10.0, Continuous, "x3");

  double x[4] = {1.0, 2.0, 5.0, 7.0};
  double g[4] = {0.0, 0.0, 0.0, 0.0};
  double gexp[4] = {1.0, 1.0, 5.4, -1.0};


  n0 = cgraph.newNode(v0);
  n1 = cgraph.newNode(v0);
  CPPUNIT_ASSERT(n0 == n1);

  n1 = cgraph.newNode(v1);
  n2 = cgraph.newNode(v1);
  CPPUNIT_ASSERT(n2 == n1);
  CPPUNIT_ASSERT(n0 != n1);

  n2 = cgraph.newNode(OpPlus, n0, n1);
  n3 = cgraph.newNode(5.4);
  n0 = cgraph.newNode(v2);
  n3 = cgraph.newNode(OpMult, n3, n0);
  n3 = cgraph.newNode(OpPlus, n2, n3); // n3 = x0 + x1 + 5.4*x2
  n0 = cgraph.newNode(v3);
  n0 = cgraph.newNode(OpMinus, n3, n0); // n0 = x0 + x1 + 5.4*x2 - x3

  cgraph.setOut(n0);
  cgraph.finalize();

  CPPUNIT_ASSERT(cgraph.numVars() == 4);
  CPPUNIT_ASSERT(cgraph.getType() == Linear);
  CPPUNIT_ASSERT(fabs(cgraph.eval(x, &error) - 23.0)<1e-10);
  CPPUNIT_ASSERT(0==error);
  cgraph.evalGradient(x, g, &error);
  CPPUNIT_ASSERT(0==error);
  for (UInt i=0; i<4; ++i) {
    CPPUNIT_ASSERT(fabs(g[i]-gexp[i])<1e-10);
  }


  x[0] = 0.0; x[1] = 2.0; x[2] = 0.0; x[3] = 7.0; 
  g[0] = 0.0; g[1] = 0.0; g[2] = 0.0; g[3] = 0.0; 
  cgraph.evalGradient(x, g, &error);
  for (UInt i=0; i<4; ++i) {
    CPPUNIT_ASSERT(fabs(g[i]-gexp[i])<1e-10);
  }
  CPPUNIT_ASSERT(0==error);
  CPPUNIT_ASSERT(fabs(n0->getVal()+5.0)<1e-10);
  CPPUNIT_ASSERT(fabs(cgraph.eval(x,&error)+5.0)<1e-10);
}


void CGraphUT::testQuad()
{

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
