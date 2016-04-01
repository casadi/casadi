// 
//    MINOTAUR -- It's only 1/2 bull
// 
//    (C)opyright 2009 - 2014 The MINOTAUR Team.
// 

#include <cmath>

#include "MinotaurConfig.h"
#include "OperationsUT.h"
#include "Operations.h"
#include "Problem.h"

CPPUNIT_TEST_SUITE_REGISTRATION(OperationsTest);
CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(OperationsTest, "OperationsUT");

using namespace Minotaur;

void OperationsTest::testGcd()
{
  CPPUNIT_ASSERT(fabs(Gcd(7.0,10.5)-3.5)<1e-9);
  CPPUNIT_ASSERT(fabs(Gcd(7.0,-10.5)-3.5)<1e-9);
  CPPUNIT_ASSERT(fabs(Gcd(-7.0,-10.5)-3.5)<1e-9);
  CPPUNIT_ASSERT(fabs(Gcd(-7.0,10.5)-3.5)<1e-9);
  CPPUNIT_ASSERT(fabs(Gcd(-7.0,0.0)-7.0)<1e-9);
  CPPUNIT_ASSERT(fabs(Gcd(0.0,0.0)-0.0)<1e-9);
  CPPUNIT_ASSERT(fabs(Gcd(0.0,-7.0)-7.0)<1e-9);
}


void OperationsTest::testToLower()
{
  std::string s = "MinOTAur MINOtaur !!**&MiNoTaUr~~ -=23";
  std::string l = "minotaur minotaur !!**&minotaur~~ -=23";
  toLowerCase(s);
  CPPUNIT_ASSERT(s==l);
}


void OperationsTest::testSortVarX()
{
  ProblemPtr p = (ProblemPtr) new Problem();
  VarVector vvec;
  double *x = new double[100];

  // test - 1
  for (int i=0; i<100; ++i) {
    vvec.push_back(p->newVariable());
    x[i] = rand() % 100;
  }

  sort(vvec, x);
  for (int i=1; i<100; ++i) {
    CPPUNIT_ASSERT(x[i-1]<=x[i]);
  }
  
  // test - 2
  for (int i=0; i<100; ++i) {
    x[i] = -1;
  }
  sort(vvec, x);
  for (int i=1; i<100; ++i) {
    CPPUNIT_ASSERT(x[i-1]<=x[i]);
  }
  
  // test - 3
  for (int i=0; i<100; ++i) {
    x[i] = vvec[i]->getIndex();
  }
  sort(vvec, x);
  for (UInt i=1; i<100; ++i) {
    CPPUNIT_ASSERT(x[i-1]<=x[i]);
    CPPUNIT_ASSERT(vvec[i]->getIndex()==i);
  }

  delete [] x;
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
