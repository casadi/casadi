// 
//     MINOTAUR -- It's only 1/2 bull
// 
//     (C)opyright 2009 - 2014 The MINOTAUR Team.
// 

#include "MinotaurConfig.h"
#include "PolyUT.h"

#include "Function.h"
#include "LinearFunction.h"
#include "PolynomialFunction.h"
#include "Types.h"
#include "Variable.h"

CPPUNIT_TEST_SUITE_REGISTRATION(PolyUT);
CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(PolyUT, "PolyUT");

using namespace Minotaur;

void PolyUT::monomial()
{
  VariablePtr x1 = (VariablePtr) new Variable(0, 0, 0., 1., Continuous, "x1");
  VariablePtr x2 = (VariablePtr) new Variable(1, 1, 0., 1., Continuous, "x2");
  MonomialFunPtr m = (MonomialFunPtr) new MonomialFunction(2., x1, 4); 
  int error=0;
  m->multiply( 5., x2, 2);
  m->multiply(-5., x2, 1);
  //m->write(std::cout);

  // -50 * x1^4 * x2^3
  double x[2];
  x[0] = 0.;
  x[1] = 1.;
  CPPUNIT_ASSERT(0 == m->eval(x, &error));
  CPPUNIT_ASSERT(0==error);

  x[0] = 1.;
  x[1] = 0.;
  CPPUNIT_ASSERT(0 == m->eval(x, &error));
  CPPUNIT_ASSERT(0==error);

  x[0] = 1.;
  x[1] = 1.;
  CPPUNIT_ASSERT(-50 == m->eval(x, &error));
  CPPUNIT_ASSERT(0==error);

  x[0] = 1.;
  x[1] = 2.;
  CPPUNIT_ASSERT(-400 == m->eval(x, &error));
  CPPUNIT_ASSERT(0==error);


  // gradient
  double g[2];
  g[0] = g[1] = 0.;
  x[0] = 0.;
  x[1] = 1.;
  m->evalGradient(x, g, &error);
  CPPUNIT_ASSERT(0 == g[0]);
  CPPUNIT_ASSERT(0 == g[1]);
  CPPUNIT_ASSERT(0==error);
  x[0] = 1.;
  x[1] = 1.;
  g[0] = g[1] = 0.;
  m->evalGradient(x, g, &error);
  //std::cout << g[0] <<  " " <<  g[1]<< std::endl;
  CPPUNIT_ASSERT(-200 == g[0]);
  CPPUNIT_ASSERT(-150 == g[1]);
  CPPUNIT_ASSERT(0==error);
  x[0] = 1.;
  x[1] = 2.;
  g[0] = g[1] = 0.;
  m->evalGradient(x, g, &error);
  CPPUNIT_ASSERT(-1600 == g[0]);
  CPPUNIT_ASSERT( -600 == g[1]);
  CPPUNIT_ASSERT(0==error);
}


void PolyUT::polynomial()
{
  VariablePtr x1 = (VariablePtr) new Variable(0, 0, 0., 1., Continuous, "x1");
  VariablePtr x2 = (VariablePtr) new Variable(1, 1, 0., 1., Continuous, "x2");
  PolyFunPtr  p1 = (PolyFunPtr) new PolynomialFunction();
  MonomialFunPtr  m;
  LinearFunctionPtr lf = (LinearFunctionPtr) new LinearFunction();
  NonlinearFunctionPtr nlf;
  QuadraticFunctionPtr qf = QuadraticFunctionPtr();
  FunctionPtr f;
  int error = 0;

  lf->addTerm(x1, 1.0);

  // x1^3 + 3x1x2
  m = (MonomialFunPtr) new MonomialFunction(1, x1, 3);
  (*p1) += m;

  m = (MonomialFunPtr) new MonomialFunction(3, x1, 1);
  m->multiply(1., x2, 1);
  (*p1) += m;

  //p1->write(std::cout);
  double x[2];
  x[0] = 0.;
  x[1] = 0.;
  CPPUNIT_ASSERT(0. == p1->eval(x, &error));
  CPPUNIT_ASSERT(0==error);

  x[0] =-1.;
  x[1] = 1.;
  CPPUNIT_ASSERT(-4. == p1->eval(x, &error));
  CPPUNIT_ASSERT(0==error);

  x[0] = 2.;
  x[1] = 0.;
  CPPUNIT_ASSERT(8. == p1->eval(x, &error));
  CPPUNIT_ASSERT(0==error);

  double g[2];
  x[0] = 0.;
  x[1] = 1.;
  g[0] = g[1] = 0.;
  p1->evalGradient(x, g, &error);
  CPPUNIT_ASSERT(3. == g[0]);
  CPPUNIT_ASSERT(0. == g[1]);
  CPPUNIT_ASSERT(0==error);

  x[0] = 2.;
  x[1] = 1.;
  g[0] = g[1] = 0.;

  f = (FunctionPtr) new Function(p1);
  CPPUNIT_ASSERT(Polynomial == f->getType());
  lf = LinearFunctionPtr(); // NULL
  f = (FunctionPtr) new Function(lf, qf, p1);
  CPPUNIT_ASSERT(Polynomial == f->getType());

  p1->evalGradient(x, g, &error);
  CPPUNIT_ASSERT(0==error);
  CPPUNIT_ASSERT(15. == g[0]);
  CPPUNIT_ASSERT( 6. == g[1]);

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
