// 
//     MINOTAUR -- It's only 1/2 bull
// 
//     (C)opyright 2009 - 2014 The MINOTAUR Team.
// 

#ifndef LINEARFUNCTIONUT_H
#define LINEARFUNCTIONUT_H

#include <string>

#include <cppunit/TestCase.h>
#include <cppunit/TestCaller.h>
#include <cppunit/TestSuite.h>
#include <cppunit/TestResult.h>
#include <cppunit/extensions/HelperMacros.h>

#include <Problem.h>
#include <Function.h>

using namespace Minotaur;

class ProblemTest : public CppUnit::TestCase {

  public:
    ProblemTest(std::string name) : TestCase(name) {}
    ProblemTest() {}

    void setUp();
    void tearDown();
    void testevalCon();
    void testConsTypes();
    void testObjTypes();  
    void testVarTypes(); 
    void testDeleteVar(); 
    void testChangeBound(); 
    void testaddToObj(); 
 
    CPPUNIT_TEST_SUITE(ProblemTest);
    CPPUNIT_TEST(testevalCon);
    CPPUNIT_TEST(testConsTypes); 
    //  CPPUNIT_TEST(testgetCons);
    CPPUNIT_TEST(testObjTypes); 
    CPPUNIT_TEST(testVarTypes);
    CPPUNIT_TEST(testDeleteVar);
    CPPUNIT_TEST(testChangeBound); 
    CPPUNIT_TEST(testaddToObj);  
    CPPUNIT_TEST_SUITE_END();

    //void testgetCons();

  private:
    ProblemPtr instance_;

};

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
