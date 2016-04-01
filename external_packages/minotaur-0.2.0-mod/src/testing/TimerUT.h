//
//    MINOTAUR -- It's only 1/2 bull
//
//    (C)opyright 2009 - 2014 The MINOTAUR Team.
//

#ifndef TIMERUT_H
#define TIMERUT_H

#include <cppunit/TestCase.h>
#include <cppunit/TestCaller.h>
#include <cppunit/TestSuite.h>
#include <cppunit/TestResult.h>
#include <cppunit/extensions/HelperMacros.h>

using namespace Minotaur;

class TimerUT : public CppUnit::TestCase {
  public:
    TimerUT(std::string name) : TestCase(name) {}
    TimerUT() {}

    void testSleep();
    void setUp();
    void tearDown();

    CPPUNIT_TEST_SUITE(TimerUT);
    CPPUNIT_TEST(testSleep);
    CPPUNIT_TEST_SUITE_END();

  private:
    TimerFactory *tFactory_;
};

#endif     // #define TIMERUT_H

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
