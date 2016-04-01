// 
//     MINOTAUR -- It's only 1/2 bull
// 
//     (C)opyright 2010 - 2014 The MINOTAUR Team.
// 

// /**
// \file EnvironmentUT.h
// \brief Header for testing the Environment class
// \author Ashutosh Mahajan, Argonne National Laboratory
// */

#ifndef ENVIRONMENTUT_H
#define ENVIRONMENTUT_H

#include <cppunit/TestCase.h>
#include <cppunit/TestCaller.h>
#include <cppunit/TestSuite.h>
#include <cppunit/TestResult.h>
#include <cppunit/extensions/HelperMacros.h>

#include <Environment.h>

using namespace Minotaur;

class EnvironmentUT : public CppUnit::TestCase {

  public:
    EnvironmentUT(std::string name) : TestCase(name) {}
    EnvironmentUT() {}

    void setUp();
    void tearDown();

    CPPUNIT_TEST_SUITE(EnvironmentUT);
    CPPUNIT_TEST(testLogger);
    CPPUNIT_TEST_SUITE_END();

    void testLogger();

  private:
    EnvPtr env_;

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
