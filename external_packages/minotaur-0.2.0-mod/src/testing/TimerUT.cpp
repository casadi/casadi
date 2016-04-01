//
//    MINOTAUR -- It's only 1/2 bull
//
//    (C)opyright 2009 - 2014 The MINOTAUR Team.
//

#include "MinotaurConfig.h"
#include "Timer.h"
#include "TimerUT.h"


CPPUNIT_TEST_SUITE_REGISTRATION(TimerUT);
CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(TimerUT, "TimerUT");

using namespace Minotaur;

void TimerUT::setUp()
{
  tFactory_ = new TimerFactory();     
}

void TimerUT::tearDown()
{
  delete tFactory_;
}

void TimerUT::testSleep()
{
  double time_used;
  Timer *timer=tFactory_->getTimer();

  // start timer
  timer->start();
  // and query immediately
  time_used = timer->query();
  CPPUNIT_ASSERT(time_used <= 0.001);

  // also query after 1 seconds
  while (timer->query() < 1) {
  }
  time_used = timer->query();
  CPPUNIT_ASSERT(time_used >= 1.0);
  CPPUNIT_ASSERT(time_used <= 2.0);
  timer->stop();

  // finally sleep test
  timer->start();
  sleep(2);
  time_used = timer->query();
  CPPUNIT_ASSERT(time_used <= 1.5);
  delete timer;
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
