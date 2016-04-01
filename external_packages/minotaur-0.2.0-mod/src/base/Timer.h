//
//    MINOTAUR -- It's only 1/2 bull
//
//    (C)opyright 2009 - 2014 The MINOTAUR Team.
//

/**
 * \file Timer.h
 * \brief Define the TimerFactory and varios Timer classes
 * \author Todd Munson, Argonne National Laboratory
 */


#ifndef MINOTAURTIMER_H
#define MINOTAURTIMER_H

#include <time.h>

#ifndef CLOCKS_PER_SEC
#  define CLOCKS_PER_SEC 1000000
#endif

#ifdef MINOTAUR_RUSAGE
#  include <sys/resource.h>

#  ifndef MICROSEC
#    define MICROSEC 1000000
#  endif
#endif

#include "Types.h"

namespace Minotaur {
  /**
   * The Timer class is used to measure time. A Timer can be queried only
   * after starting and before stopping it. A Timer should always be
   * generated through the TimerFactory::GetTimer() method because different
   * machines may need different types of timers.
   */
  class Timer {
  public:
    Timer() { };
    virtual ~Timer() { };

    virtual void start()   = 0;
    virtual void stop()    = 0;
    virtual double query() const = 0;

  private: 
    Timer (const Timer &);
    Timer & operator = (const Timer &);
  };

  /**
   * The default timer when getruage() function is not defined. Uses clock()
   * to get time.
   */
  class ClockTimer : public Timer {
  private:
    clock_t s_;
    bool is_started_;

  public:
    ClockTimer() : is_started_(false)  {  };
    ~ClockTimer() { };

    /// Start the timer.
    void start() {
      s_ = clock();
      is_started_ = true;      
      return;
    };

    /// Stop the timer. Can not query after this.
    void stop() {
      is_started_ = false;
      return;
    };

    double query() const {
      if (!is_started_) {
        throw("Some exception");
      }
      
      return ((clock() - s_) / (double) CLOCKS_PER_SEC);
    };
  };

#ifdef MINOTAUR_RUSAGE
  /**
   * The default timer when getruage() function is avaialble. Uses rusage
   * structure to get time.
   */
  class UsageTimer : public Timer {
  private:
    struct rusage s_;
    bool is_started_;

  public:
    UsageTimer() : is_started_(false)  { };
    ~UsageTimer() { };

    /// Start the timer.
    void start() {
      getrusage(RUSAGE_SELF, &s_);
      is_started_ = true;      
      return;
    };

    /// Stop the timer. Can not query after this.
    void stop() {
      is_started_ = false;
      return;
    };

    /// Get how much cpu time has been used since this timer was started.
    double query() const {
      double t1, t2, m1, m2, s1, s2, m, s;
      struct rusage e;
      if (!is_started_) {
        throw("Some exception");
      }
      getrusage(RUSAGE_SELF, &e);

      m1 = (double) e.ru_utime.tv_usec;
      m2 = (double) s_.ru_utime.tv_usec;
      m = m1 - m2;

      s1 = (double) e.ru_utime.tv_sec;
      s2 = (double) s_.ru_utime.tv_sec;
      s = s1 - s2;

      t1 = s + m / MICROSEC;
    
      m1 = (double) e.ru_stime.tv_usec;
      m2 = (double) s_.ru_stime.tv_usec;
      m = m1 - m2;

      s1 = (double) e.ru_stime.tv_sec;
      s2 = (double) s_.ru_stime.tv_sec;
      s = s1 - s2;

      t2 = s + m / MICROSEC;

      return t1 + t2;
    };
  };
#endif

  /// The TimerFactory should be used to get the approrpriate Timer.
  class TimerFactory {
  public:
    /// Default constructor.
    TimerFactory() { };

    /// Destroy.
    virtual ~TimerFactory() { };

    /// Return an appropriate Timer.
    virtual Timer *getTimer() {
#ifdef MINOTAUR_RUSAGE
      return new UsageTimer;
#else
      return new ClockTimer;
#endif
    };

  private: 
    TimerFactory (const TimerFactory &);
    TimerFactory & operator = (const TimerFactory &);
  };
}

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
