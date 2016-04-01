//
//     MINOTAUR -- It's only 1/2 bull
//
//     (C)opyright 2008 - 2014 The MINOTAUR Team.
//

#ifndef MINOTAURINTERRUPT_H
#define MINOTAURINTERRUPT_H

// This class is used to detect interrupts and exit gracefully.
// Exception CANNOT be thrown within the signal handler, which 
// is the reason for having a separate function to check if an
// interrupt has occurred.

// Questions:
//   Can we eliminate the static member functions?
//   Should we use setjmp and longjmp when an interrupt is triggered?
//   What should exception should be thrown on an interrupt?
//   If the keyboard interrupt is obtained after some limit, we exit
//   from the code.  However, we would like to use the Minotaur::Logger
//   to log the reason.  How is this done without a static pointer to
//   a Logger?


#include "Exception.h"
#include "Logger.h"

namespace Minotaur {

  class InterruptException : public Exception {
    private:
      InterruptException& operator=(const InterruptException&);

    public:
      InterruptException() : Exception() { }
      InterruptException(const InterruptException&) : Exception() { }
      ~InterruptException() { }
  };

  class Interrupt {
    public:
      Interrupt(UInt l=5) {limit_ = l; }
      virtual ~Interrupt() {}			// Destructor

      inline void SetLimit(UInt l) { limit_ = l; }
      inline UInt GetLimit() const { return limit_; }

      inline void SetLogger(LoggerPtr l) { logger_ = l; }
      inline LoggerPtr GetLogger() const { return logger_; }

      virtual void Start();
      virtual void Stop();
      virtual void Check() const throw(InterruptException);

    private:
      void (*originalHandler_)(int);

      static volatile UInt interrupts_;
      static UInt limit_;		// Maximum 
      static LoggerPtr logger_;

      static void Handler(int);

      Interrupt (const Interrupt &);
      Interrupt & operator = (const Interrupt &);
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
