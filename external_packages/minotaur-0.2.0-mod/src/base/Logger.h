// 
//     MINOTAUR -- It's only 1/2 bull
// 
//     (C)opyright 2010 - 2014 The MINOTAUR Team.
// 

/**
 * \file Logger.h
 * \brief Declare class for creating logs.
 * \author Ashutosh Mahajan, Argonne National Laboratory
 */

#ifndef MINOTAURLOGGER_H
#define MINOTAURLOGGER_H

/**
 * This class handles all output.  Attached to each message is the level and
 * only output requested by the user is generated.  Formation of C-style 
 * output is delayed until after the level is checked.  C++-style output
 * requires that a string be written and passed to the logger; the code
 * needs to check if the output would be used prior to constructing it.
 * 
 * Note GAMS will create a derived class that differentiates where to
 * send the output (log, listing, or status files) based on the level.
 * A FORTRAN derived class will create a string and send it to a wrapper 
 * that writes the string using FORTRAN calls.
 * 
 * We may want to describe the output levels in this section.
 */

#include <iostream>
#include <streambuf>

#include "Types.h"

namespace Minotaur {
  class Logger {
    public:
      /// Default constructor
      Logger(LogLevel max_level=LogInfo);

      /// Destroy
      virtual ~Logger();

      /// Do not write messages that are above maxLevel
      inline void setMaxLevel(LogLevel maxLevel) { maxLevel_ = maxLevel; }

      /// Get the maxLevel
      inline LogLevel getMaxLevel() const { return maxLevel_; }

      /// Get the stream where one can write messages.
      virtual std::ostream& msgStream(LogLevel level) const;

      /// Get the stream where one can write errors.
      std::ostream& errStream() const;

    protected:
      // Maximum output level
      LogLevel maxLevel_;

    private:
      /// Copy constructor is not allowed.
      Logger (const Logger &l);

      /// Copy by assignment is not allowed.
      Logger  & operator = (const Logger &l);

      /// nullbuf class does not print anything.
      class nullbuf : public std::streambuf {
        protected:
          /// \todo Don't know what this does.
          virtual int_type overflow(int_type c) { return c; }
      };
      /// Null stream buffer
      nullbuf nb_;

      // Null output stream
      mutable std::ostream nout_;
  };

  typedef boost::shared_ptr<Logger> LoggerPtr;
  typedef boost::shared_ptr<const Logger> ConstLoggerPtr;
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
