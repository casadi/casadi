//
// MINOTAUR -- It's only half bull!
//
// (C)opyright 2009 - 2014 The MINOTAUR Team.
//

/**
 * \file Environment.h
 * \brief Define the Environment class.
 * \author Ashutosh Mahajan, Argonne National Laboratory
 */

#ifndef ENVIRONMENT_H
#define ENVIRONMENT_H

#include "Types.h"

namespace Minotaur {

  class Interrupt;
  class Timer;
  class TimerFactory;

  /**
   * The environment is a container class that has pointers to the
   * Logger, interrupt handler, timer factory and options that need to be
   * passed to Minotaur.
   */
  class Environment {
    public:
      /// Default constructor.
      Environment();

      /// Destroy.
      ~Environment();

      /// Get the logger that is setup for the whole environment.
      LoggerPtr getLogger() const;

      /// Get the current log level
      LogLevel getLogLevel() const;

      /// Get a new timer. The calling function has to free this timer.
      Timer *getNewTimer();

      /// Get the options database.
      OptionDBPtr getOptions();

      /**
       * Get the time from the 'global timer' i.e. the total time consumed so
       * far.
       *
       * \param[out] err 0 if no error occured, positive otherwise.
       * \return The time used so far.
       */
      double getTime(int &err);

      /**
       * \brief Get the global timer.
       *
       * \return The timer used to query the time elapsed since the solver
       * started.
       */
      const Timer* getTimer();

      /// Get the version string
      std::string getVersion();

      /// Read the options using a char array that has 'argc' words in it.
      void readOptions(int argc, char **argv);

      /// Read options from a string.
      void readOptions(const std::string &str);

      /// Set the log level of the default logger
      void setLogLevel(LogLevel l);

      /**
       * \brief Start a 'global' timer that can be queried for the total time used
       * so far by the solution process.
       *
       * \param[out] err 0 if no error occured, positive otherwise.
       */
      void startTimer(int &err);
      
      /**
       * \brief Stop the 'global' timer.
       *
       * \param[out] err 0 if no error occured, positive otherwise.
       */
      void stopTimer(int &err);

      /** 
       * Write the full information about the version number and how Minotaur
       * was compiled. TODO: Fix it.
       */
      void writeFullVersion(std::ostream &out);

    private:
      /// Logger that is used for the whole environment.
      LoggerPtr logger_;

      /// For log:
      static const std::string me_;

      /// The options database
      OptionDBPtr options_;

      /// The global timer
      Timer *timer_;

      /// The generator that is used to build timers.
      TimerFactory *timerFac_;

      /// Add an option to the database.
      void convertAndAddOneOption_(BoolOptionPtr &b_option,
                                   IntOptionPtr &i_option,
                                   DoubleOptionPtr &d_option,
                                   StringOptionPtr &s_option,
                                   FlagOptionPtr &f_option,
                                   std::string &name,
                                   std::string &value,
                                   std::ostringstream &logstr);
        
      /// Create default set of options
      void createDefaultOptions_();

      /*
       * Find an option with the provided name. The option is option of that
       * type is returned and all others are NULLed.
       */
      void findOption_(const std::string &name, BoolOptionPtr &b_option,
          IntOptionPtr &i_option, DoubleOptionPtr &d_option, 
          StringOptionPtr &s_option, FlagOptionPtr &f_option);

      /**
       * Return true if the string "str" corresponds to true. ("T", "TRUE",
       * "Y", "YES", "1", or their lower case equivalents). Returns false
       * otherwise.
       */
      bool getBoolValue_(const std::string & str);

      /**
       * Reach options from config file and increment the number of
       * instances.
       */
      void readConfigFile_(std::string fname, UInt &num_p);

      /**
       * Remove leading dashes "--" or "-" from front of a string. Returns
       * the number of dashes removed. The string "name" is modified.
       */
      UInt removeDashes_(std::string &name);

      /// Remove the string "minotaur." from front of a string.
      void removeMinotaurPrepend_(std::string & fullname);

      /**
       * Remove the substring that follows '=' sign in 'name' and return the
       * substring.
       */
      std::string separateEqualToArg_(std::string &name);
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
