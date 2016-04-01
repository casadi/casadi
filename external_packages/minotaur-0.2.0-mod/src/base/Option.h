// 
//     MINOTAUR -- It's only 1/2 bull
// 
//     (C)opyright 2008 - 2014 The MINOTAUR Team.
// 

/**
 * \file Option.h
 * \brief Declare the Option class for managing options for Minotaur.
 * \author Todd Munson, Argonne National Laboratory
 * \author Ashutosh Mahajan, Argonne National Laboratory
 */

#ifndef MINOTAUR_OPTIONS
#define MINOTAUR_OPTIONS

#include <string>

#include "Types.h"

namespace Minotaur {
  /**
   * Options are used to control the behavior of Minotaur: choice of
   * algorithms, parameters, limits, outputs etc. They are also used to
   * control the behavior of interfaces and engines that are linked with
   * Minotaur: parameters for LP engine, verbosity, limits etc. This templated
   * class works for any option whose value could be integer, double or string.
   * 
   * All options in Minotaur must be one of the above three types. In the
   * command line, all option names must be entered with a preceding '--'.
   * e.g.
   * ./minotaur --assume_convex yes --node_limit=2 --time_limit 100
   * 
   * Either form '--option_name=option_value' or '--option_name option value'
   * is acceptable.
   * 
   * Any simple options without values, like for instance '-v', '-q' are
   * not allowed, but they can still be accepted if they are aliased to a
   * minotaur acceptable option. In such cases, we need to tell Minotaur, how
   * to expand small options like '-=', '-v', '-AMPL' into a proper Minotaur
   * option.  ... more to be added here.
   * 
   * Interfaces to Minotaur, like AMPL, can also add their own options to the
   * Minotaur option-database.
   * 
   * Again, we assume that all options that begin with '--' have a value
   * associated with it. We also assume that all options that begin with '-'
   * do not have any value associated with it. So, for the command,
   * ./minotaur -v 2 --assume_convex yes stub.nl
   * we will assume that '2' is the name of input file, 'yes' is the value of
   * option 'assume_convex' and stub.nl is also an input file.
   */
  template <class T> class Option {
  public:
    /// Construct the option using name, description and value.
    Option(const std::string& name, const std::string& desc, 
           bool is_known=false, T val=0);

    /// Destroy.
    virtual ~Option();

    /// Set the value of this option.
    virtual void setValue(T val) { val_ = val; return; };

    /// Get the value of option.
    virtual T getValue() { return val_; };

    /// Get the name of this option.
    virtual const std::string & getName();

    /// Get the help description of this option.
    virtual const std::string & getDesc();

    /// Check if this option was used somewhere.
    virtual bool wasEverUsed() {return everUsed_;};

    /// Set the 'used' flag. Used flag is true if
    virtual void setUsedFlag(const bool &used_flag) {everUsed_ = used_flag;};

    /// Return true if the option is known to Minotaur.
    virtual bool isKnown() {return isKnown_;};

    /// Set or unset the 'known' flag
    virtual void setKnownFlag(const bool &known_flag) 
    {isKnown_ = known_flag;};

    /// Write to the output stream.
    virtual void write(std::ostream &out) const;

  protected:
    /**
     * The name of the option should not have any spaces. All spaces should
     * be converted to underscore (_). If the name consists of a period (.),
     * then the option is meant for the external library referred to by the
     * name before the period. For instance the name:
     * ipopt.some_tolerance
     * tells us to provide Ipopt with the option "some_tolerance". Again,
     * spaces are not allowed even if the external library has spaces in the
     * name.
     */
    std::string name_;

    /**
     * A human-readable description of what the option does and what limits
     * are expected.
     */
    std::string desc_;

    /**
     * The value of the option. For now, we don't check if the user
     * specified options lie in the suggested limits.
     */
    T val_;

    /**
     * True, if this option was ever used by any of the procedures. False
     * otherwise.
     */
    bool everUsed_;

    /// True, if the option is known to minotaur.
    bool isKnown_;

  private:
    /// Copying is not allowed.
    Option (const Option<T> &o);

    /// Copy by assignment is not allowed.
    Option  & operator = (const Option<T> &o);
  };


  /**
   * OptionDB is a class to store a set of options. One use of this class 
   * is to save the user specified options.
   * The user may specify options using the command line arguments, or through
   * the API. Further, some of the options may be invalid (with typos). This
   * class can tell if the options specified by the user are legitimate
   * options.
   */
  class OptionDB {
  public:
    /// Default constructor.
    OptionDB();

    /// Destroy.
    ~OptionDB();

    /// Add a bool option to the database.
    void insert(BoolOptionPtr option, bool is_flag=false);

    /// Add an int option to the database.
    void insert(IntOptionPtr option);

    /// Add a double option to the database.
    void insert(DoubleOptionPtr option);

    /// Add a string option to the database.
    void insert(StringOptionPtr option);

    /// Find a bool option in the database.
    BoolOptionPtr findBool(const std::string &name);

    /// Find an int option in the database.
    IntOptionPtr findInt(const std::string &name);

    /// Find a double option in the database.
    DoubleOptionPtr findDouble(const std::string &name);

    /// Find a string option in the database.
    StringOptionPtr findString(const std::string &name);

    /// Find a flag option in the database.
    FlagOptionPtr findFlag(const std::string &name);

    /**
     * Remove from the database all options that have the provided name.
     * Return the number of options that were removed.
     */
    UInt remove(const std::string &name);

    /// Iterator to access the first bool option.
    BoolOptionSetIter boolBegin();

    /// Iterator to access the last bool option.
    BoolOptionSetIter boolEnd();

    /// Iterator to access the first integer option.
    IntOptionSetIter intBegin();

    /// Iterator to access the last integer option.
    IntOptionSetIter intEnd();

    /// Iterator to access the first double option.
    DoubleOptionSetIter dblBegin();

    /// Iterator to access the last double option.
    DoubleOptionSetIter dblEnd();

    /// Iterator to access the first string option.
    StringOptionSetIter strBegin();

    /// Iterator to access the last string option.
    StringOptionSetIter strEnd();

    /// Iterator to access the first flag.
    FlagOptionSetIter flagBegin();

    /// Iterator to access the last flag.
    FlagOptionSetIter flagEnd();

    /**
     * Write the database  to the output stream. It will print the option
     * name, the value and if it was ever used.
     */
    void write(std::ostream &out) const;

  private:
    /// Set of all boolean options.
    BoolOptionSet bool_ops_;

    /// Set of all integer options.
    IntOptionSet int_ops_;

    /// Set of all double options.
    DoubleOptionSet double_ops_;

    /// Set of all string options.
    StringOptionSet string_ops_;

    /// Set of all flags (options that don't need any arguments).
    FlagOptionSet flag_ops_;
  };
  typedef boost::shared_ptr <OptionDB> OptionDBPtr;
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
