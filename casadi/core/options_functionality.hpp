/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010-2014 Joel Andersson, Joris Gillis, Moritz Diehl,
 *                            K.U. Leuven. All rights reserved.
 *    Copyright (C) 2011-2014 Greg Horn
 *
 *    CasADi is free software; you can redistribute it and/or
 *    modify it under the terms of the GNU Lesser General Public
 *    License as published by the Free Software Foundation; either
 *    version 3 of the License, or (at your option) any later version.
 *
 *    CasADi is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *    Lesser General Public License for more details.
 *
 *    You should have received a copy of the GNU Lesser General Public
 *    License along with CasADi; if not, write to the Free Software
 *    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 *
 */


#ifndef CASADI_OPTIONS_FUNCTIONALITY_HPP
#define CASADI_OPTIONS_FUNCTIONALITY_HPP

#include "generic_type.hpp"
#include "shared_object.hpp"
#include <map>

namespace casadi {

  /// C++ version of Python's dictionary
  typedef GenericType::Dictionary Dictionary;

  // Forward declaration
  class OptionsFunctionalityNode;

/** \brief Provides options setting/getting functionality

  Gives a derived class the ability to set and retrieve options in a convenient way.
  It also contains error checking, making sure that the option exists
  and that the value type is correct.

  A derived class should add option names, types and default values to the corresponding vectors.


  \author Joel Andersson
  \date 2010
  Joel Andersson, K.U. Leuven 2010
  joel.andersson@esat.kuleuven.be
*/
class CASADI_CORE_EXPORT OptionsFunctionality : public SharedObject {
  public:
    /// Default constructor
    OptionsFunctionality();

    /// Destructor
    ~OptionsFunctionality();

    /// Access a member function or object
    OptionsFunctionalityNode* operator->();

    /// Const access a member function or object
    const OptionsFunctionalityNode* operator->() const;

/// \name Option Functionality
/// @{

    /** \brief  set an option.
    For a list of options, check the class documentation of this class.

    The setOptions are only considered before the init function.
    If properties changes, the init function should be called again.
    */
    void setOption(const std::string &str, const GenericType& val);

    /** \brief  set a set of options.
    For a list of options, check the class documentation of this class.

    The setOptions are only considered before the init function.
    If properties changes, the init function should be called again.
    */
    void setOption(const Dictionary& dict, bool skipUnknown = false);

    /** \brief  get an option value */
    GenericType getOption(const std::string &str) const;

    /** \brief  check if there is an option str */
    bool hasOption(const std::string &str) const;

    /** \brief  check if the user has there is an option str */
    bool hasSetOption(const std::string &str) const;

    /** \brief  Print options to a stream */
    void printOptions(std::ostream &stream=std::cout) const;

    /** \brief  Copy all options from another object*/
    void copyOptions(const OptionsFunctionality& obj, bool skipUnknown = false);

    /** \brief  Get the dictionary */
    const Dictionary& dictionary() const;

/// @}
    /** \brief Get a list of all option names */
    std::vector<std::string> getOptionNames() const;

    /** \brief Get the description of a certain option */
    std::string getOptionDescription(const std::string &str) const;

    /** \brief Get the type of a certain option */
    opt_type getOptionType(const std::string &str) const;

    /** \brief Get the type name of a certain option */
    std::string getOptionTypeName(const std::string &str) const;

    /** \brief Get the allowed values of a certain option */
    std::vector<GenericType> getOptionAllowed(const std::string &str) const;

    /// \cond INTERNAL
    /** \brief Get the index into allowed options of a certain option */
    int getOptionAllowedIndex(const std::string &name) const;

    /** \brief Set a certain option by giving its index into the allowed values */
    void setOptionByAllowedIndex(const std::string &name, int i);

    /** \brief Get the enum value corresponding to th certain option */
    int getOptionEnumValue(const std::string &name) const;

    /** \brief Set a certain option by giving an enum value */
    void setOptionByEnumValue(const std::string &name, int v);
    /// \endcond

    /** \brief Get the default of a certain option */
    GenericType getOptionDefault(const std::string &str) const;

    /// Check if a particular cast is allowed
    static bool testCast(const SharedObjectNode* ptr);
};

/// \cond INTERNAL
#ifndef SWIG

/** \brief Internal class
  \author Joel Andersson
  \date 2010
*/
class CASADI_CORE_EXPORT OptionsFunctionalityNode : public SharedObjectNode {
  friend class OptionsFunctionality;
  public:

/// Constructor, destructor
OptionsFunctionalityNode();
virtual ~OptionsFunctionalityNode();

  /** \brief  set an option.
  The setOptions are in general only considered before the init function, if any.
  If properties changes, the init function should be called again.
  (Ticket #54)
  */
  void setOption(const std::string &str, const GenericType& val);

  /** \brief  set a set of options.
  The setOptions are in general only considered before the init function, if any.
  If properties changes, the init function should be called again.
  (Ticket #54)
  */
  void setOption(const Dictionary& dict, bool skipUnknown = false);

  /** \brief Get a list of all option names */
  std::vector<std::string> getOptionNames() const;

  /** \brief Get the description of a certain option */
  std::string getOptionDescription(const std::string &str) const;

  /** \brief Get the type of a certain option */
  opt_type getOptionType(const std::string &str) const;

  /** \brief Get the type name of a certain option */
  std::string getOptionTypeName(const std::string &str) const;

  /** \brief Get the default of a certain option */
  GenericType getOptionDefault(const std::string &str) const;

  #ifndef SWIG
  /** \brief Get the allowed values of a certain option */
  std::vector<GenericType> getOptionAllowed(const std::string &str) const;
  #endif // SWIG

  /** \brief Get the index into allowed options of a certain option */
  int getOptionAllowedIndex(const std::string &name) const;

  /** \brief Set a certain option by giving its index into the allowed values */
  void setOptionByAllowedIndex(const std::string &name, int i);

  /** \brief Get the enum value corresponding to th certain option */
  int getOptionEnumValue(const std::string &name) const;

  /** \brief Set a certain option by giving an enum value */
  void setOptionByEnumValue(const std::string &name, int v);

  /** \brief  check if there is an option str */
  bool hasOption(const std::string &str) const;

  /** \brief  check if the user has there is an option str */
  bool hasSetOption(const std::string &str) const;

  /** \brief  Print options to a stream */
  void printOptions(std::ostream &stream=std::cout) const;

  /** \brief  Print all information there is to know about a certain option */
  void printOption(const std::string &name, std::ostream &stream = std::cout) const;

  /** \brief  get an option value */
  GenericType getOption(const std::string &str) const;

  /** \brief  Print description */
  virtual void print(std::ostream &stream) const = 0;

  /** \brief  Print representation */
  virtual void repr(std::ostream &stream) const;

  /** \brief  Copy all options from another object*/
  void copyOptions(const OptionsFunctionality& obj, bool skipUnknown = false);

  /** \brief  Get the dictionary */
  const Dictionary& dictionary() const;

  /** \brief Get the best suggestions for a misspelled word using a dictionary
  *
  *  \param[in] word  The word that is misspelled
  *  \param[in] dictionary  A list of correct words
  *  \param[out] suggestions The list of suggestions generated. This list will be cleared and filled.
  *  \param[in] amount Maximum number of suggestions
  *
  * \return Some metric for confidence about the best match
  */
  static double getBestMatches(const std::string & word,
                               const std::vector<std::string> &dictionary,
                               std::vector<std::string> &suggestions, int amount = 5);

  /** \brief Get th ebest suggestions of option names
  */
  double getBestMatches(const std::string & name, std::vector<std::string> &suggestions,
                        int amount = 5) const;


  /** \brief A distance metric between two words */
  static double wordDistance(const std::string &a, const std::string &b);


  void addOption(
    const std::string &str, const opt_type& type,
    const GenericType &def_val=GenericType(), const std::string& desc="n/a",
    const std::vector<GenericType> &allowed_vals = std::vector<GenericType>(),
    bool inherit = false, std::vector<int> enum_values= std::vector<int>(),
    std::vector<std::string> enum_descr= std::vector<std::string>());
  /** \brief Add an option
  *
  *  allowed_vals can take multiple forms:
  *    "foo|bar"   ->   specifies that the values "foo" and "bar" are allowed
  *    "foo:5|bar:6" -> specifies that the values "foo" and "bar" are allowed and map
  *                     to 5 and 6 respectively
  *    "foo:5:description_foo|bar:6:description_bar|" -> same as above, but specifies documentation
  *
  **/
  void addOption(const std::string &str, const opt_type& type, const GenericType &def_val,
                 const std::string& desc, const std::string &allowed_vals, bool inherit = false);

protected:

void assert_exists(const std::string &str) const;


/** \brief Sets the default value for an option without changing the current value
*/
void setDefault(const std::string &str, const GenericType &def_val);

private:

/** \brief  Allowed options  */
  std::map<std::string, opt_type> allowed_options;

/** \brief  User-set options */
  Dictionary dictionary_;

/** \brief  Option defaults */
  Dictionary defaults_;

/** \brief  Description for the options */
  std::map<std::string, std::string> description_;

/** \brief  Allowed values for the options */
  std::map<std::string, std::vector<GenericType> > allowed_vals_;

/** \brief  Enum values */
  std::map<std::string, std::vector<int> > enum_values_;

/** \brief  Enum descriptions */
  std::map<std::string, std::vector<std::string> > enum_descr_;

};

#endif // SWIG

/// \endcond

} // namespace casadi


#endif // OPTIONS_FUNCTIONALITY
