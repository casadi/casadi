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

/// \cond INTERNAL
#ifndef SWIG

  struct CASADI_EXPORT Options {
    // Base classes, whose options are also valid for the derived class
    std::vector<Options*> bases;

    // Information for a particular options entry
    struct Entry {
      TypeID type;
      std::string description;
    };

    // Lookup for options
    std::map<std::string, Entry> entries;

    // Locate an entry
    const Options::Entry* find(const std::string& name) const;
  };

  /** \brief Provides options functionality to a derived class
      \author Joel Andersson
      \date 2010-2015
  */
  class CASADI_EXPORT OptionsFunctionality {
  public:
    /** \brief Get a list of all option names */
    std::vector<std::string> optionNames() const;

    /** \brief  Print options to a stream */
    void printOptions(std::ostream &stream=casadi::userOut()) const;

    /** \brief  Print all information there is to know about a certain option */
    void printOption(const std::string &name, std::ostream &stream = userOut()) const;

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

  protected:

    void assert_exists(const std::string &str) const;

  private:

    /** \brief  Allowed options  */
    std::map<std::string, TypeID> allowed_options;

    /** \brief  Description for the options */
    std::map<std::string, std::string> description_;
  };

#endif // SWIG

  /// \endcond

} // namespace casadi


#endif // OPTIONS_FUNCTIONALITY
