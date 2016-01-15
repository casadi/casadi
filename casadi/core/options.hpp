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


#ifndef CASADI_OPTIONS_HPP
#define CASADI_OPTIONS_HPP

#include "generic_type.hpp"

namespace casadi {
/// \cond INTERNAL
#ifndef SWIG

  /** \brief Options metadata for a class
      \author Joel Andersson, Joris Gillis
      \date 2010-2016
  */
  struct CASADI_EXPORT Options {
    // Base classes, whose options are also valid for the derived class
    std::vector<Options*> bases;

    // Information for a particular options entry
    struct Entry {
      TypeID type;
      std::string description;

      // Print entry
      void print(const std::string& name, std::ostream &stream) const;
    };

    // Lookup for options
    std::map<std::string, Entry> entries;

    // Locate an entry
    const Options::Entry* find(const std::string& name) const;

    // Print all entries
    void print(std::ostream &stream) const;

    /** \brief A distance metric between two words */
    static double word_distance(const std::string &a, const std::string &b);

    /** \brief Get the best suggestions for a misspelled word */
    std::vector<std::string> suggestions(const std::string& word, int amount=5) const;

    /** \brief Find best matches */
    void best_matches(const std::string& word,
                      std::vector<std::pair<double, std::string> >& best) const;
  };

#endif // SWIG
  /// \endcond
} // namespace casadi


#endif // CASADI_OPTIONS_HPP
