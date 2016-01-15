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


#include "options_functionality.hpp"

#include "std_vector_tools.hpp"
#include "exception.hpp"
#include <algorithm>
#include <string>
#include <locale>

#include "matrix.hpp"

using namespace std;

namespace casadi {

  const Options::Entry* Options::find(const std::string& name) const {
    // Check if in one of the bases
    for (auto&& b : bases) {
      // Call recursively
      const Options::Entry* entry = b->find(name);
      if (entry) return entry;
    }

    // Lookup in this class
    auto it = entries.find(name);
    if (it!=entries.end()) {
      return &it->second;
    } else {
      return 0;
    }
  }

  double OptionsFunctionality::wordDistance(const std::string &a, const std::string &b) {
    /// Levenshtein edit distance
    if (a == b) return 0;
    int na = a.size();
    int nb = b.size();
    if (na == 0) return nb;
    if (nb == 0) return na;

    vector<int> v0(nb+1, 0);
    vector<int> v1(nb+1, 0);

    for (int i=0;i<nb+1;++i)
      v0[i] = i;

    char s;
    char t;
    std::locale loc;
    for (int i=0;i<na;i++) {
      v1[0] = i + 1;
      for (int j=0; j<nb; j++) {
        s = std::tolower(a[i], loc);
        t = std::tolower(b[j], loc);
        int cost = 0;
        if (s != t)
          cost = 1;

        v1[j+1] = min(min(v1[j] + 1, v0[j+1] + 1), v0[j] + cost);
      }

      for (int j=0; j<nb+1; j++)
        v0[j] = v1[j];
    }

    return v1[nb];
  }

  /// \cond INTERNAL
  /// A helper class to use stl::sort in OptionsFunctionality::getBestMatches
  struct mysortclass {
    bool operator()(std::pair<std::string, double> a, std::pair<std::string, double> b) {
      return (a.second<b.second);}
  } mysorter;
  /// \endcond

  double OptionsFunctionality::
  getBestMatches(const std::string & word,
                 const std::vector<std::string> &dictionary,
                 std::vector<std::string> &suggestions, int amount) {
    // Make a list of (word, score) tuples
    std::vector< std::pair<std::string, double> > candidates(dictionary.size());

    // Fill this list
    for (int i=0;i<dictionary.size();i++) {
      candidates[i].first  = dictionary[i];
      candidates[i].second = wordDistance(word, dictionary[i]);
    }

    // Sort it
    sort(candidates.begin(), candidates.end(), mysorter);

    // Put the first 'amount' of them in suggestions
    suggestions.clear();
    for (int i=0;i<amount;i++) {
      if (i<candidates.size()) {
        suggestions.push_back(candidates[i].first);
      }
    }

    return -1; // No score metric yet
  }

  double OptionsFunctionality::getBestMatches(const std::string &name,
                                                  std::vector<std::string> &suggestions,
                                                  int amount) const {
    // Work towards a vector of option names
    std::vector< std::string> dict;

    // Fill it by looping over the allowed_options map
    for (auto it=allowed_options.begin(); it!=allowed_options.end(); it++) {
      dict.push_back(it->first);
    }

    // Pass the work on to the more general method
    return getBestMatches(name, dict, suggestions, amount);
  }

  void OptionsFunctionality::assert_exists(const std::string &name) const {
    // If option contains a dot, only check what comes before
    auto dotpos = name.find('.');
    if (dotpos != string::npos) {
      return assert_exists(name.substr(0, dotpos));
    }

    // First check if the option exists
    map<string, TypeID>::const_iterator it = allowed_options.find(name);
    if (it == allowed_options.end()) {
      stringstream ss;
      ss << "Unknown option: " << name << endl;
      std::vector<std::string> suggestions;
      getBestMatches(name, suggestions, 5);
      ss << endl;
      ss << "Did you mean one of the following?" << endl;
      for (int i=0;i<suggestions.size();++i)
        printOption(suggestions[i], ss);
      ss << "Use printOptions() to get a full list of options." << endl;
      casadi_error(ss.str());
    }
  }

  void OptionsFunctionality::printOption(const std::string &name, ostream &stream) const {
    auto allowed_option_it = allowed_options.find(name);
    if (allowed_option_it!=allowed_options.end()) {

      // First print out the datatype
      stream << "> \"" << name << "\"          ["
             << GenericType::get_type_description(allowed_option_it->second)
             << "] ";

      // Print out the description on a new line.
      map<std::string, std::string>::const_iterator description_it =description_.find(name);
      if (description_it!=description_.end()) {
        if (description_it->second != "n/a")
          stream << "     \"" << description_it->second << "\""<< std::endl;
      }
    } else {
      stream << "  \"" << name << "\" does not exist.";
    }
  }

  void OptionsFunctionality::printOptions(ostream &stream) const {
    stream << "\"Option name\" [type] = value" << endl;
    for (auto it=allowed_options.begin(); it!=allowed_options.end(); ++it) {
      printOption(it->first, stream);
    }
    stream << endl;
  }

  std::vector<std::string> OptionsFunctionality::optionNames() const {
    std::vector<std::string> names;
    for (auto it=allowed_options.begin(); it!=allowed_options.end(); ++it) {
      names.push_back(it->first);
    }
    return names;
  }

} // namespace casadi
