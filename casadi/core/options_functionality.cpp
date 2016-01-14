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

  double OptionsFunctionalityNode::wordDistance(const std::string &a, const std::string &b) {
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
  /// A helper class to use stl::sort in OptionsFunctionalityNode::getBestMatches
  struct mysortclass {
    bool operator()(std::pair<std::string, double> a, std::pair<std::string, double> b) {
      return (a.second<b.second);}
  } mysorter;
  /// \endcond

  double OptionsFunctionalityNode::
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

  double OptionsFunctionalityNode::getBestMatches(const std::string &name,
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

  void setAdaptorOptions(Dict& dict, const string &name, const Dict &op) {

    // Find position of '.' separator
    std::string::size_type dotpos = name.find(".");

    // Get the adaptor name (before '.', or the entire name if '.' not found)
    std::string adaptor_name = name.substr(0, dotpos);

    // Check if adaptor name already occurs in the dictionary of options
    Dict::const_iterator it = dict.find(adaptor_name);

    if (it == dict.end()) {
      // Create an empty dictionary if not
      dict[adaptor_name] = Dict();
    } else if (!dict[adaptor_name].is_dict()) {
      // If an entry is found, make sure it is a dictionary
      casadi_error("setAdaptorOptions: Dict expected, but got " << dict[adaptor_name] << ".");
    }

    if (dotpos==std::string::npos) {
      // We reached the end of the dotted adaptor string
      Dict target = dict[adaptor_name];
      // Merge the contents of the supplied dictionary
      for (Dict::const_iterator it=op.begin(); it!=op.end(); ++it) {
        target[it->first] = it->second;
      }
      dict[adaptor_name] = target;
    } else {
      // Descend one level down
      Dict dict_adaptor = dict[adaptor_name];
      setAdaptorOptions(dict_adaptor, name.substr(dotpos+1), op);
      dict[adaptor_name] = dict_adaptor;
    }
  }

  void OptionsFunctionalityNode::assert_exists(const std::string &name) const {
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

  void OptionsFunctionalityNode::
  addOption(const string &name, const TypeID& type, const string& desc) {
    allowed_options[name] = type;
    description_[name] = desc;
  }

  void OptionsFunctionalityNode::printOption(const std::string &name, ostream &stream) const {
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

  void OptionsFunctionalityNode::printOptions(ostream &stream) const {
    stream << "\"Option name\" [type] = value" << endl;
    for (auto it=allowed_options.begin(); it!=allowed_options.end(); ++it) {
      printOption(it->first, stream);
    }
    stream << endl;
  }

  bool OptionsFunctionalityNode::hasOption(const string &str) const {
    return allowed_options.find(str) != allowed_options.end();
  }

  OptionsFunctionalityNode::OptionsFunctionalityNode() {
  }

  OptionsFunctionalityNode::~OptionsFunctionalityNode() {
  }

  std::vector<std::string> OptionsFunctionalityNode::optionNames() const {
    std::vector<std::string> names;
    for (auto it=allowed_options.begin(); it!=allowed_options.end(); ++it) {
      names.push_back(it->first);
    }
    return names;
  }

  std::string OptionsFunctionalityNode::optionDescription(const std::string &name) const {
    assert_exists(name);
    map<string, string>::const_iterator it = description_.find(name);
    if (it!=description_.end()) return it->second;
    return "N/A";
  }

  TypeID OptionsFunctionalityNode::optionType(const std::string &name) const {
    assert_exists(name);
    auto it = allowed_options.find(name);
    if (it!=allowed_options.end()) return it->second;
    return OT_UNKNOWN;
  }

  std::string OptionsFunctionalityNode::optionTypeName(const std::string &name) const {
    return GenericType::get_type_description(optionType(name));
  }

} // namespace casadi
