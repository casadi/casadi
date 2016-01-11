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

double OptionsFunctionalityNode::getBestMatches(const std::string & word,
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
  for (map<string, TypeID>::const_iterator it=allowed_options.begin();
       it!=allowed_options.end();it++) {
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
  } else if (!dict[adaptor_name].isDict()) {
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


void OptionsFunctionalityNode::setOption(const string &name, const GenericType &op) {

  assert_exists(name);

  // For options set using dictionary shorthand, no type checking
  auto dotpos = name.find('.');
  if (dotpos != string::npos) {
    dictionary_[name] = op;
    return;
  }

  // If we have an empty vector, than we are not strict about the type
  if (op.is_emptyVector()) {
    dictionary_[name] = GenericType::from_type(allowed_options[name]);
    return;
  }

  // Some typechecking
  if (!op.can_cast_to(allowed_options[name]) && !op.isNull()) {
    stringstream ss;
    ss << "Option '" << name << "' expects a '" <<
        GenericType::get_type_description(allowed_options[name]) << "' type." << endl;
    if (op.getType() == OT_BOOLEAN) {
      ss << "You supplied another type, possibly boolean." << endl;
      if (allowed_options[name]==OT_REAL || allowed_options[name]==OT_INTEGER) {
        ss << "(A common mistake is to use SX/MX instead of floats/DM in this context)"
           << endl;
      }
    } else {
      ss << "You supplied a type '" << op.get_description() << "' instead." << endl;
    }
    if (!allowed_vals_[name].empty()) {
      ss << "(Allowed values are:";
      for (std::vector<GenericType>::const_iterator it=allowed_vals_[name].begin();
           it!=allowed_vals_[name].end();it++) {
        ss << " '" << *it << "'";
      }
      ss << ")" << endl;
    }
    casadi_error(ss.str());
  }

  // If allowed values are listed, check them.
  if (!allowed_vals_[name].empty()) {
    bool found;
    GenericType problem = op;
    if (op.isStringVector()) {
      found = true;
      const std::vector<std::string> & opv = op.toStringVector();
      for (std::vector<std::string>::const_iterator it=opv.begin();it!=opv.end();it++) {
        userOut() << "checking " << *it << std::endl;
        if (std::find(allowed_vals_[name].begin(),
                      allowed_vals_[name].end(), (*it))==allowed_vals_[name].end()) {
          problem = (*it);
          found = false;
          break;
        }
      }
    } else {
      found = false;
      for (std::vector<GenericType>::const_iterator it=allowed_vals_[name].begin();
           it!=allowed_vals_[name].end();it++) {
       found = found || (*it) == op;
      }
    }
    // If supplied op is not in allowed values, raise an error.
    if (!found) {
      stringstream ss;
      ss << "Option '" << name << "' does not allow '" << problem  << "'." << endl;
      ss << "(Allowed values options are:";
      for (std::vector<GenericType>::const_iterator it=allowed_vals_[name].begin();
           it!=allowed_vals_[name].end();it++) {
        ss << " '" << *it << "'";
      }
      ss << ")" << endl;
      casadi_error(ss.str());
    }
  }
  // Save the option
  dictionary_[name] = op;
}
GenericType OptionsFunctionality::option(const string &name) const {
  return (*this)->option(name);
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

GenericType OptionsFunctionalityNode::option(const string &name) const {

  // Locate the option
  Dict::const_iterator it = dictionary_.find(name);

  // Return if found
  if (it != dictionary_.end()) {
    return GenericType(it->second);
  }

  // Check if a dictionary
  string dotname = name + ".";
  it = dictionary_.upper_bound(dotname);
  if (it!=dictionary_.end() && it->first.compare(0, dotname.size(), dotname)==0) {
    // Dictionary option
    Dict ret;
    while (it!=dictionary_.end() && it->first.compare(0, dotname.size(), dotname)==0) {
      ret[it->first.substr(dotname.size())] = it->second;
      it++;
    }
    return ret;
  }

  // Error
  stringstream ss;
  if (allowed_options.find(name)!=allowed_options.end()) {
    ss << "Option: '" << name << "' has not been set." << endl;
    printOption(name, ss);
  } else {
    ss << "Option: '" << name << "' does not exist." << endl << endl;
    std::vector<std::string> suggestions;
    getBestMatches(name, suggestions, 5);
    ss << "Did you mean one of the following?" << endl;
    for (int i=0;i<suggestions.size();++i)
      printOption(suggestions[i], ss);
    ss << "Use printOptions() to get a full list of options." << endl;
  }
  casadi_error(ss.str());
  return GenericType();
}

void OptionsFunctionalityNode::addOption(const string &name, const TypeID& type,
                                         const GenericType &def_val, const string& desc,
                                         const std::string &allowed_vals, bool inherit) {

  std::vector<GenericType> allowed_vals_vec;
  std::vector<int> enum_values;
  std::vector<std::string> enum_descr;

  std::stringstream ss(allowed_vals);
  std::string item;

  if (allowed_vals=="") allowed_vals_vec.push_back("");

  int val = -1;
  while (std::getline(ss, item, '|')) {
    std::stringstream sss(item);
    std::string name;
    std::getline(sss, name, ':');
    allowed_vals_vec.push_back(name);
    std::string descr="";
    std::getline(sss, descr, ':');
    enum_descr.push_back(descr);
    std::string value;
    if (std::getline(sss, value, ':')) {
      std::istringstream(value) >> val;
    } else {
      val++;
    }
    enum_values.push_back(val);
  }

  addOption(name, type, def_val, desc, allowed_vals_vec, inherit, enum_values, enum_descr);

}


void OptionsFunctionalityNode::addOption(
    const string &name, const TypeID& type, const GenericType &def_val,
    const string& desc, const std::vector<GenericType> &allowed_vals,
    bool inherit, std::vector<int> enum_values,
    std::vector<std::string> enum_descr) {
  // If inheriting, check if the type matches
  if (inherit && allowed_options.find(name)!=allowed_options.end()) {
     casadi_assert_message(allowed_options[name] == type,
        "The option '" << name
         << "' was indicated to inherit, but the type definition of the ancestor '"
         << GenericType::get_type_description(allowed_options[name])
         << "' conflicts with the type definition here '"
         << GenericType::get_type_description(type) << "'.");
  }

  defaults_[name] = def_val;
  allowed_options[name] = type;

  std::vector<GenericType> allowed_vals_vec;
  std::vector<int> enum_values_vec;
  std::vector<std::string> enum_descr_vec;

  // Inherit
  if (inherit && allowed_vals_.find(name)!=allowed_vals_.end()) {
    allowed_vals_vec.insert(allowed_vals_vec.end(), allowed_vals_[name].begin(),
                            allowed_vals_[name].end());
  }
  if (inherit && enum_descr_.find(name)!=enum_descr_.end()) {
    enum_descr_vec.insert(enum_descr_vec.end(), enum_descr_[name].begin(),
                          enum_descr_[name].end());
  }
  if (inherit && enum_values_.find(name)!=enum_values_.end()) {
    enum_values_vec.insert(enum_values_vec.end(), enum_values_[name].begin(),
                           enum_values_[name].end());
  }
  // Insert current allowed_vals
  allowed_vals_vec.insert(allowed_vals_vec.end(), allowed_vals.begin(), allowed_vals.end());
  enum_descr_vec.insert(enum_descr_vec.end(), enum_descr.begin(), enum_descr.end());
  enum_values_vec.insert(enum_values_vec.end(), enum_values.begin(), enum_values.end());

  if (!def_val.isNull())
    dictionary_[name] = def_val;

  // Inherit description
  std::stringstream s;
  if (inherit && description_.find(name)!=description_.end()) {
    s << description_[name];
    if (!desc.empty())
      s << std::endl;
  }
  // Insert current description
  s << desc;
  description_[name] = s.str();

  allowed_vals_[name] = allowed_vals_vec;

  enum_values_[name] = enum_values_vec;

  enum_descr_[name] = enum_descr_vec;

  casadi_assert(enum_values_vec.empty() || enum_values_vec.size() == allowed_vals_vec.size());
  casadi_assert(enum_descr_vec.empty() || enum_descr_vec.size() == allowed_vals_vec.size());

}

void OptionsFunctionalityNode::printOption(const std::string &name, ostream &stream) const {
   map<std::string, TypeID>::const_iterator allowed_option_it = allowed_options.find(name);
   if (allowed_option_it!=allowed_options.end()) {

      // First print out the datatype
      stream << "> \"" << name << "\"          ["
             << GenericType::get_type_description(allowed_option_it->second)
             << "] ";

      // Check if the option has been set, and print it's value if it is.
      Dict::const_iterator dictionary_it=dictionary_.find(name);
      if (dictionary_it==dictionary_.end())
        stream << "(not set)";
      else
        stream << "= " << dictionary_it->second;
      stream << endl;

      // Print out the description on a new line.
      map<std::string, std::string>::const_iterator description_it =description_.find(name);
      if (description_it!=description_.end()) {
        if (description_it->second != "n/a")
          stream << "     \"" << description_it->second << "\""<< std::endl;
      }

      // Print out the allowed values if applicable
      map< std::string, std::vector<GenericType> >::const_iterator allowed_it =
          allowed_vals_.find(name);
      if (allowed_it!=allowed_vals_.end()) {
        const std::vector<GenericType> & allowed = allowed_it->second;
        if (allowed.size()>0) {
          stream << "     Allowed values: ";
          for (std::vector<GenericType>::const_iterator it=allowed.begin();it!=allowed.end();it++) {
             stream << " '" << *it << "'";
          }
          stream << std::endl;
        }
      }
   } else {
     stream << "  \"" << name << "\" does not exist.";
   }
}

void OptionsFunctionalityNode::printOptions(ostream &stream) const {
  stream << "\"Option name\" [type] = value" << endl;
  for (map<string, TypeID>::const_iterator it=allowed_options.begin();
      it!=allowed_options.end(); ++it) {
    printOption(it->first, stream);
  }
  stream << endl;
}

bool OptionsFunctionalityNode::hasOption(const string &str) const {
  return allowed_options.find(str) != allowed_options.end();
}

bool OptionsFunctionalityNode::hasSetOption(const string &str) const {
  if (!hasOption(str)) casadi_error("OptionsFunctionalityNode::hasSetOption: no such option '"
                                   << str << "'");
  Dict::const_iterator it = dictionary_.find(str);
  if (it!=dictionary_.end()) return true;

  // Check if a dictionary
  string dotstr = str + ".";
  it = dictionary_.upper_bound(dotstr);
  return it!=dictionary_.end() && it->first.compare(0, dotstr.size(), dotstr)==0;
}

OptionsFunctionality::OptionsFunctionality() {
}

OptionsFunctionality::~OptionsFunctionality() {
}

OptionsFunctionalityNode* OptionsFunctionality::operator->() {
  return static_cast<OptionsFunctionalityNode*>(SharedObject::operator->());
}

const OptionsFunctionalityNode* OptionsFunctionality::operator->() const {
  return static_cast<const OptionsFunctionalityNode*>(SharedObject::operator->());
}

OptionsFunctionalityNode::OptionsFunctionalityNode() {
  addOption("defaults_recipes",    OT_STRINGVECTOR, GenericType(),
            "Changes default options according to a given recipe (low-level)");
}

OptionsFunctionalityNode::~OptionsFunctionalityNode() {
}

void OptionsFunctionality::setOption(const string &str, const GenericType& op) {
  (*this)->setOption(str, op);
}

void OptionsFunctionality::setOption(const Dict& dict, bool skipUnknown) {
  (*this)->setOption(dict, skipUnknown);
}

std::vector<std::string> OptionsFunctionality::optionNames() const {
 return (*this)->optionNames();
}

std::string OptionsFunctionality::optionDescription(const std::string &str) const {
 return (*this)->optionDescription(str);
}


TypeID OptionsFunctionality::optionType(const std::string &str) const {
 return (*this)->optionType(str);
}


std::string OptionsFunctionality::optionTypeName(const std::string &str) const {
 return (*this)->optionTypeName(str);
}


std::vector<GenericType> OptionsFunctionality::optionAllowed(const std::string &str) const {
 return (*this)->optionAllowed(str);
}

GenericType OptionsFunctionality::optionDefault(const std::string &str) const {
 return (*this)->optionDefault(str);
}

bool OptionsFunctionality::hasOption(const string &str) const {
  return (*this)->hasOption(str);
}

bool OptionsFunctionality::hasSetOption(const string &str) const {
  return (*this)->hasSetOption(str);
}

void OptionsFunctionality::printOptions(ostream &stream) const {
  (*this)->printOptions(stream);
}

bool OptionsFunctionality::test_cast(const SharedObjectNode* ptr) {
  return dynamic_cast<const OptionsFunctionalityNode*>(ptr)!=0;
}

void OptionsFunctionality::copyOptions(const OptionsFunctionality& obj, bool skipUnknown) {
  (*this)->copyOptions(obj, skipUnknown);
}

const Dict& OptionsFunctionality::dictionary() const {
  return (*this)->dictionary();
}

int OptionsFunctionality::optionAllowedIndex(const std::string &name) const {
  return (*this)->optionAllowedIndex(name);
}

void OptionsFunctionality::setOptionByAllowedIndex(const std::string &name, int i) {
  return (*this)->setOptionByAllowedIndex(name, i);
}

int OptionsFunctionality::optionEnumValue(const std::string &name) const {
  return (*this)->optionEnumValue(name);
}

void OptionsFunctionality::setOptionByEnumValue(const std::string &name, int v) {
  return (*this)->setOptionByEnumValue(name, v);
}

const Dict& OptionsFunctionalityNode::dictionary() const {
  return dictionary_;
}

void OptionsFunctionalityNode::setOption(const Dict& dict, bool skipUnknown) {
  for (Dict::const_iterator it=dict.begin(); it!=dict.end(); ++it) {
    if (!skipUnknown || hasOption(it->first)) {
      setOption(it->first, it->second);
    }
  }
}

void OptionsFunctionalityNode::copyOptions(const OptionsFunctionality& obj, bool skipUnknown) {
  setOption(obj.dictionary(), skipUnknown);
}

std::vector<std::string> OptionsFunctionalityNode::optionNames() const {
  std::vector<std::string> names;
  for (map<string, TypeID>::const_iterator it=allowed_options.begin();
      it!=allowed_options.end(); ++it) {
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
  map<string, TypeID>::const_iterator it = allowed_options.find(name);
  if (it!=allowed_options.end()) return it->second;
  return OT_UNKNOWN;
}


GenericType OptionsFunctionalityNode::optionDefault(const std::string &name) const {
  assert_exists(name);
  Dict::const_iterator it = defaults_.find(name);
  if (it!=defaults_.end()) return it->second;
  return GenericType();
}

std::string OptionsFunctionalityNode::optionTypeName(const std::string &name) const {
  return GenericType::get_type_description(optionType(name));
}

std::vector<GenericType> OptionsFunctionalityNode::optionAllowed(const std::string &name) const {
  assert_exists(name);
  map<string, std::vector<GenericType> >::const_iterator it = allowed_vals_.find(name);
  if (it!=allowed_vals_.end()) return it->second;
  return std::vector<GenericType>();
}

int OptionsFunctionalityNode::optionAllowedIndex(const std::string &name) const {
  assert_exists(name);
  casadi_assert_message(hasSetOption(name), "Option '" << name << "' has not been set.");
  map<string, std::vector<GenericType> >::const_iterator it = allowed_vals_.find(name);
  casadi_assert_message(it!=allowed_vals_.end(), "Option '" << name
                        << "' has no list of allowed values.");
  const std::vector<GenericType> &vec = it->second;
  std::vector<GenericType>::const_iterator it2 = std::find(vec.begin(), vec.end(), option(name));
  return it2-vec.begin();
}

void OptionsFunctionalityNode::setOptionByAllowedIndex(const std::string &name, int i) {
  assert_exists(name);
  map<string, std::vector<GenericType> >::const_iterator it = allowed_vals_.find(name);
  casadi_assert_message(it!=allowed_vals_.end(), "Option '" << name
                        << "' has no list of allowed values.");
  const std::vector<GenericType> &vec = it->second;
  casadi_assert_message(i>=0 && i<= vec.size()-1, "setOptionAllowedIndex('" << name
                        << "', " << i << "): index out of bounds. There are "
                        << vec.size() << " allowed values.");
  setOption(name, vec[i]);
}

int OptionsFunctionalityNode::optionEnumValue(const std::string &name) const {
  assert_exists(name);
  int i = optionAllowedIndex(name);
  map<string, std::vector<int> >::const_iterator it = enum_values_.find(name);
  casadi_assert_message(it!=enum_values_.end(), "Option '" << name
                        << "' has no list of enum values.");
  const std::vector<int> & enum_values = it->second;
  casadi_assert_message(!enum_values.empty(), "Option '" << name
                        << "' has an empty enum values list.");
  return enum_values.at(i);
}

void OptionsFunctionalityNode::setOptionByEnumValue(const std::string &name, int v) {
  assert_exists(name);
  map<string, std::vector<int> >::const_iterator it = enum_values_.find(name);
  casadi_assert_message(it!=enum_values_.end(), "Option '" << name
                        << "' has no list of enum values.");
  const std::vector<int> & enum_values = it->second;
  casadi_assert_message(!enum_values.empty(), "Option '" << name
                        << "' has an empty enum values list.");
  std::vector<int>::const_iterator it2 = std::find(enum_values.begin(), enum_values.end(), v);
  casadi_assert_message(it2!=enum_values.end(), "Option '"
                        << name << "', entry " << v << " was not found in enum value list.");
  setOptionByAllowedIndex(name, it2-enum_values.begin());
}


void OptionsFunctionalityNode::setDefault(const std::string &name, const GenericType &def_val) {
  assert_exists(name);
  defaults_[name] = def_val;
}

Dict OptionsFunctionality::addOptionRecipe(const Dict& dict, const std::string& recipe) {
  Dict ret = dict;
  Dict::const_iterator f = dict.find("defaults_recipes");
  if (f!=dict.end()) {
    std::vector<std::string> defaults_recipes = f->second;
    defaults_recipes.push_back(recipe);
    ret["defaults_recipes"] = defaults_recipes;
  } else {
    ret["defaults_recipes"] = std::vector<std::string>(1, recipe);
  }
  return ret;
}

void OptionsFunctionalityNode::setDefaultOptions() {
    if (hasSetOption("defaults_recipes")) {
      const std::vector<std::string> & recipes = option("defaults_recipes");
      setDefaultOptions(recipes);
    }
}

} // namespace casadi
