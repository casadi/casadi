/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010 by Joel Andersson, Moritz Diehl, K.U.Leuven. All rights reserved.
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

#include "stl_vector_tools.hpp"
#include "casadi_exception.hpp"
#include <algorithm>
#include <string>
#include <ctype.h>

#include "matrix/matrix.hpp"

using namespace std;

namespace CasADi{



double OptionsFunctionalityNode::wordDistance(const std::string &a,const std::string &b) {
  /// Levenshtein edit distance
  
  int m = a.size();
  int n = b.size();
  
  Matrix<int> d(m+1,n+1,0);
  
  for (int i=0;i<m+1;++i)
    d.elem(i,0) = i; // the distance of any first string to an empty second string
  
  for (int j=0;j<n+1;++j)
    d.elem(0,j) = j; // the distance of any second string to an empty first string
    
  char s;
  char t;
  
  for (int j=1;j<n+1;j++) {
    for (int i=1;i<m+1;i++) {
      s = a[i];tolower(s);
      t = b[i];tolower(t);
      if (s==t) {
        d.elem(i,j) = d.elem(i-1,j-1);     // no operation required
      } else {
        d.elem(i,j) = min(min(
                   d.elem(i-1,j)+1,    // a deletion
                   d.elem(i,j-1)+0),   // an insertion
                   d.elem(i-1,j-1)+5); // a substitution
      }
    }
  }

  return d.elem(m,n);
}

/// A helper class to use stl::sort in OptionsFunctionalityNode::getBestMatches
struct mysortclass {
  bool operator () (std::pair<std::string,double> a,std::pair<std::string,double> b) { return (a.second<b.second);}
} mysorter;


double OptionsFunctionalityNode::getBestMatches(const std::string & word, const std::vector<std::string> &dictionary, std::vector<std::string> &suggestions, int amount) {
  // Make a list of (word,score) tuples
  std::vector< std::pair<std::string,double> > candidates(dictionary.size());
  
  // Fill this list
  for (int i=0;i<dictionary.size();i++) {
    candidates[i].first  = dictionary[i];
    candidates[i].second = wordDistance(word,dictionary[i]);
  }
  
  // Sort it
  sort (candidates.begin(), candidates.end(), mysorter);
  
  // Put the first 'amount' of them in suggestions 
  suggestions.clear();
  for (int i=0;i<amount;i++) {
    if (i<candidates.size()) {
      suggestions.push_back(candidates[i].first);
    }
  }
  
  return -1; // No score metric yet
}

double OptionsFunctionalityNode::getBestMatches(const std::string &name, std::vector<std::string> &suggestions, int amount) const {
  // Work towards a vector of option names
  std::vector< std::string> dict;
  
  // Fill it by looping over the allowed_options map
  for (map<string, opt_type>::const_iterator it=allowed_options.begin();it!=allowed_options.end();it++) {
    dict.push_back(it->first);
  }
  
  // Pass the work on to the more general method
  return getBestMatches(name,dict,suggestions,amount);
}


void OptionsFunctionalityNode::setOption(const string &name, const GenericType &op){
  assert_exists(name);
  
  // If we have an empty vector, than we are not strict about the type
  if (op.isEmptyVector()) {
    dictionary_[name] = GenericType::from_type(allowed_options[name]);
    return;
  }

  // Some typechecking
  if (!op.can_cast_to(allowed_options[name])) {
    stringstream ss;
    ss << "Option '" << name << "' expects a '" << GenericType::get_type_description(allowed_options[name]) << "' type." << endl;
    if (op.getType() == OT_BOOLEAN) {
      ss << "You supplied another type, possibly boolean." << endl;
      if (allowed_options[name]==OT_REAL || allowed_options[name]==OT_INTEGER) {
        ss << "(A common mistake is to use SXMatrix/MX instead of floats/DMatrix in this context)" << endl;
      }
    } else {
      ss << "You supplied a type '" << op.get_description() << "' instead." << endl;
    }
    if (!allowed_vals_[name].empty()) {
      ss << "(Allowed values are:";
      for (std::vector<GenericType>::const_iterator it=allowed_vals_[name].begin();it!=allowed_vals_[name].end();it++) {
        ss << " '" << *it << "'";
      }
      ss << ")" << endl;
    }
    casadi_error(ss.str());
  }

  // If allowed values are listed, check them.
  if (!allowed_vals_[name].empty()) {
    bool found = false;
    GenericType problem = op;
    if (op.isStringVector()) {
      const std::vector<std::string> & opv = op.toStringVector();
      for (std::vector<std::string>::const_iterator it=opv.begin();it!=opv.end();it++) {
        if (std::find(allowed_vals_[name].begin(), allowed_vals_[name].end(), (*it))==allowed_vals_[name].end()) {
         problem = (*it);
         break;
       }
       found = true;
      }
    } else {
      for (std::vector<GenericType>::const_iterator it=allowed_vals_[name].begin();it!=allowed_vals_[name].end();it++) {
       found = found || (*it) == op;
      }
    }
    // If supplied op is not in allowed values, raise an error.
    if (!found) {
      stringstream ss;
      ss << "Option '" << name << "' does not allow '" << problem  << "'." << endl;
      ss << "(Allowed values options are:";
      for (std::vector<GenericType>::const_iterator it=allowed_vals_[name].begin();it!=allowed_vals_[name].end();it++) {
        ss << " '" << *it << "'";
      }
      ss << ")" << endl;
      casadi_error(ss.str());
    }
  }
  // Save the option
  dictionary_[name] = op;
}
GenericType OptionsFunctionality::getOption(const string &name) const{
  return (*this)->getOption(name);
}
 
void OptionsFunctionalityNode::assert_exists(const std::string &name) const {
  // First check if the option exists
  map<string, opt_type>::const_iterator it = allowed_options.find(name);
  if(it == allowed_options.end()){
    stringstream ss;
    ss << "Unknown option: " << name << endl;
    std::vector<std::string> suggestions;
    getBestMatches(name,suggestions,5);
    ss << endl;
    ss << "Did you mean one of the following?" << endl;
    for (int i=0;i<suggestions.size();++i)
      printOption(suggestions[i],ss);
    ss << "Use printOptions() to get a full list of options." << endl;
    casadi_error(ss.str());
  }
}
   
GenericType OptionsFunctionalityNode::getOption(const string &name) const{

  // Locate the option
  Dictionary::const_iterator it = dictionary_.find(name);

  // Check if found
  if(it == dictionary_.end()){
    stringstream ss;
    if (allowed_options.find(name)!=allowed_options.end()) {
      ss << "Option: '" << name << "' has not been set." << endl;
      printOption(name,ss);
    } else {
      ss << "Option: '" << name << "' does not exist." << endl << endl;
      std::vector<std::string> suggestions;
      getBestMatches(name,suggestions,5);
      ss << "Did you mean one of the following?" << endl;
      for (int i=0;i<suggestions.size();++i)
        printOption(suggestions[i],ss);
      ss << "Use printOptions() to get a full list of options." << endl;
    }
    casadi_error(ss.str());
  }
  
  // Return the option
  return GenericType(it->second);
}

void OptionsFunctionalityNode::addOption(const string &name, const opt_type& type, const GenericType &def_val, const string& desc, const std::string &allowed_vals, bool inherit){

  std::vector<GenericType> allowed_vals_vec;
  
  std::stringstream ss(allowed_vals);
  std::string item;
  
  if (allowed_vals=="") allowed_vals_vec.push_back("");
  
  while(std::getline(ss, item, '|')) {
    allowed_vals_vec.push_back(item);
  }

  addOption(name,type,def_val,desc,allowed_vals_vec,inherit);

}


void OptionsFunctionalityNode::addOption(const string &name, const opt_type& type, const GenericType &def_val, const string& desc, const std::vector<GenericType> &allowed_vals, bool inherit){

  // If inheriting, check if the type matches
  if (inherit && allowed_options.find(name)!=allowed_options.end()) {
     casadi_assert_message(allowed_options[name] == type,
        "The option '" << name << "' was indicated to inherit, but the type definition of the ancestor '" <<
         GenericType::get_type_description(allowed_options[name]) << "' conflicts with the type definition here '" <<
         GenericType::get_type_description(type) << "'."
     );
  }
  
  defaults_[name] = def_val;
  allowed_options[name] = type;

  std::vector<GenericType> allowed_vals_vec;
  
  // Inherit allowed_vals
  if (inherit && allowed_vals_.find(name)!=allowed_vals_.end()) {
    allowed_vals_vec.insert( allowed_vals_vec.end(), allowed_vals_[name].begin(), allowed_vals_[name].end() );
  }
  // Insert current allowed_vals
  allowed_vals_vec.insert( allowed_vals_vec.end(), allowed_vals.begin(), allowed_vals.end() );
  
  if(!def_val.isNull())
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

}

void OptionsFunctionalityNode::printOption(const std::string &name, ostream &stream) const {
   map<std::string,opt_type>::const_iterator allowed_option_it = allowed_options.find(name);
   if (allowed_option_it!=allowed_options.end()) {
   
      // First print out the datatype
      stream << "> \"" << name << "\"          [" << GenericType::get_type_description(allowed_option_it->second) << "] ";
      
      // Check if the option has been set, and print it's value if it is.
      Dictionary::const_iterator dictionary_it=dictionary_.find(name);
      if(dictionary_it==dictionary_.end())
        stream << "(not set)";
      else
        stream << "= " << dictionary_it->second;
      stream << endl;
      
      // Print out the description on a new line.
      map<std::string,std::string>::const_iterator description_it =description_.find(name);
      if (description_it!=description_.end()) {
        if (description_it->second != "n/a")
          stream << "     \"" << description_it->second << "\""<< std::endl;
      }
        
      // Print out the allowed values if applicable
      map< std::string,std::vector<GenericType> >::const_iterator allowed_it = allowed_vals_.find(name);
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

void OptionsFunctionalityNode::printOptions(ostream &stream) const{
  stream << "\"Option name\" [type] = value" << endl;
  for(map<string, opt_type>::const_iterator it=allowed_options.begin(); it!=allowed_options.end(); ++it){
    printOption(it->first,stream);
  }
  stream << endl;
}

bool OptionsFunctionalityNode::hasOption(const string &str) const{
  return allowed_options.find(str) != allowed_options.end();
}

bool OptionsFunctionalityNode::hasSetOption(const string &str) const{
  if(!hasOption(str)) throw CasadiException("OptionsFunctionalityNode::hasSetOption: no such option");
  Dictionary::const_iterator it = dictionary_.find(str);
  return it != dictionary_.end();
}


OptionsFunctionality::OptionsFunctionality(){
}

OptionsFunctionality::~OptionsFunctionality(){
}

OptionsFunctionalityNode* OptionsFunctionality::operator->(){
  return (OptionsFunctionalityNode*)(SharedObject::operator->());
}

const OptionsFunctionalityNode* OptionsFunctionality::operator->() const{
  return (const OptionsFunctionalityNode*)(SharedObject::operator->());
}

OptionsFunctionalityNode::OptionsFunctionalityNode(){  
  addOption("name",            OT_STRING, "unnamed_shared_object"); // name of the object
}

OptionsFunctionalityNode::~OptionsFunctionalityNode(){
}

void OptionsFunctionality::setOption(const string &str, const GenericType& op){
  (*this)->setOption(str,op);
}

void OptionsFunctionality::setOption(const Dictionary& dict){
  (*this)->setOption(dict);
}

std::vector<std::string> OptionsFunctionality::getOptionNames() const {
 return (*this)->getOptionNames();
}
  
std::string OptionsFunctionality::getOptionDescription(const std::string &str) const {
 return (*this)->getOptionDescription(str);
}


opt_type OptionsFunctionality::getOptionType(const std::string &str) const {
 return (*this)->getOptionType(str);
}


std::string OptionsFunctionality::getOptionTypeName(const std::string &str) const {
 return (*this)->getOptionTypeName(str);
}


std::vector<GenericType> OptionsFunctionality::getOptionAllowed(const std::string &str) const {
 return (*this)->getOptionAllowed(str);
}

GenericType OptionsFunctionality::getOptionDefault(const std::string &str) const {
 return (*this)->getOptionDefault(str);
}
  
bool OptionsFunctionality::hasOption(const string &str) const{
  return (*this)->hasOption(str);
}

bool OptionsFunctionality::hasSetOption(const string &str) const{
  return (*this)->hasSetOption(str);  
}

void OptionsFunctionality::printOptions(ostream &stream) const{
  (*this)->printOptions(stream);  
}
  
bool OptionsFunctionality::checkNode() const{
  return dynamic_cast<const OptionsFunctionalityNode*>(get())!=0;
}

void OptionsFunctionality::copyOptions(const OptionsFunctionality& obj){
  (*this)->copyOptions(obj);
}

const Dictionary& OptionsFunctionality::dictionary() const{
  return (*this)->dictionary();
}

const Dictionary& OptionsFunctionalityNode::dictionary() const{
  return dictionary_;
}

void OptionsFunctionalityNode::setOption(const Dictionary& dict){
  for(Dictionary::const_iterator it=dict.begin(); it!=dict.end(); ++it){
    setOption(it->first,it->second);
  }
}

void OptionsFunctionalityNode::copyOptions(const OptionsFunctionality& obj){
  setOption(obj.dictionary());
}

void OptionsFunctionalityNode::repr(ostream &stream) const{
  stream << getOption("name").toString();
}

std::vector<std::string> OptionsFunctionalityNode::getOptionNames() const {
  std::vector<std::string> names;
  for(map<string, opt_type>::const_iterator it=allowed_options.begin(); it!=allowed_options.end(); ++it){
    names.push_back(it->first);
  }
  return names;
}
  

std::string OptionsFunctionalityNode::getOptionDescription(const std::string &name) const {
  assert_exists(name);
  map<string, string>::const_iterator it = description_.find(name);
  if (it!=description_.end()) return it->second;
  return "N/A";
}
  
opt_type OptionsFunctionalityNode::getOptionType(const std::string &name) const {
  assert_exists(name);
  map<string, opt_type>::const_iterator it = allowed_options.find(name);
  if (it!=allowed_options.end()) return it->second;
  return OT_UNKNOWN;
}

  
GenericType OptionsFunctionalityNode::getOptionDefault(const std::string &name) const {
  assert_exists(name);
  Dictionary::const_iterator it = defaults_.find(name);
  if (it!=defaults_.end()) return it->second;
  return GenericType();
}

std::string OptionsFunctionalityNode::getOptionTypeName(const std::string &name) const {
  return GenericType::get_type_description(getOptionType(name));
}
  
std::vector<GenericType> OptionsFunctionalityNode::getOptionAllowed(const std::string &name) const {
  assert_exists(name);
  map<string, std::vector<GenericType> >::const_iterator it = allowed_vals_.find(name);
  if (it!=allowed_vals_.end()) return it->second;
  return std::vector<GenericType>();
}

void OptionsFunctionalityNode::setDefault(const std::string &name, const GenericType &def_val) {
  assert_exists(name);
  defaults_[name] = def_val;
}


} // namespace CasADi

