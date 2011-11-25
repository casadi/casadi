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

using namespace std;

namespace CasADi{

void OptionsFunctionalityNode::setOption(const string &name, const GenericType &op){
  // First check if the option exists
  map<string, opt_type>::const_iterator it = allowed_options.find(name);
  if(it == allowed_options.end()){
    stringstream ss;
    ss << "Unknown option: " << name << endl;
    ss << "(Available options are:";
    for (map<string, opt_type>::const_iterator it=allowed_options.begin();it!=allowed_options.end();it++) {
      ss << " " << it->first;
    }
    ss << ")" << endl;
    casadi_error(ss.str());
  }
  
  // If we have an empty vector, than we are not strict about the type
  if (op.isEmptyVector()) {
    dictionary_[name] = GenericType::from_type(allowed_options[name]);
    return;
  }

  // Some typechecking
  if (!op.can_cast_to(allowed_options[name])) {
    stringstream ss;
    ss << "Option '" << name << "' expects a '" << GenericType::get_type_description(allowed_options[name]) << "' type." << endl;
    ss << "You supplied a type '" << op.get_description() << "' instead." << endl;
    if (!allowed_vals_[name].empty()) {
      ss << "(Allowed values options are:";
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
 
GenericType OptionsFunctionalityNode::getOption(const string &name) const{

  // Locate the option
  Dictionary::const_iterator it = dictionary_.find(name);

  // Check if found
  if(it == dictionary_.end()){
    stringstream ss;
    ss << "Option: " << name << " has not been set." << endl;
    ss << "(Available options are:";
    for (map<string, opt_type>::const_iterator it=allowed_options.begin();it!=allowed_options.end();it++) {
      ss << " " << it->first;
    }
    ss << ")" << endl;
    casadi_error(ss.str());
  }
  
  // Return the option
  return GenericType(it->second);
}

void OptionsFunctionalityNode::addOption(const string &name, const opt_type& type, const GenericType &def_val, const string& desc, const std::string &allowed_vals, bool inherit){

  std::vector<GenericType> allowed_vals_vec;
  
  std::stringstream ss(allowed_vals);
  std::string item;
  
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
    s << description_[name] << std::endl;
  }
  // Insert current description
  s << desc;
  description_[name] = s.str();

  allowed_vals_[name] = allowed_vals_vec;

}

void OptionsFunctionalityNode::printOptions(ostream &stream) const{
  // Print allowed options
  stream << "\"Option name\" [type] = value" << endl;
  for(map<string, opt_type>::const_iterator it=allowed_options.begin(); it!=allowed_options.end(); ++it){
    stream << "  \"" << it->first << "\" [" << GenericType::get_type_description(it->second) << "] ";
    
    // Check if it is has been set
    Dictionary::const_iterator j=dictionary_.find(it->first);
    if(j==dictionary_.end())
      stream << "(not set)";
    else
      stream << "= " << j->second;
    
    stream << endl;
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


} // namespace CasADi

