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

const char *opt_type_name[] =  { "boolean", "integer", "real", "string", "integervector", "realvector", "Dictionary", "NLPSolver", "LinearSolver", "Integrator", "QPSolver", "ImplicitSolver", "JacobianGenerator", "SparsityGenerator", "void pointer"};

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
  
  // Some typechecking
  //if (allowed_options[name] == OT_STRING && op)

  // If allowed values are listed, check them.
  if (!allowed_vals_[name].empty()) {
    bool found = false;
    for (std::vector<GenericType>::const_iterator it=allowed_vals_[name].begin();it!=allowed_vals_[name].end();it++) {
     found = found || (*it) == op;
    }
    // If supplied op is not in allowed values, raise an error.
    if (!found) {
      stringstream ss;
      ss << "Option '" << name << "' does not allow '" << op  << "'." << endl;
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

void OptionsFunctionalityNode::addOption(const string &name, const opt_type& type, const GenericType &def_val, const string& desc, const std::string &allowed_vals){
  std::vector<GenericType> allowed_vals_vec;
  
  std::stringstream ss(allowed_vals);
  std::string item;
  
  while(std::getline(ss, item, '|')) {
      allowed_vals_vec.push_back(item);
  }

  addOption(name,type,def_val,desc,allowed_vals_vec);

}


void OptionsFunctionalityNode::addOption(const string &name, const opt_type& type, const GenericType &def_val, const string& desc, const std::vector<GenericType> &allowed_vals){
  allowed_options[name] = type;

  if(!def_val.isNull())
    dictionary_[name] = def_val;

  description_[name] = desc;
  
  allowed_vals_[name] = allowed_vals;
}

void OptionsFunctionalityNode::printOptions(ostream &stream) const{
  // Print allowed options
  stream << "\"Option name\" [type] = value" << endl;
  for(map<string, opt_type>::const_iterator it=allowed_options.begin(); it!=allowed_options.end(); ++it){
    stream << "  \"" << it->first << "\" [" << opt_type_name[it->second] << "] ";
    
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

