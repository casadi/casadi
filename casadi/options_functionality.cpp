#include "options_functionality.hpp"
#include <cassert>
#include "stl_vector_tools.hpp"
#include "casadi_exception.hpp"

const char *opt_type_name[] =  { "boolean", "integer", "real", "string", "integervector", "realvector" };

using namespace std;

namespace CasADi{

void OptionsFunctionalityNode::setOption(const string &name, const Option &op){
  // First check if the option exists
  map<string, opt_type>::const_iterator it = allowed_options.find(name);
  if(it == allowed_options.end()){
    stringstream ss;
    ss << "Unknown option: " << name << endl;
    throw CasadiException(ss.str());
  }

  // Save the option
  options[name] = op;
}
Option OptionsFunctionality::getOption(const string &name) const{
  return (*this)->getOption(name);
}
 
Option OptionsFunctionalityNode::getOption(const string &name) const{

  // Locate the option
  map<string, Option>::const_iterator it = options.find(name);

  // Check if found
  if(it == options.end()){
    stringstream ss;
    ss << "Option: " << name << " has not been set." << endl;
    throw CasadiException(ss.str());
  }
  
  // Return the option
  return Option(it->second);
}

void OptionsFunctionalityNode::addOption(const string &name, const opt_type& type, const Option &def_val){
  allowed_options[name] = type;

  if(!def_val.isNull())
    options[name] = Option(def_val);

}

void OptionsFunctionalityNode::printOptions(ostream &stream) const{
  // Print allowed options
  stream << "\"Option name\" [type] = value" << endl;
  for(map<string, opt_type>::const_iterator it=allowed_options.begin(); it!=allowed_options.end(); ++it){
    stream << "  \"" << it->first << "\" [" << opt_type_name[it->second] << "] ";
    
    // Check if it is has been set
    map<string, Option>::const_iterator j=options.find(it->first);
    if(j==options.end())
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
  map<string, Option>::const_iterator it = options.find(str);
  return it != options.end();
}


OptionsFunctionality::OptionsFunctionality(){
}

OptionsFunctionality::~OptionsFunctionality(){
}

OptionsFunctionalityNode* OptionsFunctionality::get(){
  return (OptionsFunctionalityNode*)(SharedObject::get());
}

const OptionsFunctionalityNode* OptionsFunctionality::get() const{
  return (const OptionsFunctionalityNode*)(SharedObject::get());
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

void OptionsFunctionality::setOption(const string &str, const Option& op){
  (*this)->setOption(str,op);
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
  
void OptionsFunctionality::assertNode() const{
  if(!dynamic_cast<const OptionsFunctionalityNode*>(get()))
    throw CasadiException("OptionsFunctionality::assertNode");
}

} // namespace CasADi

