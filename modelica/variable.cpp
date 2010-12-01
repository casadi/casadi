#include "variable.hpp"
#include <cassert>
#include "../casadi/expression_tools.hpp"
#include "../casadi/casadi_exception.hpp"

using namespace std;
namespace CasADi{
namespace Modelica{
  
const char* VariableNode::typenames[] = {"time","differential state","algebraic state","control","parameter","dependent","derivative","not set"};
Variable::Variable(){
}

Variable::Variable(const string& name){
  assignNode(new VariableNode(name));
}

void Variable::setType(VarType type){
  (*this)->type = type;
}

Variable::~Variable(){
}

VariableNode* Variable::operator->(){
  return (VariableNode*)(SharedObject::operator->());
}

const VariableNode* Variable::operator->() const{
  return (const VariableNode*)(SharedObject::operator->());
}
  
VariableNode::VariableNode(const string& name) : name_(name){
  sx_ = SX(name);

  variability_ = CONTINUOUS;
  causality_ = INTERNAL;
  alias_ = NO_ALIAS;
  description_ = "";
  valueReference_ -1; //?
  min_ = -numeric_limits<double>::infinity();
  max_ = numeric_limits<double>::infinity();
  nominal_ = 1.0;
  start_ = 0.0;
  unit_ = "";
  displayUnit_ = "";

  type = TYPE_NOT_SET;
}

void VariableNode::init(){
  // Determine the type
  if(variability_ == PARAMETER){
    type = TYPE_PARAMETER;
  } else if(variability_ == CONTINUOUS) {
      if(causality_ == INTERNAL){
        type = TYPE_ALGEBRAIC;
      } else if(causality_ == INPUT){
        type = TYPE_CONTROL;
      }
      
      // Continuous variables have derivatives
      if(d.isNull()){
        stringstream ss;
        ss << "der(" << name_ << ")";
        d = Variable(ss.str());
        d.setType(TYPE_DERIVATIVE);
      }
      
    } else if (variability_ == CONSTANT){
      type = TYPE_DEPENDENT;
    }
}
  
VariableNode::~VariableNode(){
}

Variable::operator SX() const{
  return sx();  
}

SX Variable::der() const{
  return (*this)->der();  
}

SX Variable::sx() const{
  return (*this)->sx_;  
}



const string& VariableNode::getName() const{
  return name_;
}

string VariableNode::getType() const{
  return typenames[type];
}


bool Variable::isTime() const{
  return (*this)->type == TYPE_TIME;
}

bool Variable::isDifferentialState() const{
  return (*this)->type == TYPE_STATE;
}

bool Variable::isAlgebraicState() const{
  return (*this)->type == TYPE_ALGEBRAIC;  
}

bool Variable::isControl() const{
  return (*this)->type == TYPE_CONTROL;     
}

bool Variable::isParameter() const{
  return (*this)->type == TYPE_PARAMETER;    
}

bool Variable::isDependent() const{
  return (*this)->type == TYPE_DEPENDENT;     
}

SX VariableNode::der() const{
  if(type!=TYPE_STATE){
    stringstream ss;
    ss << "VariableNode::der(): Cannot take derivative of " << sx_;
    throw CasadiException(ss.str());
  } else {
    return d->sx_;
  }
}

double Variable::getValue() const{
  return (*this)->val;
}

void Variable::setValue(double val){
  (*this)->val = val;
}
  
Variable Variable::operator()(const string& name) const{
  // try to locate the variable
  map<string, int>::const_iterator it = (*this)->name_part.find(name);

  // check if the variable exists
  if(it==(*this)->name_part.end()){
    stringstream ss;
    ss << "Variable::operator(): No such variable: " << name;
    throw CasadiException(ss.str());
  }

  // return the variable
  return (*this)->col.at(it->second);  
}

Variable Variable::operator[](int ind) const{
  try{
    int base = 1;
    return (*this)->col.at(ind-base);
  } catch(exception& ex){
    throw CasadiException(string("Variable::operator[] failed <= ") + ex.what());
  }
}

int VariableNode::add(const Variable& var, const string& namepart){
  col.push_back(var);
  int ind = col.size()-1;
  name_part[namepart] = ind;
  return ind;
}

bool VariableNode::has(const string& name) const{
  // try to locate the variable
  map<string, int>::const_iterator it = name_part.find(name);

  // check if the variable exists
  return it!=name_part.end();
}

void VariableNode::repr(ostream &stream) const{
  stream << name_;
}


void VariableNode::print(ostream &stream) const{
  print(stream,0);
}

void VariableNode::print(ostream &stream, int indent) const{
  // Add indentation
  for(int i=0; i<indent; ++i) stream << "  ";
  
  // Print name
  stream << name_ << ": " << getType() << endl;

  if(!col.empty()){
    for(vector<Variable>::const_iterator it = col.begin(); it!=col.end(); ++it){
      (*it)->print(stream,indent+2);
    }
  }
}

void VariableNode::getAll(vector<Variable>& vars) const{
  if(col.empty()){
    Variable temp;
    temp.assignNode(const_cast<VariableNode*>(this)); // ugly trick
    vars.push_back(temp);
  } else {
    for(vector<Variable>::const_iterator it=col.begin(); it!=col.end(); ++it)
      if(!it->isNull())
        (*it)->getAll(vars);
  }
}

Variable::operator vector<Variable>() const{
  if(isNull())
    return vector<Variable>();
  else{
    vector<Variable> ret;
    (*this)->getAll(ret);
    return ret;
  }
}

const string& Variable::getName() const{
  return (*this)->name_;
}

void Variable::setName(const string& name){
  (*this)->name_ = name;
}

Variability Variable::getVariability() const{
  return (*this)->variability_;
}

void Variable::setVariability(Variability variability){
  (*this)->variability_ = variability;
}

Causality Variable::getCausality() const{
  return (*this)->causality_;
}

void Variable::setCausality(Causality causality){
  (*this)->causality_ = causality;
}
    
Alias Variable::getAlias() const{
  return (*this)->alias_;
}

void Variable::setAlias(Alias alias){
  (*this)->alias_ = alias;
}
    
const string& Variable::getDescription() const{
  return (*this)->description_;
}

void Variable::setDescription(const string& description){
  (*this)->description_ = description;
}
    
int Variable::getValueReference() const{
  return (*this)->valueReference_;
}

void Variable::setValueReference(int valueReference){
  (*this)->valueReference_ = valueReference;
}
    
double Variable::getMin() const{
  return (*this)->min_;
}

void Variable::setMin(double min){
  (*this)->min_ = min;
}
    
double Variable::getMax() const{
  return (*this)->max_;
}

void Variable::setMax(double max){
  (*this)->max_ = max;
}
    
double Variable::getNominal() const{
  return (*this)->nominal_;
}

void Variable::setNominal(double nominal){
  (*this)->nominal_ = nominal;
}
    
double Variable::getStart() const{
  return (*this)->start_;
}

void Variable::setStart(double start){
  (*this)->start_ = start;
}
    
const string& Variable::getUnit() const{
  return (*this)->unit_;
}

void Variable::setUnit(const string& unit){
  (*this)->unit_ = unit;
}
    
const string& Variable::getDisplayUnit() const{
  return (*this)->displayUnit_;
}

void Variable::setDisplayUnit(const string& displayUnit){
  (*this)->displayUnit_ = displayUnit;
}



} // namespace Modelica
} // namespace CasADi

