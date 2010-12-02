/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A minimalistic computer algebra system with automatic differentiation 
 *    and framework for dynamic optimization.
 *    Copyright (C) 2010 by Joel Andersson, Moritz Diehl et al., K.U.Leuven. All rights reserved.
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

#include "variable.hpp"
#include <cassert>
#include "../casadi/expression_tools.hpp"
#include "../casadi/casadi_exception.hpp"

using namespace std;
namespace CasADi{
namespace Modelica{
  
Variable::Variable(){
}

Variable::Variable(const string& name){
  assignNode(new VariableNode(name));
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
  // No expression by default
  sx_ = SX::nan;
  
  // Not differentable by default
  dx_ = SX::nan;
  
  // Binding equation undefined by default
  be_ = SX::nan;
  
  // Differential equation undefined by default
  de_ = SX::nan;

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

  // Update the type
  init();
}

void VariableNode::init(){
  // Set the type to unknown
  type_ = TYPE_UNKNOWN;
  
  // Try to determine the type
  if(!sx_->isNan()){
    if(!be_->isNan()){
      type_ = TYPE_DEPENDENT;
    } else {
      if(variability_ == PARAMETER){
        type_ = TYPE_PARAMETER;
      } else if(variability_ == CONTINUOUS) {
        if(causality_ == INTERNAL){
          type_ = !dx_->isNan() ? TYPE_STATE : TYPE_ALGEBRAIC;
        } else if(causality_ == INPUT){
          type_ = TYPE_CONTROL;
        }
      } else if(variability_ == CONSTANT){
        type_ = TYPE_CONSTANT;
      }
    }
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
  return (*this)->sx();  
}



const string& VariableNode::getName() const{
  return name_;
}

string VariableNode::getTypeName() const{
  return typenames[type_];
}

VarType Variable::getType() const{
  return (*this)->type_;
}

SX VariableNode::der() const{
  if(de_->isNan())
    return dx_;
  else
    return de_;
}

SX VariableNode::sx() const{
  if(be_->isNan())
    return sx_;
  else
    return be_;
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
  stream << name_ << ": " << getTypeName() << endl;

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
  init();
}

Causality Variable::getCausality() const{
  return (*this)->causality_;
}

void Variable::setCausality(Causality causality){
  (*this)->causality_ = causality;
  init();
}
    
Alias Variable::getAlias() const{
  return (*this)->alias_;
}

void Variable::setAlias(Alias alias){
  (*this)->alias_ = alias;
  init();
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

void Variable::setExpression(const SX& sx){
  (*this)->sx_ = sx;
  init();
}

const SX& Variable::getExpression() const{
  return (*this)->sx_;
}

void Variable::setDerivative(const SX& dx){
  (*this)->dx_ = dx;
  init();
}

const SX& Variable::getDerivative() const{
  return (*this)->dx_;
}

void Variable::setBindingEquation(const SX& be){
  (*this)->be_ = be;
  init();
}
    
const SX& Variable::getBindingEquation() const{
  return (*this)->be_;
}
    
void Variable::setDifferentialEquation(const SX& de){
  (*this)->de_ = de;
  init();
}
    
const SX& Variable::getDifferentialEquation() const{
  return (*this)->de_;
}



} // namespace Modelica
} // namespace CasADi

