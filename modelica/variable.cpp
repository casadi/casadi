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

#include "variable_internal.hpp"
#include "../casadi/casadi_exception.hpp"

using namespace std;
namespace CasADi{
namespace Modelica{
  
Variable::Variable(){
}

Variable::Variable(const string& name){
  assignNode(new VariableInternal(name));
}


Variable::~Variable(){
}

VariableInternal* Variable::operator->(){
  return (VariableInternal*)(SharedObject::operator->());
}

const VariableInternal* Variable::operator->() const{
  return (const VariableInternal*)(SharedObject::operator->());
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

VarType Variable::getType() const{
  return (*this)->type_;
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

void Variable::setIndependent(bool independent){
  (*this)->independent_ = independent;
  init();
}
    
bool Variable::getIndependent() const{
  return (*this)->independent_;
}

bool Variable::checkNode() const{
  return dynamic_cast<const VariableInternal*>(get());
}

} // namespace Modelica
} // namespace CasADi

