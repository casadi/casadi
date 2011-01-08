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
  
VariableInternal::VariableInternal(const string& name) : name_(name){
  // No expression by default
  sx_ = SX::nan;
  
  // Not differentable by default
  dx_ = SX::nan;
  
  // Binding equation undefined by default
  be_ = SX::nan;
  
  // Differential equation undefined by default
  de_ = SX::nan;

  independent_ = false;
  variability_ = CONTINUOUS;
  causality_ = INTERNAL;
  alias_ = NO_ALIAS;
  description_ = "";
  valueReference_ = -1; //?
  min_ = -numeric_limits<double>::infinity();
  max_ = numeric_limits<double>::infinity();
  nominal_ = 1.0;
  start_ = 0.0;
  unit_ = "";
  displayUnit_ = "";

  // Update the type
  init();
}

void VariableInternal::init(){
  // Set the type to unknown
  type_ = TYPE_UNKNOWN;
  
  // Try to determine the type
  if(independent_){
    type_ = TYPE_INDEPENDENT;
  } else {
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
}
  
VariableInternal::~VariableInternal(){
}

const string& VariableInternal::getName() const{
  return name_;
}

string VariableInternal::getTypeName() const{
  return typenames[type_];
}

SX VariableInternal::der() const{
  if(de_->isNan())
    return dx_;
  else
    return de_;
}

SX VariableInternal::sx() const{
  if(be_->isNan())
    return sx_;
  else
    return be_;
}

int VariableInternal::add(const Variable& var){
  return add(var,var.getName());
}

int VariableInternal::add(const Variable& var, const string& namepart){
  col.push_back(var);
  int ind = col.size()-1;
  name_part[namepart] = ind;
  return ind;
}

bool VariableInternal::has(const string& name) const{
  // try to locate the variable
  map<string, int>::const_iterator it = name_part.find(name);

  // check if the variable exists
  return it!=name_part.end();
}

void VariableInternal::repr(ostream &stream) const{
  stream << name_;
}


void VariableInternal::print(ostream &stream) const{
  print(stream,0);
}

void VariableInternal::print(ostream &stream, int indent) const{
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

void VariableInternal::getAll(vector<Variable>& vars) const{
  if(col.empty()){
    Variable temp;
    temp.assignNode(const_cast<VariableInternal*>(this)); // ugly trick
    vars.push_back(temp);
  } else {
    for(vector<Variable>::const_iterator it=col.begin(); it!=col.end(); ++it)
      if(!it->isNull())
        (*it)->getAll(vars);
  }
}

} // namespace Modelica
} // namespace CasADi

