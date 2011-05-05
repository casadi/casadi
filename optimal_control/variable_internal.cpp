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
namespace OptimalControl{
  
VariableInternal::VariableInternal(const string& name) : name_(name){
  // No expression by default
  sx_ = casadi_limits<SX>::nan;
  
  // Not differentable by default
  dx_ = casadi_limits<SX>::nan;
  
  // Equations undefined by default
  lhs_ = casadi_limits<SX>::nan;
  rhs_ = casadi_limits<SX>::nan;
  
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

  index_ = -1;
  
  // Update the type
  init();
}

void VariableInternal::init(){
  // Set the type to unknown
  type_ = TYPE_UNKNOWN;
  
  // Try to determine the type
  if(!sx_->isNan()){
    if(lhs_.isEqual(var())){
      type_ = TYPE_DEPENDENT;
    } else {
      if(variability_ == PARAMETER){
        type_ = TYPE_PARAMETER;
      } else if(variability_ == CONTINUOUS) {
        if(causality_ == INTERNAL){
          type_ = TYPE_STATE;
        } else if(causality_ == INPUT){
          type_ = TYPE_CONTROL;
        }
      } else if(variability_ == CONSTANT){
        type_ = TYPE_CONSTANT;
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

SX VariableInternal::der(bool allocate) const{
  casadi_assert(!allocate);
  return const_cast<VariableInternal*>(this)->der(false);
}

SX VariableInternal::der(bool allocate){
  if(dx_.isNan()){
    if(allocate){
      // Create derivative
      stringstream varname;
      varname << "der_" << var();
      dx_ = SX(varname.str());
    } else {
      stringstream msg;
      repr(msg);
      msg << " has no derivative expression";
      throw CasadiException(msg.str());
    }
  }
  return dx_;
}

SX VariableInternal::var() const{
  return sx_;
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

SX VariableInternal::atTime(double t, bool allocate) const{
  casadi_assert(!allocate);
  return const_cast<VariableInternal*>(this)->atTime(t,false);
}

SX VariableInternal::atTime(double t, bool allocate){
  // Find an existing element
  map<double,SX>::const_iterator it = timed_sx_.find(t);
  
  // If not found
  if(it==timed_sx_.end()){
    if(allocate){
      // Create a timed variable
      stringstream ss;
      ss << var() << ".atTime(" << t << ")";
      SX tvar(ss.str());
      
      // Save to map
      timed_sx_[t] = tvar;
      
      // Return the expression
      return tvar;
    } else {
      stringstream ss;
      repr(ss);
      ss << " has no timed variable with t = " << t << ".";
      throw CasadiException(ss.str());
    }
    
  } else {
    // Return the expression
    return it->second;
  }
}



} // namespace OptimalControl
} // namespace CasADi

