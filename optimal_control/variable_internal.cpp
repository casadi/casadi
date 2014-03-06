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
#include "../symbolic/casadi_exception.hpp"

using namespace std;
namespace CasADi{
  
  VariableInternal::VariableInternal(const string& name) : name_(name){
    // No expression by default
    var_ = casadi_limits<SXElement>::nan;
  
    // Not differentable by default
    der_ = casadi_limits<SXElement>::nan;
    
    // Not binding expressions by default
    binding_ = casadi_limits<SXElement>::nan;
    der_binding_ = casadi_limits<SXElement>::nan;
    
    variability_ = CONTINUOUS;
    causality_ = INTERNAL;
    category_ = CAT_UNKNOWN;
    alias_ = NO_ALIAS;
    description_ = "";
    valueReference_ = -1; //?
    min_ = -numeric_limits<double>::infinity();
    max_ = numeric_limits<double>::infinity();
    initial_guess_ = 0;
    nominal_ = 1.0;
    start_ = 0.0;
    derivative_start_ = 0.0;
    unit_ = "";
    displayUnit_ = "";
    free_ = false;
    is_differential_ = false;

    index_ = -1;
  }

  VariableInternal::~VariableInternal(){
  }

  const string& VariableInternal::getName() const{
    return name_;
  }

  void VariableInternal::repr(ostream &stream) const{
    stream << name_;
  }


  void VariableInternal::print(ostream &stream) const{
    stream << name_;
  }

  SXElement VariableInternal::atTime(double t, bool allocate) const{
    casadi_assert(!allocate);
    return const_cast<VariableInternal*>(this)->atTime(t,false);
  }

  SXElement VariableInternal::atTime(double t, bool allocate){
    // Find an existing element
    map<double,SXElement>::const_iterator it = timed_sx_.find(t);
  
    // If not found
    if(it==timed_sx_.end()){
      if(allocate){
        // Create a timed variable
        stringstream ss;
        ss << var_ << ".atTime(" << t << ")";
        SXElement tvar = SXElement::sym(ss.str());
      
        // Save to map
        timed_sx_[t] = tvar;
      
        // Return the expression
        return tvar;
      } else {
        casadi_error(" has no timed variable with t = " << t << ".");
      }
    
    } else {
      // Return the expression
      return it->second;
    }
  }

} // namespace CasADi

