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

#include "../symbolic/casadi_exception.hpp"
#include "variable.hpp"

using namespace std;
namespace CasADi{

  Variable::Variable(){
    this->variability = CONTINUOUS;
    this->causality = INTERNAL;
    this->category = CAT_UNKNOWN;
    this->alias = NO_ALIAS;
    this->description = "";
    this->valueReference = -1; //?
    this->value = numeric_limits<double>::quiet_NaN();
    this->min = -numeric_limits<double>::infinity();
    this->max = numeric_limits<double>::infinity();
    this->initialGuess = 0;
    this->nominal = 1.0;
    this->start = 0.0;
    this->derivativeStart = 0.0;
    this->unit = "";
    this->displayUnit = "";
    this->free = false;
    this->index = -1;
  }

  string Variable::name() const{
    return this->v.getName();
  }

  SX Variable::atTime(double t, bool allocate) const{
    casadi_assert(!allocate);
    return const_cast<Variable*>(this)->atTime(t,false);
  }

  SX Variable::atTime(double t, bool allocate){
    // Find an existing element
    map<double,SX>::const_iterator it = timed_.find(t);
  
    // If not found
    if(it==timed_.end()){
      if(allocate){
        // Create a timed variable
        stringstream ss;
        ss << name() << ".atTime(" << t << ")";
        SX tvar = SX::sym(ss.str());
      
        // Save to map
        timed_[t] = tvar;
      
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

  void Variable::repr(ostream &stream) const{
    stream << name();
  }

  void Variable::print(ostream &stream) const{
    stream << name();
  }

} // namespace CasADi
