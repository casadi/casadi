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

#ifndef VARIABLE_INTERNAL_HPP
#define VARIABLE_INTERNAL_HPP

#include "variable.hpp"

namespace CasADi{

    /// Internal node class
  class VariableInternal : public SharedObjectNode{
    friend class Variable;
    public:

      // Constructor only available to the smart pointer class!
      VariableInternal(const std::string& name);
            
      // Destructor
      virtual ~VariableInternal();

      // Clone
      virtual VariableInternal* clone() const{ return new VariableInternal(*this);}

      // Timed variable (never allocate)
      SX atTime(double t, bool allocate) const;

      // Timed variable (allocate if necessary)
      SX atTime(double t, bool allocate);

      // Print
      virtual void repr(std::ostream &stream) const;
      virtual void print(std::ostream &stream) const;

    protected:
      
      // Attributes
      Variability variability_;
      Causality causality_;
      Category category_;
      Alias alias_;
      std::string description_;
      int valueReference_;
      bool free_;
      
      double min_, max_, nominal_, start_, derivative_start_, initial_guess_;
      std::string unit_, displayUnit_;
      
      // variable expression
      SX var_; 

      // Derivative expression
      SX der_;
          
      // Binding expression
      SX binding_;
          
      // Binding expression for the derivative
      SX der_binding_;
          
      // Timed variables
      std::map<double,SX> timed_sx_;
            
      // Index
      int index_;
      
      // Does the expression appear differentiated
      bool is_differential_;
  };
  
  
} // namespace CasADi


#endif // VARIABLE_INTERNAL_HPP
