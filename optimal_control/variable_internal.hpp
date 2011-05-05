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
  namespace OptimalControl{

    /// Internal node class
  class VariableInternal : public SharedObjectNode{
    friend class Variable;
    public:

    // Constructor only available to the smart pointer class!
    VariableInternal(const std::string& name);
          
    // Destructor
    virtual ~VariableInternal();

    // Get name
    const std::string& getName() const;

    // Variable/binding equation
    SX var() const;  

    // Derivative/differential equation (never allocate)
    SX der(bool allocate) const;
    
    // Derivative/differential equation (allocate if necessary)
    SX der(bool allocate);
    
    // Timed variable (never allocate)
    SX atTime(double t, bool allocate) const;

    // Timed variable (allocate if necessary)
    SX atTime(double t, bool allocate);

    // Print
    virtual void repr(std::ostream &stream) const;
    virtual void print(std::ostream &stream) const;

    protected:

    // Attributes
    std::string name_;
    Variability variability_;
    Causality causality_;
    Alias alias_;
    std::string description_;
    int valueReference_;
    double min_;
    double max_;
    double nominal_;
    double start_;
    std::string unit_;
    std::string displayUnit_;
    
    // variable expression
    SX sx_; 

    // Derivative expression
    SX dx_;

    // Left hand side of the defining equation
    SX lhs_;
    
    // Right hand side of the defining equation
    SX rhs_;
        
    // Timed variables
    std::map<double, SX> timed_sx_;
    
    // Numerical value
    double val;
    
    // Index
    int index_;
    
    
  };
  
} // namespace OptimalControl
} // namespace CasADi


#endif // VARIABLE_INTERNAL_HPP