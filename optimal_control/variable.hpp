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

#ifndef VARIABLE_HPP
#define VARIABLE_HPP

#include <iostream>
#include "../symbolic/fx/sx_function.hpp"
#include "../symbolic/mx/mx.hpp"

namespace CasADi{
    
  /// Time variability of a variable (see Fritzon page 89)
  enum Variability{CONSTANT,PARAMETER,DISCRETE,CONTINUOUS};

  /// Causality of a variable
  enum Causality{INPUT,OUTPUT,INTERNAL};
    
  /// Dynamics of the variable
  enum Dynamics{ALGEBRAIC,DIFFERENTIAL};
    
  /// Dynamics of the variable
  enum Alias{NO_ALIAS,ALIAS,NEGATED_ALIAS};
    
  /// Variable category
  enum Category{
    /** Unknown, not set */
    CAT_UNKNOWN,
    /** A state derivative */
    CAT_DERIVATIVE,
    /** A differential state, i.e. a variable that appears differentiated in the model */
    CAT_STATE, 
    /** An independent constant: "constant Real c1 = 3" */
    CAT_DEPENDENT_CONSTANT,
    /** A dependent constant "constant Real c2=c1*3". */
    CAT_INDEPENDENT_CONSTANT,
    /** A dependent parameter "parameter Real p1=p2"*/
    CAT_DEPENDENT_PARAMETER,
    /** An independent parameter "parameter Real p2=3"*/
    CAT_INDEPENDENT_PARAMETER,
    /** An algebraic variabel or input */
    CAT_ALGEBRAIC
  };

  // Forward declaration
  class VariableInternal;

  /** \brief Smart pointer class to a Variable
   *  
   *  A Variable is an SX expression with meta-data attached.
   */
  class Variable : public SharedObject{
  public:
    
    /// Default (empty) constructor
    Variable();

    /// Create a new variable
    explicit Variable(const std::string& name);
    
    /// Destructor
    virtual ~Variable();
        
    /// Variable name
    std::string name() const;

    /// Variable expression
    
    /// Differential expression

    /// Timed variable (never allocate)
    SX atTime(double t, bool allocate=false) const;

    /// Timed variable (allocate if necessary)
    SX atTime(double t, bool allocate=false);
    
    /// Variable index
        
    /// Access functions of the node
    VariableInternal* operator->();

    /// Const access functions of the node
    const VariableInternal* operator->() const;

    /// Variability (see Fritzon)

    /// Causality (see Fritzon)

    /// Variable category

    /// Is the variable is an alias variable?
    
    /// Description

    /// Variable reference (XML)
    
    /// Nominal value

    /// Value at time 0

    /// Lower bound

    /// Upper bound

    /// Derivative at time 0

    /// Unit

    /// Display unit
    
    /// Expression

    /// Derivative expression
                
    /// Variable index

    /// Free attribute
        
    /// Check if the node is pointing to the right type of object
    virtual bool checkNode() const;
  };
} // namespace CasADi

#ifdef SWIG
// Template instantiations
%template(VariableVector) std::vector<CasADi::Variable>;
#endif // SWIG  


#endif // VARIABLE_HPP

