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

  /** \brief Holds expressions and meta-data corresponding to a physical quantity evolving in time
      \date 2012-2014
      \author Joel Andersson
   */
  struct Variable : public PrintableObject{
    
    /// Default constructor
    Variable();
    
    /// Variable name
    std::string name() const;

    /// Variable expression
    SX v; 

    /// Derivative expression
    SX d;            

    /// (Current) value of the variable
    double value;

    /// Nominal value
    double nominal;
    
    /// Value at time 0
    double start;

    /// Lower bound
    double min;

    /// Upper bound
    double max;

    /// Derivative at time 0
    double derivativeStart;

    /// Initial guess
    double initialGuess;

    /// Variability (see Fritzon)
    Variability variability;

    /// Causality (see Fritzon)
    Causality causality;

    /// Variable category
    Category category;

    /// Is the variable is an alias variable?
    Alias alias;
    
    /// Description
    std::string description;

    /// Variable reference (XML)
    int valueReference;
        
    /// Unit
    std::string unit;

    /// Display unit
    std::string displayUnit;
    
    /// Variable index
    int index;

    /// Free attribute
    bool free;

    /// Timed variable (never allocate)
    SX atTime(double t, bool allocate=false) const;

    /// Timed variable (allocate if necessary)
    SX atTime(double t, bool allocate=false);
        
  private:
#ifndef SWIG
    // Timed variables
    std::map<double,SX> timed_;

    // Print
    virtual void repr(std::ostream &stream) const;
    virtual void print(std::ostream &stream) const;
#endif // SWIG

  };
} // namespace CasADi

#endif // VARIABLE_HPP

