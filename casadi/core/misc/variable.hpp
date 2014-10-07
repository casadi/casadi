/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010-2014 Joel Andersson, Joris Gillis, Moritz Diehl,
 *                            K.U. Leuven. All rights reserved.
 *    Copyright (C) 2011-2014 Greg Horn
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


#ifndef CASADI_VARIABLE_HPP
#define CASADI_VARIABLE_HPP

#include <iostream>
#include "../function/sx_function.hpp"
#include "../mx/mx.hpp"

namespace casadi {

  /// Time variability of a variable (see Fritzon page 89)
  enum Variability {CONSTANT, PARAMETER, DISCRETE, CONTINUOUS};

  /// Causality of a variable
  enum Causality {INPUT, OUTPUT, INTERNAL};

  /// Dynamics of the variable
  enum Dynamics {ALGEBRAIC, DIFFERENTIAL};

  /// Dynamics of the variable
  enum Alias {NO_ALIAS, ALIAS, NEGATED_ALIAS};

  /// Variable category
  enum Category {
    /** Unknown, not set */
    CAT_UNKNOWN,
    /** A state derivative */
    CAT_DERIVATIVE,
    /** A differential state, i.e. a variable that appears differentiated in the model */
    CAT_STATE,
    /** An independent constant: <tt>constant Real c1 = 3</tt> */
    CAT_DEPENDENT_CONSTANT,
    /** A dependent constant <tt>constant Real c2=c1*3</tt>. */
    CAT_INDEPENDENT_CONSTANT,
    /** A dependent parameter <tt>parameter Real p1=p2</tt>*/
    CAT_DEPENDENT_PARAMETER,
    /** An independent parameter <tt>parameter Real p2=3</tt>*/
    CAT_INDEPENDENT_PARAMETER,
    /** An algebraic variable or input */
    CAT_ALGEBRAIC
  };

  /** \brief Holds expressions and meta-data corresponding to a physical quantity evolving in time
      \date 2012-2014
      \author Joel Andersson
   */
  struct CASADI_CORE_EXPORT Variable : public PrintableObject<Variable> {

    /// Default constructor
    Variable();

    /// Variable name
    std::string name() const;

    /// Set the variable name (and corresponding expressions)
    void setName(const std::string& name);

    /// Variable expression
    SXElement v;

    /// Derivative expression
    SXElement d;

    /// Binding equation. Equal to "v" if unknown
    SXElement beq;

    /** \brief Derivative binding equation, i.e. ordinary differential equation (ODE).
     * Equal do "d" if unknown
     */
    SXElement ode;

    /// Nominal value
    double nominal;

    /// Value at time 0
    double start;

    /// Lower bound
    SXElement min;

    /// Upper bound
    SXElement max;

    /// Initial guess
    SXElement initialGuess;

    /// Derivative at time 0
    double derivativeStart;

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
    SXElement atTime(double t, bool allocate=false) const;

    /// Timed variable (allocate if necessary)
    SXElement atTime(double t, bool allocate=false);

    /// Print a description of the object
    void print(std::ostream &stream=std::cout, bool trailing_newline=true) const;

    /// Print a representation of the object
    void repr(std::ostream &stream=std::cout, bool trailing_newline=true) const;

#ifndef SWIG
  private:
    // Timed variables
    std::map<double, SXElement> timed_;
#endif // SWIG

  };
} // namespace casadi

#endif // CASADI_VARIABLE_HPP

