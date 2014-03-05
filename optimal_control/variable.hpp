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
  *  In a sense, a Variable is an SX expression with meta-data attached.
  */
  class Variable : public SharedObject{
    public:
    
    /// Default (empty) constructor
    Variable();

    /// Create a new variable
    explicit Variable(const std::string& name, bool create_expression = true);
    
    /// Destructor
    virtual ~Variable();
    
    /// Get the variable expression
    SXElement var() const;
    
    /// Get differential expression
    SXElement der() const;
    
    /// Get the binding expression for the variable or its derivative
    SXElement binding(bool derivative=false) const;
    
    /// Get the highest order derivative (i.e. der() or var())
    SXElement highest() const;
    
    /// Timed variable (never allocate)
    SXElement atTime(double t, bool allocate=false) const;

    /// Timed variable (allocate if necessary)
    SXElement atTime(double t, bool allocate=false);
    
    /// Get the variable index
    int index() const;
        
    /// Access functions of the node
    VariableInternal* operator->();

    /// Const access functions of the node
    const VariableInternal* operator->() const;

    /// Get variable name
    const std::string& getName() const;

    /// Set variable name
    void setName(const std::string& name);

    /// Get the variability (see Fritzon)
    Variability getVariability() const;

    /// Set the variability (see Fritzon)
    void setVariability(Variability variability);

    /// Get the causality (see Fritzon)
    Causality getCausality() const;

    /// Set the causality (see Fritzon)
    void setCausality(Causality causality);

    /// Get the variable category
    Category getCategory() const;

    /// Set the variable category
    void setCategory(Category category);

    /// Check if the variable is an alias variable
    Alias getAlias() const;

    /// Set if the variable is an alias variable
    void setAlias(Alias alias);
    
    /// Get the description
    const std::string& getDescription() const;
    
    /// Set the description
    void setDescription(const std::string& description);
    
    /// Get the variable reference (XML)
    int getValueReference() const;

    /// Set the variable reference (XML)
    void setValueReference(int valueReference);
    
    /// Get the lower bound
    double getMin() const;

    /// Set the lower bound
    void setMin(double min);
    
    /// Access the lower bound
    double& min();

    /// Get the upper bound
    double getMax() const;

    /// Set the upper bound
    void setMax(double max);
    
    /// Access the upper bound
    double& max();

    /// Get the nominal value of the variable
    double getNominal() const;

    /// Set the nominal value of the variable
    void setNominal(double nominal);
    
    /// Access the nominal value of the variable
    double& nominal();

    /// Get the value at time 0
    double getStart() const;

    /// Set the value at time 0
    void setStart(double start);

    /// Access the value at time 0
    double& start();
        
    /// Get the lower bound
    double getInitialGuess() const;

    /// Set the lower bound
    void setInitialGuess(double initial_guess);

    /// Access the lower bound
    double& initialGuess();

    /// Get the derivative at time 0
    double getDerivativeStart() const;

    /// Set the derivative at time 0
    void setDerivativeStart(double start);
    
    /// Access the derivative at time 0
    double& derivativeStart();

    /// Get the unit
    const std::string& getUnit() const;

    /// Set the unit
    void setUnit(const std::string& unit);
    
    /// Access the unit
    std::string& unit();
    
    /// Get the display unit
    const std::string& getDisplayUnit() const;

    /// Set the display unit
    void setDisplayUnit(const std::string& displayUnit);

    /// Get the display unit
    std::string& displayUnit();
    
    /// Set the expression
    void setExpression(const SXElement& v);

    /// Set the derivative expression
    void setDerivative(const SXElement& d);
                
    /// Set the binding expression for the variable or its derivative
    void setBinding(const SXElement& binding, bool derivative=false);
                
    /// Set the variable index
    void setIndex(int ind);
    
    /// Is differential?
    bool isDifferential() const;
    
    /// Set differential
    void setDifferential(bool is_differential);
    
    /// Get the the free attribute
    bool getFree() const;

    /// Set the the free attribute
    void setFree(bool free);

    /// Access the the free attribute
    bool& free();
        
    /// Check if the node is pointing to the right type of object
    virtual bool checkNode() const;
  };
} // namespace CasADi

#ifdef SWIG
// Template instantiations
%template(VariableVector) std::vector<CasADi::Variable>;
#endif // SWIG  


#endif // VARIABLE_HPP

