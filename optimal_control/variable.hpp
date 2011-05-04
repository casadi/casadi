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
#include "../casadi/fx/sx_function.hpp"
#include "../casadi/mx/mx.hpp"

namespace CasADi{
  namespace OptimalControl{
    
    /// Time variability of a variable (see Fritzon page 89)
    enum Variability{CONSTANT,PARAMETER,DISCRETE,CONTINUOUS};

    /// Causality of a variable
    enum Causality{INPUT,OUTPUT,INTERNAL};
    
    /// Dynamics of the variable
    enum Dynamics{ALGEBRAIC,DIFFERENTIAL};
    
    /// Dynamics of the variable
    enum Alias{NO_ALIAS,ALIAS,NEGATED_ALIAS};
    
    /// Variable types
    enum VarType{TYPE_STATE,TYPE_ALGEBRAIC,TYPE_CONTROL,TYPE_PARAMETER,TYPE_CONSTANT,TYPE_DEPENDENT,TYPE_UNKNOWN};    

    // Forward declaration
    class VariableInternal;

  // Smart pointer class
  class Variable : public SharedObject{
    public:
    
    /// Default (empty) constructor
    Variable();

    /// Create a new variable
    explicit Variable(const std::string& name, bool create_expression = true);
    
    /// Destructor
    virtual ~Variable();
        
    /// Get the variable expression
    SX var() const;
    
    /// Get differential expression (never allocate)
    SX der(bool allocate=false) const;
    
    /// Get differential expression (allocate if necessary)
    SX der(bool allocate=false);
    
    /// Left hand side of the defining (differential) equation
    SX lhs() const;
    
    /// Right hand side of the defining (differential) equation
    SX rhs() const;
    
    /// Timed variable (never allocate)
    SX atTime(double t, bool allocate=false) const;

    /// Timed variable (allocate if necessary)
    SX atTime(double t, bool allocate=false);
    
    /// Get the variable index
    int index() const;
    
    /// Access a sub-collection by name
    Variable operator()(const std::string& name) const;

#ifndef SWIG
    /// Access a sub-collection by index
    Variable operator[](int ind) const;

    /// Type conversion to SX
    //operator SX() const;

    /// Type conversion to variable vector
    operator std::vector<Variable>() const;
#endif // SWIG
    
    /// Access functions of the node
    VariableInternal* operator->();

    /// Const access functions of the node
    const VariableInternal* operator->() const;

    /// Get type
    VarType getType() const;
    
    /// Get numerical value
    double getValue() const;

    /// Set numerical value
    void setValue(double val);    

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
    
    /// Get the upper bound
    double getMax() const;

    /// Set the upper bound
    void setMax(double max);
    
    /// Get the nominal value of the variable
    double getNominal() const;

    /// Set the nominal value of the variable
    void setNominal(double nominal);
    
    /// Get the value at time 0
    double getStart() const;

    /// Set the value at time 0
    void setStart(double start);
    
    /// Get the unit
    const std::string& getUnit() const;

    /// Set the unit
    void setUnit(const std::string& unit);
    
    /// Get the display unit
    const std::string& getDisplayUnit() const;

    /// Set the display unit
    void setDisplayUnit(const std::string& displayUnit);
    
    /// Set the expression
    void setExpression(const SX& sx);

    /// Set the derivative expression
    void setDerivative(const SX& dx);
    
    /// Set the left hand side of the defining equation
    void setLHS(const SX& ex);
    
    /// Set the right hand side of the defining equation
    void setRHS(const SX& ex);
    
    /// Set equation
    void setEquation(const SX& l, const SX& r);
            
    /// Set the variable index
    void setIndex(int ind);
    
    /// Check if the node is pointing to the right type of object
    virtual bool checkNode() const;
  };
} // namespace OptimalControl
} // namespace CasADi

  
#ifdef SWIG  
%extend CasADi::OptimalControl::Variable {
  // Not inherited: Bug?
  std::string __repr__()  { return $self->getRepresentation(); }
  
  // Not automatically translated
  CasADi::OptimalControl::Variable __getitem__(int ind) const { return (*$self)[ind];}
}

// Template instantiations
%template(vector_variable) std::vector<CasADi::OptimalControl::Variable>;
#endif // SWIG  


#endif // VARIABLE_HPP