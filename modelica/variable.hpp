/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A minimalistic computer algebra system with automatic differentiation 
 *    and framework for dynamic optimization.
 *    Copyright (C) 2010 by Joel Andersson, Moritz Diehl et al., K.U.Leuven. All rights reserved.
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
  namespace Modelica{
    
    /// Time variability of a variable (see Fritzon page 89)
    enum Variability{CONSTANT,PARAMETER,DISCRETE,CONTINUOUS};

    /// Causality of a variable
    enum Causality{INPUT,OUTPUT,INTERNAL};
    
    /// Dynamics of the variable
    enum Dynamics{ALGEBRAIC,DIFFERENTIAL};
    
    /// Dynamics of the variable
    enum Alias{NO_ALIAS,ALIAS,NEGATED_ALIAS};
    
    /// Variable types
    enum VarType{TYPE_INDEPENDENT, TYPE_STATE,TYPE_ALGEBRAIC,TYPE_CONTROL,TYPE_PARAMETER,TYPE_CONSTANT,TYPE_DEPENDENT,TYPE_UNKNOWN};    

    /// Names of the variable types
    static const char* typenames[] = {"INDEPENDENT","STATE","ALGEBRAIC","CONTROL","PARAMETER","CONSTANT","DEPENDENT","UNKNOWN"};

    
    // Forward declaration
    class VariableNode;

  // Smart pointer class
  class Variable : public SharedObject{
    public:
    
    /// Default (empty) constructor
    Variable();

    /// Create a new variable
    explicit Variable(const std::string& name);
    
    /// Destructor
    ~Variable();
        
    /// Get the scalar expression or binding expression
    SX sx() const;
    
    /// Get the time derivative or differential equation
    SX der() const;  
       
    /// Access a sub-collection by name
    Variable operator()(const std::string& name) const;

    /// Access a sub-collection by index
    Variable operator[](int ind) const;

    /// Type conversion to SX
    operator SX() const;

    /// Type conversion to variable vector
    operator std::vector<Variable>() const;
    
    /// Access functions of the node
    VariableNode* operator->();
    const VariableNode* operator->() const;

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

    /// Get the expression
    const SX& getExpression() const;

    /// Set the derivative expression
    void setDerivative(const SX& dx);
    
    /// Get the derivative expression
    const SX& getDerivative() const;
    
    /// Set the binding equation
    void setBindingEquation(const SX& be);
    
    /// Get the binding equation
    const SX& getBindingEquation() const;
    
    /// Set the differential equation
    void setDifferentialEquation(const SX& de);
    
    /// Get the differential equation
    const SX& getDifferentialEquation() const;
    
    /// Mark the variable as independent/not independent (time)
    void setIndependent(bool independent);
    
    /// Check if variable is independent (time)
    bool getIndependent() const;
    
  };

  /// Internal node class
  class VariableNode : public SharedObjectNode{
    friend class Variable;
    public:

    // Constructor only available to the smart pointer class!
    VariableNode(const std::string& name);
          
    // Destructor
    virtual ~VariableNode();

    // Get type name
    std::string getTypeName() const;
        
    // Get name
    const std::string& getName() const;

    // Variable/binding equation
    SX sx() const;  

    // Derivative/differential equation
    SX der() const;  
    
    // Update the type
    virtual void init();
    
    // Print
    virtual void repr(std::ostream &stream) const;
    virtual void print(std::ostream &stream) const;
    void print(std::ostream &stream, int indent) const;

    // Get all the variables
    void getAll(std::vector<Variable>& vars) const;
    
    // Add a subcollection
    int add(const Variable& var, const std::string& namepart);

    // Add a subcollection, default name
    int add(const Variable& var);

    // Check if a subcollection exists
    bool has(const std::string& name) const;

    // Set of sub-collections in the current collection
    std::vector<Variable> col;
  
    // Names of contained collections
    std::map<std::string,int> name_part;

    protected:

    // Variable type
    VarType type_;
    
    // Is the variable independent
    bool independent_;

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

    // Binding equation
    SX be_;
    
    // Differential equation
    SX de_;
    
    // Numerical value
    double val;
    


    
  };
  
} // namespace Modelica
} // namespace CasADi


#endif // VARIABLE_HPP