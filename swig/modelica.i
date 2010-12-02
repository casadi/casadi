%{
#include "modelica/fmi_parser.hpp"
#include "modelica/optimica_ocp.hpp"
#include "modelica/variable.hpp"
%}

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

class Variable : public OptionsFunctionality{
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

%extend Variable {
  // Access a sub-collection by name
  Variable __getattr__(const std::string& name)  { return $self->operator()(name); }
  
  // Access a sub-collection by index
  Variable __getitem__(int ind){ return $self->operator[](ind); }
  
  // Print (why is this not inherited?)
  std::string __repr__()  { return $self->getRepresentation(); }
}

} // namespace Modelica
} // namespace CasADi

// Template instantiations
namespace std {
%template(vector_variable) vector<CasADi::Modelica::Variable>;
} // namespace std;

namespace CasADi{
  namespace Modelica{

/** Symbolic, object oriented representation of an optimal control problem (OCP) */
class OCPVariables : public PrintableObject{
  public:
    /// Constructor (automatic type conversion allowed)
    OCPVariables(const Variable& var);

    /// Time
    Variable t;
    
    /// Differential states
    std::vector<Variable> x;

    /// Algebraic states
    std::vector<Variable> z;
    
    /// Controls
    std::vector<Variable> u;
    
    /// Free parameters
    std::vector<Variable> p;

    /// Constants
    std::vector<Variable> c;

    /// Dependent
    std::vector<Variable> d;
};

%extend OCPVariables {
  // Print (why is this not inherited?)
  std::string __repr__()  { return $self->getRepresentation(); }
}

class OCP : public PrintableObject{
  public:    
    /// OCP
    OCP();

    /// Access the variables in a class hierarchy -- public data member
    Variable variables;

    /// Differential algebraic equations
    std::vector<SX> dae;
    
    /// Algebraic equations
    std::vector<SX> ae;

    /// Initial equations
    std::vector<SX> initeq;

    /// Path constraint function with upper and lower bounds
    std::vector<SX> cfcn, cfcn_lb, cfcn_ub;

    /// Mayer terms
    std::vector<SX> mterm;
    
    /// Mayer time time points
    std::vector<double> mtp;
        
    /// Initial time
    double t0;
    
    /// Initial time is free
    bool t0_free;
    
    /// Final time
    double tf;
    
    /// Final time is free
    bool tf_free;
};

%extend OCP{
  // Print (why is this not inherited?)
  std::string __repr__()  { return $self->getRepresentation(); }
}


class FMIParser{
public:
  FMIParser(const std::string& filename);    // constructor
  virtual ~FMIParser(); // destructor
  OCP& parse();
};

%extend FMIParser {
std::string __str__()  { return $self->getDescription(); }


}


} // namespace Modelica
} // namespace CasADi

