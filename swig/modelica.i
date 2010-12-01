%{
#include "modelica/fmi_parser.hpp"
#include "modelica/optimica_ocp.hpp"
#include "modelica/variable.hpp"
%}

namespace CasADi{
namespace Modelica{

    // Time variability of a variable (see Fritzon page 89)
    enum Variability{CONSTANT,PARAMETER,DISCRETE,CONTINUOUS};

    // Causality of a variable
    enum Causality{INPUT,OUTPUT,INTERNAL};
    
    // Dynamics of the variable
    enum Dynamics{ALGEBRAIC,DIFFERENTIAL};
    
    // Dynamics of the variable
    enum Alias{NO_ALIAS,ALIAS,NEGATED_ALIAS};
    
    // Variable types (REMOVE)
    enum VarType{TYPE_TIME,TYPE_STATE,TYPE_ALGEBRAIC,TYPE_CONTROL,TYPE_PARAMETER,TYPE_DEPENDENT,TYPE_DERIVATIVE, TYPE_NOT_SET};    

class Variable : public OptionsFunctionality{
  public:
    /** \brief Default constructor */
    Variable();
    
    /** \brief Construct an expression */
    explicit Variable(const std::string& name);
    
    /// Set and get value
    double getValue() const;
    void setValue(double val);    
    
    /** \brief Get the scalar expression */
    SX sx() const;
    
    /** \brief Time derivative */
    SX der() const;  
   
    /// Check if the variable has a certain type
    bool isTime() const;
    bool isDifferentialState() const;
    bool isAlgebraicState() const;
    bool isParameter() const;
    bool isControl() const;
    bool isDependent() const;

    const std::string& getName() const;
    void setName(const std::string& name);

    Variability getVariability() const;
    void setVariability(Variability variability);

    Causality getCausality() const;
    void setCausality(Causality causality);
    
    Alias getAlias() const;
    void setAlias(Alias alias);
    
    const std::string& getDescription() const;
    void setDescription(const std::string& description);
    
    int getValueReference() const;
    void setValueReference(int valueReference);
    
    double getMin() const;
    void setMin(double min);
    
    double getMax() const;
    void setMax(double max);
    
    double getNominal() const;
    void setNominal(double nominal);
    
    double getStart() const;
    void setStart(double start);
    
    const std::string& getUnit() const;
    void setUnit(const std::string& unit);
    
    const std::string& getDisplayUnit() const;
    void setDisplayUnit(const std::string& displayUnit);
    
    
};

%extend Variable {
  // Access a sub-collection by name
  Variable __getattr__(const std::string& name)  { return $self->operator()(name); }
  
  // Access a sub-collection by index
  Variable __getitem__(int ind){ return $self->operator[](ind); }
  
  // Print (why is this not inherited?)
  std::string __repr__()  { return $self->getRepresentation(); }
}



class OCP : public PrintableObject{
  public:    
    OCP();
    void sortVariables();
    void makeExplicit();
    void makeSemiExplicit();

    /// Variables in a class hierarchy
    Variable variables;

    
    
    
    // Time variable(s)
    std::vector<SX> t;
    
    // Differential states appearing implicitly
    std::vector<SX> x;

    // Time derivatiev of the differential states appearing implicitly
    std::vector<SX> xdot;
 
    // Differential states
    std::vector<SX> xd;
 
    // Algebraic states
    std::vector<SX> xa;
    
    // Controls
    std::vector<SX> u;
    
    // Parameters
    std::vector<SX> p;

    // Fully implicit equations
    std::vector<SX> dyneq;
    
    // Explicit differential equations
    std::vector<SX> diffeq;

    // Algebraic equations
    std::vector<SX> algeq;
    
    // Initial equations
    std::vector<SX> initeq;
    
    // Definition of dependent variables
    std::vector<SX> depdef;

    // OBJECTIVE
    // Mayer terms
    std::vector<SX> mterm;
    
    // Mayer time time point
    std::vector<double> mtp;
    
    // Constraint function with upper and lower bounds
    std::vector<SX> cfcn, cfcn_lb, cfcn_ub;

    // Initial time
    double t0;
    
    // Initial time is free
    bool t0_free;
    
    // Final time
    double tf;
    
    // Final time is free
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

