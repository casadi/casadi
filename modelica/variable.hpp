#ifndef VARIABLE_HPP
#define VARIABLE_HPP

#include <iostream>
#include "../casadi/fx/sx_function.hpp"
#include "../casadi/mx/mx.hpp"

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

    // Forward declaration
    class VariableNode;

  // Smart pointer class
  class Variable : public SharedObject{
    public:
    
    /** \brief Default constructor */
    Variable();

    /** \brief Construct an expression */
    explicit Variable(const std::string& name);
    
    /** \brief Set the type */
    void setType(VarType type);
    
    /** \brief  Destructor */
    ~Variable();
        
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
    
    // Access a sub-collection by name
    Variable operator()(const std::string& name) const;

    // Access a sub-collection by index
    Variable operator[](int ind) const;

    /** \brief Type conversion to SX */
    operator SX() const;

    /** \brief Type conversion to variable vector */
    operator std::vector<Variable>() const;
    
    /** \brief  Access functions of the node */
    VariableNode* operator->();
    const VariableNode* operator->() const;

    /** \brief  Getters and setters: Public data members instead? */
    //@{
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
    //@}
    
  };

  /// Internal node class
  class VariableNode : public SharedObjectNode{
    friend class Variable;
    public:

    // Constructor only available to the smart pointer class!
    VariableNode(const std::string& name);
          
    // Destructor
    virtual ~VariableNode();

    // Get type
    virtual std::string getType() const;
        
    // Get name
    virtual const std::string& getName() const;

    // Derivative        
    virtual SX der() const;  
    
    // Initialize
    virtual void init();    
    
    // Print
    virtual void repr(std::ostream &stream) const;
    virtual void print(std::ostream &stream) const;
    void print(std::ostream &stream, int indent) const;

    static const char* typenames[];
        
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
    
    SX sx_; // variable expression

    VarType type;
    
    // Derivative
    Variable d;
    
    // Numerical value
    double val;
    
    // Set of sub-collections in the current collection
    std::vector<Variable> col;
  
    // Names of contained collections
    std::map<std::string,int> name_part;

    // Get all the variables
    void getAll(std::vector<Variable>& vars) const;
    
    // Add a subcollection
    int add(const Variable& var, const std::string& namepart);

    // Check if a subcollection exists
    bool has(const std::string& name) const;

    
  };
  
} // namespace Modelica
} // namespace CasADi


#endif // VARIABLE_HPP