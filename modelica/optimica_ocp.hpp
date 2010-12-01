#ifndef OPTIMICA_OCP_HPP
#define OPTIMICA_OCP_HPP

#include "casadi/printable_object.hpp"
#include "variable.hpp"

namespace CasADi{
  namespace Modelica{

/** Symbolic, object oriented representation of an optimal control problem (OCP) */
class OCP : public PrintableObject{
  public:    
    /** \brief OCP */
    OCP();

    /// Print a destription of the object
    virtual void print(std::ostream &stream=std::cout) const;

    /// Sort variables
    void sortVariables();

    /// Try to make explicit by symbolically solving for xdot (experimental, only small systems)
    void makeExplicit();

    /// Replace all state derivatives by algebraic variables with the same name
    void makeSemiExplicit();

    /// Access the variables in a class hierarchy -- public data member
    Variable variables;
    
  public:

    /// Time variable(s) TODO: Should not be a vector
    std::vector<SX> t;
    
    /// Differential states appearing implicitly
    std::vector<SX> x;

    /// Time derivative of the differential states appearing implicitly
    std::vector<SX> xdot;
 
    /// Differential states
    std::vector<SX> xd;
 
    /// Algebraic states
    std::vector<SX> xa;
    
    /// Controls
    std::vector<SX> u;
    
    /// Parameters
    std::vector<SX> p;

    /// Dependent variables and constants
    std::vector<Variable> d;
    
    // EQUATIONS

    /// Fully implicit equations
    std::vector<SX> dyneq;
    
    /// Explicit differential equations
    std::vector<SX> diffeq;

    /// Algebraic equations
    std::vector<SX> algeq;
    
    /// Initial equations
    std::vector<SX> initeq;
    
    /// Definition of dependent variables
    std::vector<SX> depdef;

    // OBJECTIVE
    /// Mayer terms
    std::vector<SX> mterm;
    
    /// Mayer time time point
    std::vector<double> mtp;
    
    /// Constraint function with upper and lower bounds
    std::vector<SX> cfcn, cfcn_lb, cfcn_ub;

    /// Lagrange terms (symbolic/numeric)
    // std::vector<LagrangeTerm> lterm;
    
    /// Least squares terms (symbolic/numeric)
    // std::vector<LagrangeTerm> lsqterm;
    
    /// Initial time
    double t0;
    
    /// Initial time is free
    bool t0_free;
    
    /// Final time
    double tf;
    
    /// Final time is free
    bool tf_free;

};

  } // namespace Modelica
} // namespace CasADi

#endif // OPTIMICA_OCP_HPP


