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

#ifndef FLAT_OCP_HPP
#define FLAT_OCP_HPP

#include "variable.hpp"

namespace CasADi{
namespace OptimalControl{

// Forward declaration
class FlatOCPInternal;

/** \brief A flat OCP representation coupled to an XML file

 <H3>Variables:  </H3>
  \verbatim
   t :     time
   x :     differential and algebraic states defined by a fully-implicit DAE
   xd:     differential states defined by an explicit ODE
   xa:     algebraic states defined by an algebraic equation
   q :     quadrature states
   y :     dependent variables
   p :     independent parameters
   u :     control signals
  \endverbatim 
  
  <H3>Equations:  </H3>
  \verbatim
  fully implicit DAE:       0 = dae(t,x,\dot{x},xd,xa,u,p)
  explicit ODE:      \dot{xd} = ode(t,x,xd,xa,u,p)
  quadratures:        \dot{q} = quad(t,x,xd,xa,u,p)
  algebraic equations:      0 = alg(t,x,xd,xa,u,p)
  dependent equations:      y = dep(t,x,xd,xa,u,p)
  initial equations:        0 = initial(t,x,\dot{x},xd,xa,u,p)
  \endverbatim 

  Note that when parsed, all dynamic states, differential and algebraic, end up in the category "x" 
  and all dynamic equations end up in the implicit category "dae". At a later state, the DAE can be
  reformulated, for example in semi-explicit form, possibly in addition to a set of quadrature states.
 
  Also note that division of the states into three categories for states defined by a DAE, states
  defined by an ODE and states defined by an algebraic equation. The category "xd" does thus _not_
  include differential states that are implicitly defined by the DAE.

  The functions for reformulation is are provided as member functions to this class or as independent
  functions located in the header file "ocp_tools.hpp".

  <H3>Usage skeleton:</H3>
  
  1. Call constructor with an FMI conformant XML file or pass an empty string ("") to build up the OCP from scratch <BR>
  > FlatOCP ocp(xml_file_name)
  
  2. Set options <BR>
  > ocp.setOption(...,...)

  3. Initialize and parse XML <BR>
  > ocp.init()

  4. Modify/add variables, equations, optimization <BR>
  > ...
  
  When the optimal control problem is in a suitable form, it is possible to either generate functions
  for numeric/symbolic evaluation or exporting the OCP formulation into a new FMI conformant XML file.
  The latter functionality is not yet available.

  \date 2011
  \author Joel Andersson
*/
class FlatOCP : public OptionsFunctionality{
  public:
    /// Default (empty) constructor
    FlatOCP();
    
    /// Create an FMI parser instance given the filename
    FlatOCP(const std::string& filename);

    /// Parse from XML to C++ format
    void parse();

    /// Add a variable
    void addVariable(const std::string& name, const Variable& var);
    
    /// Access a variable by name
    Variable& variable(const std::string& name);

    /// Make a differential state algebraic by replacing its time derivative by 0
    void makeAlgebraic(const std::string& name);
    
    /// Access to the internal class
    FlatOCPInternal* operator->();

    /// Const access to the internal class
    const FlatOCPInternal* operator->() const;

    /// Check if the node is pointing to the right type of object
    virtual bool checkNode() const;

    /** @name Variables
    *  Get all variables of a particular type 
    */
    //@{
    /// Time
    SX t() const;
    
    /// Differential and algebraic states defined by a fully-implicit DAE (length == dae().size())
    std::vector<Variable>& x();
    
    /// Differential states defined by an explicit ODE (length == ode().size())
    std::vector<Variable>& xd();
    
    /// Algebraic states defined by an algebraic equation (length == alg().size())
    std::vector<Variable>& xa();
    
    /// All states, differential and algebraic (includes x, xd and xa)
    std::vector<Variable> x_all() const;
    
    /// Quadrature states (length == quad().size())
    std::vector<Variable>& q();
    
    /// Dependent variables (length == dep().size())
    std::vector<Variable>& y();
    
    /// Independent parameters
    std::vector<Variable>& p();
    
    /// Control signals
    std::vector<Variable>& u();
    //@}
    
    /** @name Equations
    *  Get all equations of a particular type 
    */
    //@{
    /// Fully implicit DAE (length == x().size())
    std::vector<SX>& dae();
    
    /// Explicit ODE  (length == xd().size())
    std::vector<SX>& ode();
    
    /// Algebraic equations (length == xa().size())
    std::vector<SX>& alg();
    
    /// Quadrature states (length == q().size())
    std::vector<SX>& quad();
    
    /// Dependent equations (length == y().size())
    std::vector<SX>& dep();
    
    /// Initial equations (remove?)
    std::vector<SX>& initial();
    //@}

    /** @name Optimization
    *  Formulate the dynamic optimization problem. Note that the variable bounds are located inside
       the respective variable.
    */
    //@{
    /// Mayer terms in the objective
    std::vector<SX>& mterm();
    
    /// Lagrange terms in the objective
    std::vector<SX>& lterm();

    /// Path constraint functions
    std::vector<SX>& path();
    
    /// Path constraint functions upper bounds
    std::vector<double>& path_min();

    /// Path constraint functions upper bounds
    std::vector<double>& path_max();

    /// Interval start time
    double t0() const;
    
    /// Interval final time
    double tf() const;
    
    /// Interval start time is free
    bool t0_free() const;
    
    /// Interval final time is free
    bool tf_free() const;
    
    /// Set interval start time
    void set_t0(double t);
    
    /// Set interval final time
    void set_tf(double t);
    
    /// Set interval start time as free or fixed
    void set_t0_free(bool free);
    
    /// Set interval final time as free or fixed
    void set_tf_free(bool free);
    //@}
    
    /** @name Manipulation
    *  Reformulate the dynamic optimization problem.
    */
    //@{
    
    /// Eliminate interdependencies in the dependent equations
    void eliminateInterdependencies();
    
    /// Eliminate dependent equations, by default sparing the dependent variables with upper or lower bounds
    void eliminateDependent();

    /// Sort the DAE equations and variables
    void sortDAE();

    /// Transform the fully implicit DAE to a explicit or semi-explicit form
    void makeExplicit();

    /// Get the DAE input arguments
    std::vector<SXMatrix> daeArg() const;
    
    /// Substitute the dependents from a set of expressions
    std::vector<SXMatrix> substituteDependents(const std::vector<SXMatrix>& x) const;
    
    /** \brief Get the ODE/DAE right hand side function
    * The returned FX has the following input/output scheme:
    * @copydoc scheme_DAEInput 
    * @copydoc scheme_DAEOutput
    */
    FX daeFcn() const;
    
    /// Generate a MUSCOD-II compatible DAT file
    void generateMuscodDatFile(const std::string& filename, const Dictionary& mc2_ops=Dictionary()) const;
    
    //@}
};

} // namespace OptimalControl
} // namespace CasADi

#endif //FLAT_OCP_HPP
