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
*
*  <H3>Variables:  </H3>
*    \f$ t \f$:     time <BR>
*    \f$ s \f$:     implicitly defined states (differential or algebraic) <BR>
*   \f$ x \f$:     differential states <BR>
*    \f$ z \f$:     algebraic states <BR>
*    \f$ q \f$:     quadrature states <BR>
*    \f$ y \f$:     dependent variables <BR>
*    \f$ p \f$:     independent parameters <BR>
*    \f$ u \f$:     control signals <BR>
*  
*  <H3>Equations:  </H3>
*  fully implicit DAE: \f$ 0 = \text{dae}(t,s,\dot{s},x,z,u,p)  \f$  <BR>
*  explicit ODE:        \f$  \dot{x} = \text{ode}(t,s,x,z,u,p)       \f$ <BR>
*  quadratures:         \f$  \dot{q} = \text{qua}(t,s,x,z,u,p)       \f$ <BR>
*  algebraic equations: \f$     0 = \text{alg}(t,s,x,z,u,p)       \f$ <BR>
*  dependent equations: \f$     y = \text{dep}(t,s,x,z,u,p)       \f$ <BR>
*  initial equations:   \f$     0 = \text{ieq}(t,s,\dot{s},x,z,u,p)  \f$ <BR>
*
*  Note that when parsed, all dynamic states, differential and algebraic, end up in the category "s" 
*  and all dynamic equations end up in the implicit category "dae". At a later state, the DAE can be
*  reformulated, for example in semi-explicit form, possibly in addition to a set of quadrature states.
*
*  The functions for reformulation is are provided as member functions to this class or as independent
*  functions located in the header file "ocp_tools.hpp".
*
*  <H3>Usage skeleton:</H3>
*  
*  1. Call constructor with an FMI conformant XML file or pass an empty string ("") to build up the OCP from scratch <BR>
*  > FlatOCP ocp(xml_file_name)
*  
*  2. Set options <BR>
*  > ocp.setOption(...,...)
*
*  3. Initialize and parse XML <BR>
*  > ocp.init()
*
*  4. Modify/add variables, equations, optimization <BR>
*  > ...
*  
*  When the optimal control problem is in a suitable form, it is possible to either generate functions
*  for numeric/symbolic evaluation or exporting the OCP formulation into a new FMI conformant XML file.
*  The latter functionality is not yet available.
*
*  \date 2011
*  \author Joel Andersson
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
    
    /// Implicitly defined states
    std::vector<Variable>& s();
    
    /// Differential states
    std::vector<Variable>& x();
    
    /// Algebraic states
    std::vector<Variable>& z();
    
    /// Quadrature states
    std::vector<Variable>& q();
    
    /// Dependent variables
    std::vector<Variable>& y();
    
    /// Independent parameters
    std::vector<Variable>& p();
    
    /// Control signals
    std::vector<Variable>& u();
    //@}
    
    #if 0
    /** @name Equations
    *  Get all equations of a particular type 
    */
    //@{
    /// Fully implicit DAE (length == s().size())
    std::vector<SX>& dae();
    
    /// Explicit ODE  (length == x().size())
    std::vector<SX>& ode();
    
    /// Algebraic equations (length == z().size())
    std::vector<SX>& alg();
    
    /// Quadrature states (length == q().size())
    std::vector<SX>& qua();
    
    /// Dependent equations (length == y().size())
    std::vector<SX>& dep();
    
    /// Initial equations (remove?)
    std::vector<SX>& ieq();
    //@}
    #endif

    /** @name Optimization
    *  Formulate the dynamic optimization problem. Note that the variable bounds are located inside
       the respective variable.
    */
    //@{
    /// Mayer terms in the objective
    std::vector<SX>& mterm();
    
    /// Lagrange terms in the objective
    std::vector<SX>& lterm();
    
    /// Interval start time
    double& t0();
    
    /// Interval final time
    double& tf();
    
    /// Interval start time is free
    bool& t0_free();
    
    /// Interval final time is free
    bool& tf_free();
    //@}
};

#ifdef SWIG
%extend FlatOCP {
  // Not inherited
  std::string __repr__()  { return $self->getRepresentation(); }
}
#endif // SWIG

} // namespace OptimalControl
} // namespace CasADi

#endif //FLAT_OCP_HPP
