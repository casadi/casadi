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
  Variables:
  t:     time
  x:     differential states
  z:     algebraic states
  y:     dependent variables
  p:     independent parameters
  u:     control signals
  xdot:  state derivatives
  
  Equations
  fully implicit DAE:       0 == dae(t,x,xdot,z,u,p)
  explicit ODE:          xdot == ode(t,x,z,u,p)
  quadratures:           qdot == qua(t,x,z,u,p)
  algebraic equations:      0 == alg(t,x,z,u,p)
  dependent equations:      y == dep(t,x,z,u,p)
  initial equations:        0 == ieq(t,x,xdot,z,u,p)


  Usage skeleton (starting with an XML file):
  
  ** 1. Call constructor (pass an empty string ("") to start with an empty file
  > FlatOCP ocp(xml_file_name)
  
  ** 2. Set options
  > ocp.setOption(...,...)

  ** 3. Initialize and parse XML
  > ocp.init()

  ** 4. Modify/add variables, equations, optimization
  > ...
  
  ** 5. Export FMI XML (not implemented)
  > ocp.exportFMI()

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
