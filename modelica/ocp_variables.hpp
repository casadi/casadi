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

#ifndef OCP_VARIABLES_HPP
#define OCP_VARIABLES_HPP

#include "casadi/printable_object.hpp"
#include "variable.hpp"

namespace CasADi{
  namespace Modelica{

/** Symbolic, object oriented representation of an optimal control problem (OCP) */
class OCPVariables : public PrintableObject{
  public:    
    /// Constructor (automatic type conversion allowed)
    OCPVariables(const Variable& var);

#ifndef SWIG    
    /// Print a destription of the object
    virtual void print(std::ostream &stream=std::cout) const;
#endif // SWIG

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

#ifdef SWIG
%extend OCPVariables {
  // Print (why is this not inherited?)
  std::string __repr__()  { return $self->getRepresentation(); }
}
#endif

  } // namespace Modelica
} // namespace CasADi

#endif // OCP_VARIABLES_HPP


