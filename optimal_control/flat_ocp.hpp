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

#include "symbolic_ocp.hpp"
#include "xml_node.hpp"
#include "variable.hpp"

namespace CasADi{
namespace OptimalControl{

  /// Tree structure for storing variables
  class VariableTree{
    public:
      /// Access a sub-collection by name
      VariableTree& subByName(const std::string& name, bool allocate=false);

      /// Access a sub-collection by index
      VariableTree& subByIndex(int ind, bool allocate=false);
      
      /// Get all variables
      void getAll(std::vector<Variable>& v) const;
  
      /// Get all names
      std::vector<std::string> getNames() const;
      
      /// Print node
      #ifndef SWIG
      void print(std::ostream &stream, int indent=0) const;
      #endif // SWIG

      /// Variable
      Variable var_;
      
      /// Children nodes
      std::vector<VariableTree> children_;
      
      /// Names of children
      std::map<std::string,int> name_part_;
              
  };

// Forward declaration
class FlatOCPInternal;
  
class FlatOCP : public OptionsFunctionality{
  public:
    /// Default (empty) constructor
    FlatOCP();
    
    /// Create an FMI parser instance given the filename
    FlatOCP(const std::string& filename);

    /// Parse from XML to C++ format
    void parse();

    /// Get the OCP
    SymbolicOCP& ocp();

    /// Get the OCP (const ref)
    const SymbolicOCP& ocp() const;

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
