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

#ifndef FMI_PARSER_HPP
#define FMI_PARSER_HPP

#include "optimica_ocp.hpp"
#include "xml_node.hpp"

namespace CasADi{
namespace Modelica{

// Forward declaration
class FMIParserInternal;
  
class FMIParser : public SharedObject{
  public:
    /// Default (empty) constructor
    FMIParser();
    
    /// Create an FMI parser instance given the filename
    FMIParser(const std::string& filename);

    /// Parse from XML to C++ format
    OCP& parse();

    /// Get the OCP
    OCP& ocp();

    /// Get the OCP (const ref)
    const OCP& ocp() const;

    /// Access to the internal class
    FMIParserInternal* operator->();

    /// Const access to the internal class
    const FMIParserInternal* operator->() const;

    /// Check if the node is pointing to the right type of object
    virtual bool checkNode() const;
};

#ifdef SWIG
%extend FMIParser {
  // Not inherited
  std::string __repr__()  { return $self->getRepresentation(); }
}
#endif // SWIG

} // namespace Modelica
} // namespace CasADi

#endif //FMI_PARSER_HPP
